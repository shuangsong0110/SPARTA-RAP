###########################################################
#      SPARTA - Step 3.2: Fit the null model
###########################################################
# For each phenotype, SPARTA combines an *observed* phenotype (available for a
# subset of samples) with an ML/AI-*predicted* phenotype (available for ALL
# samples) into a debiased "pseudo phenotype", and fits a STAAR null model on it.
#
# For each phenotype the script saves one .rda containing:
#   obj_nullmodel : STAAR null model fit on the pseudo phenotype
#   resid_ratio   : ratio of residual variance in observed (R==1) vs
#                   unobserved (R==0) samples  (used by the SPA correction)
#   pi1           : fraction of samples with an observed phenotype
###########################################################

library(data.table)
library(STAAR)
library(STAARpipeline)

###########################################################
#           User Input
###########################################################
## Phenotype field id(s) / name(s). One null model is fit per element; set a
## single value for a single trait, or several for a PheWAS.
field_all <- c("26526_i2")

## Predicted-phenotype file for a given `field` (see README Step 2). A headed
## data.frame with columns:
##   ID       sample id
##   pheno    observed phenotype (NA if unobserved)
##   s_pheno  predicted phenotype for ALL samples (from your AI/ML model)
pheno_file <- function(field) paste0("./phenotype/participant.p", field, "_norm_rf100_cali_v5.txt")

## Covariate file: a headed data.frame keyed by sample id (column `IID`),
## containing the covariates listed in `pcs` below.
covar_file <- "./phenotype/covar.txt"

## Covariates to adjust for in the null model. age2/agesex/age2sex are derived
## automatically from `age` and `sex` if they are not already columns.
pcs <- c("sex", "age", "age2", "agesex", "age2sex", paste0("pc", 1:20))

## Sequenced samples passing QC (one sample id per line, no header). Restricts
## the null model to samples that also have sequencing data, e.g.
##   samid <- seqGetData(seqOpen(gds), "sample.id")
samid_file <- "./phenotype/sequenced_sample_ids.txt"

## Sparse GRM (OPTIONAL). Set grm_file <- NULL if the samples are unrelated;
## the null model is then fit without a GRM. Generate the sparse GRM with
## FastSparseGRM (see README Step 3.1). The .RData must provide a matrix `sGRM`.
grm_file <- "/path/to/grmout_nRandomSNPs_0.sGRM.RData"   # or NULL

## Output directory
output_path <- "./step1_output/"

###########################################################
#           Main function
###########################################################
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

## load the sparse GRM once (shared across phenotypes)
if (!is.null(grm_file)) {
  load(grm_file)   # provides `sGRM`
  ## strip any suffix on GRM ids (e.g. "<id>_<id>") so they match sample ids
  grm_id <- unlist(lapply(strsplit(colnames(sGRM), "_"), function(x) x[1]))
  colnames(sGRM) <- rownames(sGRM) <- grm_id
}

## QC-passing sequenced sample ids
samid <- fread(samid_file, header = FALSE)[[1]]

fit_one <- function(field) {
  message("Fitting null model for phenotype: ", field)

  ## ---- read predicted phenotype, build observation indicator R ----
  dat <- fread(pheno_file(field))
  dat$R <- ifelse(is.na(dat$pheno), 0, 1)     # 1 = observed, 0 = unobserved
  dat$pheno[is.na(dat$pheno)] <- 0

  ## ---- attach phenotype to covariates ----
  covar <- fread(covar_file)
  covar$s_pheno <- dat$s_pheno[match(covar$IID, dat$ID)]
  covar$pheno   <- dat$pheno[match(covar$IID, dat$ID)]
  covar1 <- covar[!is.na(covar$s_pheno), ]    # keep samples with a prediction
  colnames(covar1)[which(colnames(covar1) == "IID")] <- "ID"

  ## ---- derive covariates if not already present ----
  if (!"age2"    %in% colnames(covar1)) covar1$age2    <- covar1$age^2
  if (!"agesex"  %in% colnames(covar1)) covar1$agesex  <- covar1$age  * covar1$sex
  if (!"age2sex" %in% colnames(covar1)) covar1$age2sex <- covar1$age2 * covar1$sex

  ## ---- restrict to sequenced, QC-passing (and GRM) samples ----
  covar1 <- covar1[which(covar1$ID %in% samid), ]
  if (!is.null(grm_file)) covar1 <- covar1[which(covar1$ID %in% grm_id), ]

  ## ---- SPARTA pseudo phenotype (prediction-powered debiased outcome) ----
  covar1$R <- dat$R[match(covar1$ID, dat$ID)]
  covar1$pheno_pseudo <- covar1$s_pheno +
    (1 / mean(covar1$R)) * covar1$R * (covar1$pheno - covar1$s_pheno)

  formula <- as.formula(paste("pheno_pseudo ~", paste(pcs, collapse = " + ")))
  covar1$R <- as.factor(covar1$R)

  ## ---- fit the STAAR null model ----
  if (!is.null(grm_file)) {
    covar1 <- covar1[which(covar1$ID %in% colnames(sGRM)), ]
    obj_nullmodel <- fit_nullmodel(formula, data = covar1, kins = sGRM,
                                   use_sparse = TRUE, id = "ID",
                                   family = gaussian(link = "identity"),
                                   verbose = TRUE)
  } else {
    obj_nullmodel <- fit_nullmodel(formula, data = covar1, kins = NULL,
                                   id = "ID",
                                   family = gaussian(link = "identity"),
                                   verbose = TRUE)
  }

  ## ---- SPARTA scaling factors consumed by the SPA correction ----
  ## pi1        = fraction of samples with an observed phenotype
  ## resid_ratio= residual variance in observed (R==1) vs unobserved (R==0)
  obj_nullmodel$R <- covar1$R
  tmp1 <- mean(obj_nullmodel$residuals[which(covar1$R == 1)]^2)
  tmp2 <- mean(obj_nullmodel$residuals[which(covar1$R == 0)]^2)
  resid_ratio <- tmp1 / tmp2
  pi1 <- mean(as.numeric(as.character(covar1$R)))
  message(sprintf("  resid_ratio = %.4f,  pi1 = %.4f", resid_ratio, pi1))

  ## ---- save ----
  save(obj_nullmodel, resid_ratio, pi1,
       file = paste0(output_path, "obj.STAAR_500K_", field,
                     "_norm_pseudo_rf100_v5.rda"))
}

for (field in field_all) fit_one(field)
