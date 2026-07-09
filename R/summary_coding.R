###########################################################
#      SPARTA - Step 5: Summarize gene-centric CODING results
###########################################################
# Aggregates the per-job coding results from Step 4.1 into genome-wide tables,
# significant-mask CSVs, and Manhattan + Q-Q plots. One summary is produced per
# phenotype. Run from the repository root:  Rscript R/summary_coding.R
###########################################################

rm(list = ls())
gc()

library(STAAR)
library(STAARpipeline)
library(STAARpipelineSummary)
library(parallel)
library(data.table)

## SPARTA coding summary function (reads the single-cohort STAAR-O / STAAR-B columns)
source("R/Gene_Centric_Coding_Results_Summary_anno.R")

###########################################################
#           User Input
###########################################################
## Phenotype id(s) to summarize (one output folder per phenotype).
field_all <- c("26526_i2")

## Root holding the Step 4.1 output. For phenotype <field>, the per-job result
## files (Coding_<kk>.Rdata, kk = 1..gene_centric_coding_jobs_num) live in
##   <results_root>/<field>/
## and the summary is written to
##   <results_root>/<field>/Coding_summary/
results_root <- "./step3_output/"

## Number of gene-centric jobs genome-wide (chromosome x 100-gene blocks).
gene_centric_coding_jobs_num <- 195
## Prefix of the per-job result files, i.e. <gene_centric_results_name>_<kk>.Rdata
## (matches output_file_name in SPARTA_coding.R).
gene_centric_results_name <- "Coding"
## Genome-wide significance threshold and minimum number of variants per mask.
alpha <- 5e-07
SNV_cutoff <- 30
## Parallel workers across phenotypes.
n_cores <- 1

###########################################################
#           Main
###########################################################
get_res <- function(field) {
  input_path  <- paste0(results_root, field, "/")
  output_path <- paste0(results_root, field, "/Coding_summary/")
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  Gene_Centric_Coding_Results_Summary2(
    gene_centric_coding_jobs_num = gene_centric_coding_jobs_num,
    input_path = input_path, output_path = output_path,
    gene_centric_results_name = gene_centric_results_name,
    alpha = alpha, manhattan_plot = TRUE, QQ_plot = TRUE,
    SNV_cutoff = SNV_cutoff)
}
invisible(mclapply(field_all, get_res, mc.cores = n_cores))

## (optional) combine the significant coding masks across phenotypes into one table
if (length(field_all) > 1) {
  res_all <- NULL
  for (field in field_all) {
    f <- paste0(results_root, field, "/Coding_summary/coding_sig.csv")
    if (file.exists(f)) {
      res <- fread(f)
      if (nrow(res) > 0) { res$field <- field; res_all <- rbind(res_all, res, fill = TRUE) }
    }
  }
  if (!is.null(res_all))
    write.csv(res_all, paste0(results_root, "coding_sig_all.csv"),
              row.names = FALSE, quote = FALSE)
}
