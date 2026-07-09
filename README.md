# SPARTA (UK Biobank RAP / DNANexus)
A powerful and computationally efficient framework for rare variant association testing that leverages machine-learning-predicted phenotypes to boost power in large-scale biobank studies.

This is the **UK Biobank Research Analysis Platform (DNANexus)** version of SPARTA: the gene-centric steps run as DNANexus array jobs, pulling inputs from your RAP project with `dx download` and pushing results back with `dx upload`. For the local-filesystem version, see [`SPARTA`](https://github.com/shuangsong0110/SPARTA). See that README for the method overview and the predicted-phenotype format, which are identical here.


## Table of contents

* [Installation & Authentication](#package-step-1-installation--authentication)
* [Prepare the Predicted Phenotype](#dna-step-2-prepare-the-predicted-phenotype)
* [Fit the Null Model](#scroll-step-3-fit-the-null-model)
* [Gene-centric Rare Variant Tests on the RAP](#rocket-step-4-gene-centric-rare-variant-tests-on-the-rap)
* [Summary and Figures](#framed_picture-step-5-summary-and-figures)


## :package: Step 1: Installation & Authentication

Inside the RAP worker, the gene-centric scripts install the DNANexus toolkit and authenticate at run time. **Never hard-code your token in the script** — export it as an environment variable, and the script reads it with `Sys.getenv("DX_TOKEN")`:

```
export DX_TOKEN=xxxxxxxxxxxxxxxxxxxxxxxx     # your DNANexus API token
```

The STAAR ecosystem must also be available on the worker. Use the STAARpipeline Docker image (all required packages, including `STAARpipelinePheWAS`, are pre-installed):
```
docker pull zilinli/staarpipeline:0.9.9
```
Run each job inside this image, e.g. `docker run --entrypoint Rscript zilinli/staarpipeline:0.9.9 --vanilla SPARTA_coding.R ${arrayid}`. See the [`SPARTA`](https://github.com/shuangsong0110/SPARTA#package-step-1-installation) README for the full R dependency list if you build your own image.


## :dna: Step 2: Prepare the predicted phenotype

Identical to the local version — a `data.frame` with columns `ID`, `pheno` (observed; `NA` if unobserved), and `s_pheno` (predicted for all samples). See the [`SPARTA`](https://github.com/shuangsong0110/SPARTA#dna-step-2-prepare-the-predicted-phenotype) README.


## :scroll: Step 3: Fit the null model

The null model is fit off-platform (e.g. on your local cluster) exactly as in the local version — see <a href="R/SPARTA_null_model.R">**SPARTA_null_model.R**</a> and the [`SPARTA`](https://github.com/shuangsong0110/SPARTA#scroll-step-3-fit-the-null-model) README (Steps 3.1 sparse GRM and 3.2 null model). Upload the resulting per-phenotype `.rda` files (each storing `obj_nullmodel`, `resid_ratio`, `pi1`) to your RAP project so the gene-centric jobs can `dx download` them.


## :rocket: Step 4: Gene-centric rare variant tests on the RAP

Both scripts partition the genome into **195 jobs** (blocks of 100 genes per chromosome); the array index `arrayid` (1..195) selects one block. Each job authenticates with `$DX_TOKEN`, downloads its inputs (aGDS, null models, annotation catalog; the noncoding jobs also download the per-gene noncoding variant-id files) from the RAP project, runs the SPARTA test, and uploads the per-job `.Rdata` result. Hardy-Weinberg filtering is assumed to be handled upstream in variant QC. Edit the `dx` project paths in the `User Input` block to match your project.

### Step 4.1 Gene-centric coding regions

#### Script: <a href="R/SPARTA_coding.R">**SPARTA_coding.R**</a>
SPARTA gene-centric coding analysis across five functional categories (putative loss of function; missense; disruptive missense; putative loss of function + disruptive missense; synonymous).
#### Output: 195 `Coding_<kk>.Rdata` files uploaded to your RAP project.

### Step 4.2 Gene-centric non-coding regions

#### Script: <a href="R/SPARTA_noncoding.R">**SPARTA_noncoding.R**</a>
SPARTA gene-centric noncoding analysis across seven regulatory categories (promoter-CAGE; promoter-DHS; enhancer-CAGE; enhancer-DHS; UTR; upstream; downstream).
#### Output: 195 `Noncoding_<kk>.Rdata` files uploaded to your RAP project.

**Launching the array on the RAP.** Submit one job per `arrayid` (1..195), passing it as the command-line argument, e.g. with the Swiss Army Knife app or a custom applet:
```
for arrayid in $(seq 1 195); do
  dx run swiss-army-knife \
    -iin="R/SPARTA_coding.R" \
    -icmd="Rscript SPARTA_coding.R ${arrayid}" \
    --name "SPARTA_coding_${arrayid}" -y
done
```
(Set `DX_TOKEN` in the applet environment, or authenticate the worker as your project requires.)


## :framed_picture: Step 5: Summary and figures

#### Script: <a href="R/summary_coding.R">**summary_coding.R**</a> and <a href="R/summary_noncoding.R">**summary_noncoding.R**</a>
After downloading the per-job results from the RAP project, summarize them exactly as in the local version — see the [`SPARTA`](https://github.com/shuangsong0110/SPARTA#framed_picture-step-5-summary-and-figures) README.
**Output:** a Manhattan plot, a Q-Q plot, and CSVs of the significant coding / noncoding masks.


## Data Availability
The whole-genome functional annotation data and the precomputed annotation principal components are available at the [Functional Annotation of Variant - Online Resource (FAVOR)](https://favor.genohub.org) site and [FAVOR Essential Database](https://doi.org/10.7910/DVN/1VGTJI).


## Citation
If you use SPARTA, please cite:

Shuang Song, Ruoyu Wang, Haoyu Yang, Max Elliott, Xihao Li, Jingnan Du, Randy Lee Buckner, Xihong Lin. *Augmented synthetic phenotype analysis enhances rare variant association discovery in biobank-scale whole genome sequencing studies.*


## :busts_in_silhouette: Maintainer
Please contact Shuang Song (shuangsong@hsph.harvard.edu) if there are any problems or questions.


## References
[STAAR](https://github.com/xihaoli/STAAR/)
[STAARpipeline](https://github.com/xihaoli/STAARpipeline/)
[STAARpipelineSummary](https://github.com/xihaoli/STAARpipelineSummary/)
[FastSparseGRM](https://github.com/rounakdey/FastSparseGRM/)
