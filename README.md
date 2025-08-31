# wrapperQTL

`wrapperQTL` provides an end-to-end workflow for preparing, running, and post-processing QTL analyses.  
It standardizes SNP preparation, expression/covariate matrix generation, Matrix eQTL execution, and downstream operations.

---

## 📦 Installation

You can install the latest version of `wrapperQTL` directly from GitHub using the `devtools` package:

```r
# Install devtools if you don't have it
install.packages("devtools")

# Install wrapperQTL from GitHub
devtools::install_github("drowsygoat/wrapperQTL")
```

## Installation & Setup

Scripts to be run from the CLI are installed as system files.  

Add the package’s `inst/aux` directory to your `PATH`:

```bash
export PATH="$(Rscript -e "cat(normalizePath(system.file('aux', package='wrapperQTL')))"):$PATH"
```

> This step is necessary even if using Docker/Singularity, as SLURM cannot be called from within a container.

---

## Typical Workflow

### 1. Prepare SNP Data from VCF

Convert combined VCF files into PLINK and Matrix eQTL–compatible formats with:

```bash
module_plink_for_MQTL.sh
```

This script performs:

1. **VCF → PLINK binary conversion** (`.bed/.bim/.fam`)  
2. **Allele frequency calculation**  
3. **Recode to transposed raw format** (Matrix eQTL input)  
4. **Genotype PCA** (for covariates)  
5. **Final formatting** into:
   - `SNP.txt` – genotype matrix  
   - `Covariates.txt` – genotype covariates (PCs)

---

### 2. Clean and Organize Inputs

- **Standardize sample IDs:**

  ```r
  CleanIDs(".../SNP.txt")
  ```

  Sample IDs should be in the form `IDXX`, where `XX` is a sample number with leading zeros.  
  This function attempts to correct IDs automatically.

- **(Optional) Extract sample and cluster metadata.**  
  In typical cases, column names are in the format `cluster_x_sample`:

  ```r
  colData(se_atac_full)$cluster   <- gsub("_x_.*$", "", colnames(se_atac_full))
  colData(se_atac_full)$sample_id <- gsub("^.*_x_", "", colnames(se_atac_full))
  ```

- **(Optional) Inspect the "depth" covariate** using:

  ```bash
  PlotClusterCovariateSummary.R
  ```

---

### 3. Infer Sex Covariate from MosDepth Data

Generate inputs suitable for Matrix eQTL with:

```r
infer_sex_from_mosdepth()
```

This produces a `covariate_sex.txt` file.  
Optionally, fix IDs with:

```r
CleanIDs(".../covariate_sex.txt")
```

---

### 4. Prepare Matrix eQTL Input Files

Generate inputs suitable for Matrix eQTL with:

```r
prepareMatrixEQTLInputs()
```

This step:

- Filters expression, SNP, and covariate matrices  
- Accepts `SummarizedExperiment` or raw matrices  
- Optionally annotates features via GTF  
- Splits large datasets into chunks  
- Saves processed files to disk  

See function documentation for details.

---

### 5. Run Matrix eQTL

Execute the core analysis with:

```bash
module_R_MEQTL.sh
```

- Calls `matrixEQTLrun.R` internally  
- Produces results in a dedicated directory per group/cluster:  

  ```
  group_{clusterID}_results
  ```

For many clusters (e.g., C1 through C49):

```bash
for i in {1..49}; do module_R_MEQTL.sh group_C${i}_results; done
```

---

### 6. Summarize Test Counts

Collect the number of tests per run (used for FDR calculation):

```bash
gather_ntests.R
```

Example:

```bash
find . -type d -name "*MEQTL_res*" -print0 | xargs -0 -I{} gather_ntests.R {}
```

---

### 7. Post-Processing & Annotation

Finalize and annotate results with:

```r
run_qtl_workflow()
```

This performs:

- FDR recalculation (optional)  
- Peak range detection (cis/trans)  
- QTL type summaries  
- (Optional) Annotation with:
  - Feature metadata (`rowData`)  
  - Nearest gene by TSS (via GTF)  

**Outputs include:**

- Filtered QTL tables  
- Cis/trans peak ranges  
- Annotated results (nearest genes if enabled)  
- CSV and PDF outputs (if `write_outputs = TRUE`)  

---

## Workflow Summary

1. **Prepare SNPs** (VCF → PLINK → Matrix eQTL input)  
2. **Clean IDs & metadata**  
3. **Prepare expression/SNP/covariate matrices** with `prepareMatrixEQTLInputs()`  
4. **Run Matrix eQTL** via `module_R_MEQTL.sh`  
5. **Summarize test counts** with `gather_ntests.R`  
6. **Post-process** results with `run_qtl_workflow()`  

---

## Documentation

- Each function and script has detailed documentation (see source files).  
- This README provides only a **general workflow overview**.
