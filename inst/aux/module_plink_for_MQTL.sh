#!/bin/bash
## PLINK2 → Matrix eQTL prep

# This script converts a VCF into PLINK2 binaries, computes allele frequencies & PCs, and formats files for **Matrix eQTL**.  
# Runs PLINK2 from a **Docker container** converted to **Singularity** (`bioinfo_toolkit` image via `sing.sh`).

# ### Steps
# 1. **VCF → PLINK2 binary** (`--make-bed`, `--max-alleles 2`, `--set-all-var-ids @:#`)
# 2. **Allele freq**: `plink2 --freq` → `genotype_freq.afreq`
# 3. **SNP matrix**: `plink2 --recode A-transpose` → `.traw` → R reformats:
#    - `GeneID` = `SNP_ALT`
#    - Genotypes 0/1/2, missing = `NA`
#    - Sample IDs cleaned & made unique
# 4. **PC covariates**: `plink2 --pca 10` → eigenvec → transposed for Matrix eQTL

# ### Outputs (`${JOB_NAME}/matrix_eqtl_input/`)
# - `SNP.txt` – SNP genotypes (rows: SNPs, cols: samples)
# - `Covariates.txt` – PCs (rows: PC1..PC10, cols: samples)

# ### Notes
# - Missing genotypes remain `NA` (no dosage imputation).
# - Adjust `--chr-set` if using non-mouse genome.
# - Requires: `plink2`, R (`data.table`), Singularity image from Docker (`bioinfo_toolkit`)

# replace sing.sh with docker_sing.sh !!!!
# 
# MODULE PARAMETERS
RUN_COMMAND="run_shell_command.sh"
JOB_NAME="plink2_COVs_GENs_with_missing"
PARTITION="shared"
NODES=1
TIME="100"
TASKS=1
CPUS=40
DRY="with_eval"

VCF="/cfs/klemming/projects/snic/sllstore2017078/kaczma-workingdir/sarek/AIL_haplotypecaller/run1_with_haplotypecaller/variant_calling/haplotypecaller/joint_variant_calling/joint_germline.vcf.gz"

process_file() {

    local vcf=$1
    local outdir=$2

    export vcf outdir

    $RUN_COMMAND -J "${JOB_NAME}" -p "$PARTITION" -n "$TASKS" -t "$TIME" -N "$NODES" -c "$CPUS" -d "$DRY" -- \
    '
        mkdir -p "${outdir}"
        mkdir -p "${outdir}/matrix_eqtl_input"

        CYAN="\033[0;36m"
        NC="\033[0m" # No Color

        mkdir -p "${outdir}"
        mkdir -p "${outdir}/matrix_eqtl_input"

        echo -e "${CYAN}🔧 PART 1: Convert VCF to PLINK binary${NC}"
        sing.sh -B "$(dirname "${vcf}")" bioinfo_toolkit plink2 \
        --vcf "${vcf}" \
        --chr-set 50 no-xy \
        --make-bed \
        --max-alleles 2 \
        --set-all-var-ids @:# \
        --out "${outdir}/genotype" \
        --threads 16

        echo -e "${CYAN}📊 PART 2: Calculate allele frequencies${NC}"
        sing.sh -B "$(dirname "${vcf}")" bioinfo_toolkit plink2 \
        --bfile "${outdir}/genotype" \
        --chr-set 50 no-xy \
        --freq \
        --out "${outdir}/genotype_freq"

        echo -e "${CYAN}📁 PART 3: Recode to raw format for Matrix eQTL${NC}"
        sing.sh -B "$(dirname "${vcf}")" bioinfo_toolkit plink2 \
        --bfile "${outdir}/genotype" \
        --chr-set 50 no-xy \
        --recode A-transpose \
        --out "${outdir}/matrix_eqtl_input/genotype_matrix"

        echo -e "${CYAN}🧬 PART 4: Run PCA for genotype covariates${NC}"
        sing.sh -B "$(dirname "${vcf}")" bioinfo_toolkit plink2 \
        --bfile "${outdir}/genotype" \
        --chr-set 50 no-xy \
        --read-freq "${outdir}/genotype_freq.afreq" \
        --pca 10 \
        --out "${outdir}/matrix_eqtl_input/genotype_pcs"

        echo -e "${CYAN}🧾 PART 5: Format PLINK .raw and eigenvec for Matrix eQTL${NC}"
        Rscript - <<EOF
        
library(data.table)

# File paths
raw_file <- "${outdir}/matrix_eqtl_input/genotype_matrix.traw"
pcs_file <- "${outdir}/matrix_eqtl_input/genotype_pcs.eigenvec"
snp_out  <- "${outdir}/matrix_eqtl_input/SNP.txt"
cov_out  <- "${outdir}/matrix_eqtl_input/Covariates.txt"

message("=== [1] START: Loading genotype matrix from: ", raw_file)
message("File size: ", file.info(raw_file)\$size / 1024^2, " MB")

# Load the transposed raw matrix
geno <- fread(raw_file, data.table = FALSE)

# Generate GeneID: SNP_POSITION_ALTALLELE
gene_ids <- paste0(geno\$SNP, "_", geno[[6]])  # 6th column is ALT

# Drop metadata columns (first 6), keep genotype values
genotypes <- geno[, -(1:6)]

# Clean sample IDs
clean_sample_ids <- gsub("^0_", "", colnames(genotypes))
clean_sample_ids <- make.unique(clean_sample_ids, sep = "_")
colnames(genotypes) <- clean_sample_ids

# Final genotype matrix: GeneID + genotype values
genotypes <- cbind(GeneID = gene_ids, genotypes)

# Check for duplicates
if (anyDuplicated(colnames(genotypes)) > 0) {
stop("❌ Duplicate column names found in SNP matrix!")
}

# Write SNP matrix
message("=== [2] Writing SNP matrix to: ", snp_out)
fwrite(genotypes, file = snp_out, sep = "\t", quote = FALSE, na = "NA")
message("✅ SNP matrix saved.")

# ========== Covariates (e.g., PCs) ==========
message("=== [3] Loading PCA covariates from: ", pcs_file)
pcs <- read.table(pcs_file, header = FALSE)

ids <- gsub("^0_", "", pcs[, 2])  # Clean sample IDs (match genotype cols)
pcs_data <- pcs[, -(1:2)]
colnames(pcs_data) <- paste0("PC", 1:ncol(pcs_data))
rownames(pcs_data) <- ids

# Transpose: Matrix eQTL wants samples as columns
pcs_t <- t(pcs_data)

# Write covariates
message("=== [4] Writing covariates to: ", cov_out)
write.table(pcs_t, file = cov_out, sep = "\t", quote = FALSE, col.names = NA)
message("✅ Covariates saved.")

message("🏁 Done.")
EOF

        echo "✅ Matrix eQTL input files created:"
        echo "   - SNP matrix: ${outdir}/matrix_eqtl_input/SNP.txt"
        echo "   - Covariates: ${outdir}/matrix_eqtl_input/Covariates.txt"
    '
}

process_file $VCF $JOB_NAME
    
# without A-transpose missing genotypes are replaced in with probable dosage