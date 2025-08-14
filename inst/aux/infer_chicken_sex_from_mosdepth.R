#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
mosdepth_dir <- args[1]  # e.g. "path/to/mosdepth"

# Check directory
if (!dir.exists(mosdepth_dir)) {
  stop("Directory does not exist: ", mosdepth_dir)
}

# Find all summary files
summary_files <- list.files(mosdepth_dir, pattern = "summary\\.txt$", recursive = TRUE, full.names = TRUE)

# Collect coverage values
results <- data.frame(Sample = character(), Z = numeric(), W = numeric(), ZW_ratio = numeric())

for (file in summary_files) {
  sample_name <- basename(dirname(file))
  df <- read.table(file, header = TRUE)
  
  z_cov <- df$mean[df$chrom == "Z"]
  w_cov <- df$mean[df$chrom == "W"]
  
  if (length(z_cov) == 1 && length(w_cov) == 1 && w_cov > 0) {
    results <- rbind(results, data.frame(
      Sample = sample_name,
      Z = z_cov,
      W = w_cov,
      ZW_ratio = z_cov / w_cov
    ))
  }
}

# ✅ Corrected: ZW → female, ZZ → male
results$Sex <- ifelse(results$ZW_ratio > 1.5, "female", "male")
results$Sex_numeric <- ifelse(results$Sex == "female", 1, 0)  # For Matrix eQTL: female = 1, male = 0

# Save Z:W barplot
pdf("ZW_coverage_ratio.pdf", width = 8, height = 5)

ggplot(results, aes(x = Sample, y = ZW_ratio, fill = Sex)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "red") +
  labs(title = "Z:W Chromosome Coverage Ratio (female = ZW)",
       x = "Sample",
       y = "Z / W Mean Coverage Ratio") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("female" = "orange", "male" = "steelblue"))

dev.off()

# Save Matrix eQTL covariates
covariate_df <- data.frame(
  row.names = "sex",
  t(results$Sex_numeric)
)
colnames(covariate_df) <- results$Sample

write.table(covariate_df, file = "covariates_sex.txt", sep = "\t", quote = FALSE, col.names = NA)

cat("✅ Done. Output files:\n - ZW_coverage_ratio.pdf\n - covariates_sex.txt\n")
