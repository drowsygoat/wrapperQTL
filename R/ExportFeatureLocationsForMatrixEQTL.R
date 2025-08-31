#' Export gene feature locations for Matrix eQTL analysis
#'
#' This function extracts and formats genomic feature location data from a GTF file for a specified set of genes,
#' and writes it in a format compatible with Matrix eQTL. It supports gene identifiers from Seurat objects,
#' `SummarizedExperiment`s, or character vectors.
#'
#' @param gtf_file Path to a GTF file containing gene annotations.
#' @param genes A Seurat object, a `SummarizedExperiment` object, or a character vector of gene identifiers.
#' @param output_dir Directory where the output file will be saved. Defaults to the current directory.
#' @param use_gene_name Logical. If `TRUE`, match genes by `gene_name` first, falling back to `gene_id`. 
#' If `FALSE`, match only by `gene_id`.
#' @param add_chr_prefix Logical. If `TRUE`, adds `"chr"` prefix to chromosome names (e.g., `1` -> `chr1`).
#' @param filename Name of the output file. Defaults to `"feature_locations.txt"`.
#'
#' @return (Invisibly) the path to the output file containing the feature locations.
#'
#' @details The function filters GTF entries to keep only those of type `"gene"`, and retrieves
#' their chromosome, start, and end coordinates. If gene identifiers are duplicated between 
#' `gene_name` and `gene_id`, matches by `gene_name` are prioritized.
#'
#' @importFrom rtracklayer import
#' @importFrom dplyr select filter mutate distinct bind_rows
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @seealso \code{\link{convert_chromosome_ids}} for chromosome name formatting
#'
#' @examples
#' \dontrun{
#' ExportFeatureLocationsForMatrixEQTL(
#'   gtf_file = "annotation.gtf",
#'   genes = c("GeneA", "GeneB", "GeneC"),
#'   output_dir = "eqtl_files",
#'   use_gene_name = TRUE,
#'   add_chr_prefix = TRUE
#' )
#' }
#'
#' @export
ExportFeatureLocationsForMatrixEQTL <- function(
  gtf_file,
  genes,
  output_dir = ".",
  use_gene_name = TRUE,
  add_chr_prefix = FALSE,
  filename = "feature_locations.txt"
) {
  # Load required libraries
  require(rtracklayer)
  require(dplyr)
  require(SummarizedExperiment)
  
  # Extract gene list depending on input type
  if (inherits(genes, "Seurat")) {
    gene_list <- rownames(genes)
  } else if (inherits(genes, "SummarizedExperiment")) {
    gene_list <- rownames(genes)
  } else if (is.character(genes)) {
    gene_list <- genes
  } else {
    stop("`genes` must be a Seurat object, SummarizedExperiment, or a character vector of gene names.")
  }
  
  # Import GTF
  gtf <- rtracklayer::import(gtf_file)
  
  # Filter for gene entries
  genes_gtf <- gtf[gtf$type == "gene"]
  gene_df <- as.data.frame(genes_gtf)
  
  # Ensure required columns
  if (!"gene_id" %in% colnames(gene_df)) stop("`gene_id` column not found in GTF file.")
  if (!"gene_name" %in% colnames(gene_df)) gene_df$gene_name <- NA_character_
  
  # Prepare base feature table

  if (add_chr_prefix) {
    feature_df <- gene_df %>%
      dplyr::select(gene_id, gene_name, seqnames, start, end) %>%
      dplyr::distinct() %>%
      dplyr::mutate(chr = convert_chromosome_ids(seqnames, direction = "add"))
  } else {
    feature_df <- gene_df %>%
      dplyr::select(gene_id, gene_name, seqnames, start, end) %>%
      dplyr::distinct() %>%
      dplyr::mutate(chr = as.character(seqnames))
  }
  
  # Match genes using both gene_name and gene_id
  if (use_gene_name) {
    matched_gene_name <- feature_df %>%
      dplyr::filter(gene_name %in% gene_list) %>%
      dplyr::mutate(geneid = gene_name)
    
    matched_gene_id <- feature_df %>%
      dplyr::filter(gene_id %in% gene_list) %>%
      dplyr::mutate(geneid = gene_id)
    
    # Combine, preferring gene_name matches if duplicated
    feature_loc_filtered <- bind_rows(matched_gene_name, matched_gene_id) %>%
      dplyr::distinct(geneid, .keep_all = TRUE)
  } else {
    # Use only gene_id for matching
    feature_loc_filtered <- feature_df %>%
      dplyr::filter(gene_id %in% gene_list) %>%
      dplyr::mutate(geneid = gene_id)
  }

  # Final formatting
  feature_loc_filtered <- feature_loc_filtered %>%
    dplyr::mutate(s1 = start, s2 = end) %>%
    dplyr::select(geneid, chr, s1, s2)

  # Ensure output directory exists
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Write to file
  output_path <- file.path(output_dir, filename)
  write.table(
    feature_loc_filtered,
    file = output_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  
  message("✅ Feature location file saved to: ", output_path)
  return(invisible(output_path))
}
