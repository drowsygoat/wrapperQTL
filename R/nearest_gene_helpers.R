#' Nearest gene (by TSS) from a GTF for given chromosome/position(s)
#' (patched to safely harmonize seqlevels)
#'
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges GRanges seqnames start end strand distanceToNearest
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols queryHits subjectHits
#' @importFrom GenomeInfoDb renameSeqlevels seqlevels seqlevelsStyle keepSeqlevels
#' @export
nearest_gene_from_gtf <- function(gtf_file, chr, pos, output = c("gene","distance","both"), tss_gr = NULL) {
  output <- match.arg(output)
  stopifnot(length(chr) == length(pos))

  ## 1) Build/accept TSS (width=1)
  if (is.null(tss_gr)) {
    stopifnot(is.character(gtf_file), length(gtf_file) == 1L, file.exists(gtf_file))
    gr <- rtracklayer::import(gtf_file)
    if (!"type" %in% names(S4Vectors::mcols(gr))) stop("GTF import lacks a 'type' column.")
    gr_gene <- gr[S4Vectors::mcols(gr)$type == "gene"]
    if (length(gr_gene) == 0L) stop("No 'gene' features found in the GTF.")
    gene_name <- S4Vectors::mcols(gr_gene)$gene_name
    if (is.null(gene_name)) gene_name <- S4Vectors::mcols(gr_gene)$gene_id
    if (is.null(gene_name)) gene_name <- rep(NA_character_, length(gr_gene))
    tss_pos <- ifelse(as.character(GenomicRanges::strand(gr_gene)) == "+",
                      GenomicRanges::start(gr_gene),
                      GenomicRanges::end(gr_gene))
    tss_gr <- GenomicRanges::GRanges(
      seqnames = GenomicRanges::seqnames(gr_gene),
      ranges   = IRanges::IRanges(start = tss_pos, width = 1),
      strand   = GenomicRanges::strand(gr_gene)
    )
    S4Vectors::mcols(tss_gr)$gene_name <- as.character(gene_name)
    if ("gene_id" %in% names(S4Vectors::mcols(gr_gene)))
      S4Vectors::mcols(tss_gr)$gene_id <- as.character(S4Vectors::mcols(gr_gene)$gene_id)
  } else {
    stopifnot(inherits(tss_gr, "GRanges"))
  }

  ## 2) Queries -> GRanges points
  q <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges   = IRanges::IRanges(start = as.integer(pos), width = 1),
    strand   = "*"
  )

  ## 3) Harmonize chromosome style SAFELY (use renameSeqlevels)
  tss_chr_levels <- as.character(GenomeInfoDb::seqlevels(tss_gr))
  q_chr_levels   <- as.character(GenomeInfoDb::seqlevels(q))

  tss_has_chr <- any(startsWith(tss_chr_levels, "chr"))
  q_has_chr   <- any(startsWith(q_chr_levels,   "chr"))

  if (tss_has_chr && !q_has_chr) {
    # map "1"->"chr1", "2"->"chr2", ...
    map <- setNames(paste0("chr", q_chr_levels), q_chr_levels)
    q <- GenomeInfoDb::renameSeqlevels(q, map)
  } else if (!tss_has_chr && q_has_chr) {
    # map "chr1"->"1", ...
    map <- setNames(sub("^chr", "", q_chr_levels), q_chr_levels)
    q <- GenomeInfoDb::renameSeqlevels(q, map)
  }
  # After renaming, keep only shared levels
  common <- intersect(GenomeInfoDb::seqlevels(q), GenomeInfoDb::seqlevels(tss_gr))
  if (!length(common)) {
    if (output == "gene")     return(rep(NA_character_, length(chr)))
    if (output == "distance") return(rep(NA_integer_,  length(chr)))
    return(data.frame(
      query_chr          = chr,
      query_pos          = as.integer(pos),
      gene_name          = NA_character_,
      gene_strand        = NA_character_,
      tss_pos            = NA_integer_,
      distance_bp        = NA_integer_,
      signed_distance_bp = NA_integer_
    ))
  }
  q       <- GenomeInfoDb::keepSeqlevels(q,       common, pruning.mode = "coarse")
  tss_use <- GenomeInfoDb::keepSeqlevels(tss_gr,  common, pruning.mode = "coarse")

  ## 4) Nearest by TSS
  hits <- GenomicRanges::distanceToNearest(q, tss_use, ignore.strand = TRUE)
  q_idx <- S4Vectors::queryHits(hits)
  s_idx <- S4Vectors::subjectHits(hits)
  dist  <- S4Vectors::mcols(hits)$distance

  # Prepare outputs aligned to original queries
  n <- length(chr)
  out_gene   <- rep(NA_character_, n)
  out_dist   <- rep(NA_integer_,  n)
  out_sdist  <- rep(NA_integer_,  n)
  out_strand <- rep(NA_character_, n)
  out_tsspos <- rep(NA_integer_,  n)

  # Map back by position in original vector
  orig_idx <- seq_len(n)            # q preserves order of input vectors
  gene_strand <- as.character(GenomicRanges::strand(tss_use))[s_idx]
  tss_pos_sel <- GenomicRanges::start(tss_use)[s_idx]
  q_pos_sel   <- GenomicRanges::start(q)[q_idx]
  signed      <- ifelse(gene_strand == "+", q_pos_sel - tss_pos_sel, tss_pos_sel - q_pos_sel)

  out_gene[orig_idx[q_idx]]   <- S4Vectors::mcols(tss_use)$gene_name[s_idx]
  if ("gene_id" %in% names(S4Vectors::mcols(tss_use))) {
    miss <- is.na(out_gene)
    out_gene[miss] <- S4Vectors::mcols(tss_use)$gene_id[s_idx][match(which(miss), orig_idx[q_idx][miss])]
  }
  out_dist[orig_idx[q_idx]]   <- as.integer(dist)
  out_sdist[orig_idx[q_idx]]  <- as.integer(signed)
  out_strand[orig_idx[q_idx]] <- gene_strand
  out_tsspos[orig_idx[q_idx]] <- as.integer(tss_pos_sel)

  if (output == "gene")     return(out_gene)
  if (output == "distance") return(out_dist)

  data.frame(
    query_chr          = chr,
    query_pos          = as.integer(pos),
    gene_name          = out_gene,
    gene_strand        = out_strand,
    tss_pos            = out_tsspos,
    distance_bp        = out_dist,
    signed_distance_bp = out_sdist,
    stringsAsFactors = FALSE
  )
}


#' Build (and reuse) a TSS GRanges from a GTF (helper for speed)
#'
#' @param gtf_file Path to GTF.
#' @return A \code{GRanges} of width-1 TSS with \code{gene_name} (and possibly \code{gene_id})
#'   in \code{mcols}. Pass this as \code{tss_gr=} to \code{nearest_gene_from_gtf()} for repeated queries.
#' @export
nearest_gene_from_gtf_build_tss <- function(gtf_file) {
  stopifnot(is.character(gtf_file), length(gtf_file) == 1L, file.exists(gtf_file))
  gr <- rtracklayer::import(gtf_file)
  if (!"type" %in% names(S4Vectors::mcols(gr))) {
    stop("Imported GTF lacks a 'type' column; cannot locate gene features.")
  }
  gr_gene <- gr[S4Vectors::mcols(gr)$type == "gene"]
  if (length(gr_gene) == 0L) stop("No 'gene' features found in the GTF.")

  gene_name <- S4Vectors::mcols(gr_gene)$gene_name
  if (is.null(gene_name)) gene_name <- S4Vectors::mcols(gr_gene)$gene_id
  if (is.null(gene_name)) gene_name <- rep(NA_character_, length(gr_gene))

  tss_pos <- ifelse(as.character(GenomicRanges::strand(gr_gene)) == "+",
                    GenomicRanges::start(gr_gene),
                    GenomicRanges::end(gr_gene))
  tss_gr <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(gr_gene),
    ranges   = IRanges::IRanges(start = tss_pos, width = 1),
    strand   = GenomicRanges::strand(gr_gene)
  )
  S4Vectors::mcols(tss_gr)$gene_name <- as.character(gene_name)
  if ("gene_id" %in% names(S4Vectors::mcols(gr_gene)))
    S4Vectors::mcols(tss_gr)$gene_id <- as.character(S4Vectors::mcols(gr_gene)$gene_id)
  tss_gr
}
