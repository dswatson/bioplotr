#' Transform Counts
#'
#' This convenient wrapper function converts raw counts (or count-like data) to the
#' log2-counts per million scale, following RLE normalization.
#'
#' @importFrom edgeR DGEList calcNormFactors cpm
#'

lcpm <- function(mat) {
  mat <- DGEList(mat)
  mat <- calcNormFactors(mat, method = 'RLE')
  mat <- cpm(mat, log = TRUE, prior.count = 1L)
}
