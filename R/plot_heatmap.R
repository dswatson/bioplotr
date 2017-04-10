#' Omic Heatmap
#'
#' This function visualizes a probe by sample matrix as a heatmap.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples.
#' @param anno Optional character, factor, numeric, or logical vector of length
#'   equal to sample size. Alternatively, a data frame or list of such vectors,
#'   optionally named. Values are used to color one or several annotation tracks
#'   atop the heatmap.
#' @param dist Distance measure to be used. Currently supports any method available
#'   in \code{\link[stats]{dist}} or \code{\link[stats]{cor}}.
#' @param hclustfun The agglomeration method to be used for hierarchical clustering.
#'   See \code{\link[stats]{hclust}} for available options.
#' @param col Color palette to use for heatmap tiles. Preset options include \code{
#'   "RdBu"} for red to blue gradient, \code{"GrRd"} for green to red gradient, and
#'   \code{"BuYl"} for blue to yellow gradient. Alternatively, any user-supplied
#'   color palette is acceptable.
#' @param title Optional plot title.
#'
#' @details
#' Heatmaps are a common and intuitive way to display the values of an omic data
#' matrix, especially after top probes have been selected for closer investigation.
#' Hierarchical clustering dendrograms cluster both the rows and the columns,
#' revealing latent structure in the data. Annotation tracks atop the figure may be
#' used to investigate associations with phenotypic features.
#'
#' @examples
#' mat <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
#' grp <- rep(c("A", "B"), each = 5)
#' plot_heatmap(mat, anno = grp)
#'
#' @export
#' @importFrom purrr map_lgl
#' @importFrom edgeR calcNormFactors cpm
#' @importFrom DESeq2 sizeFactors normalizationFactors estimateSizeFactors
#' @importFrom SummarizedExperiment assay
#' @importFrom limma getEAWP
#' @importFrom NMF aheatmap
#' @import RColorBrewer
#'

plot_heatmap <- function(dat,
                         anno = NULL,
                         dist = 'pearson',
                    hclustfun = 'average',
                          col = 'RdBu',
                        title = NULL) {

  # Preliminaries
  if (is.data.frame(anno)) {
    anno <- as.list(anno)
  } else if (!is.list(anno)) {
    anno <- list('Variable' = anno)
  } else {
    if (is.null(names(anno))) {
      if (length(anno) == 1L) {
        names(anno) <- 'Variable'
      } else {
        names(anno) <- paste('Variable', seq_along(anno))
      }
    }
  }
  if (any(map_lgl(seq_along(anno), function(j) {
    length(anno[[j]]) != ncol(dat)
  }))) {
    stop('anno length must match number of samples in dat.')
  }
  if (any(map_lgl(seq_along(anno), function(j) {
    if (is.numeric(anno[[j]])) var(anno[[j]]) == 0L
    else length(unique(anno[[j]])) == 1L
  }))) {
    stop('anno is invariant.')
  }
  if (!dist %in% c('euclidean', 'pearson', 'MI', 'KLD')) {
    stop('dist must be one of "euclidean", "pearson", "MI", or "KLD".')
  }
  if (!hclustfun %in% c('ward.D', 'ward.D2', 'single', 'complete', 'average',
                        'mcquitty', 'median', 'centroid')) {
    stop('hclustfun must be one of "ward.D", "ward.D2", "single", "complete", ',
         '"average", "mcquitty", "median", or "centroid". See ?hclust.')
  }
  if (col == 'RdBu') {
    col <- colorRampPalette(brewer.pal(10L, 'RdBu'))(n = 256L)
  } else if (col == 'GrRd') {
    col <- colorRampPalette(c('green', 'black', 'red'))(n = 256L)
  } else if (col == 'BuYl') {
    col <- colorRampPalette(c('blue', 'grey', 'yellow'))(n = 256L)
  }
  if (is.null(title)) {
    title <- 'Omic Heatmap'
  }

  # Tidy data
  if (is(dat, 'DGEList')) {
    keep <- rowSums(dat$counts) > 0L             # Minimal count filter
    dat <- dat[keep, ]
    if (is.null(dat$samples$norm.factors) |      # Calculate size factors
        all(dat$samples$norm.factors == 1L)) {
      dat <- calcNormFactors(dat)
    }
    dat <- cpm(dat, log = TRUE, prior.count = 1L)
    warning('Transforming raw counts to log2-CPM scale.')
  } else if (is(dat, 'DESeqDataSet')) {
    if (is.null(sizeFactors(dat)) & is.null(normalizationFactors(dat))) {
      dat <- estimateSizeFactors(dat)            # Normalize counts
    }
    dat <- counts(dat, normalized = TRUE)
    keep <- rowMeans(dat) > 0L                   # Minimal count filter
    dat <- dat[keep, , drop = FALSE]
    dat <- cpm(dat, log = TRUE, prior.count = 1L)
    warning('Transforming raw counts to log2-CPM scale.')
  } else if (is(dat, 'DESeqTransform')) {
    dat <- assay(dat)
  } else {
    dat <- getEAWP(dat)$expr
    keep <- rowSums(is.finite(dat)) == ncol(dat)
    dat <- dat[keep, , drop = FALSE]
  }

  # Plot
  if (is.null(anno)) {
    aheatmap(dat, distfun = dist, scale = 'row', col = col,
             hclustfun = hclustfun, main = title)
  } else {
    aheatmap(dat, distfun = dist, scale = 'row', col = col,
             hclustfun = hclustfun, main = title, annCol = anno)
  }

}


# Replace aheatmap with pheatmap?
# Rows should be centered, but not scaled, I think?

