#' MDS Plot
#'
#' This function plots a low-dimensional projection of an omic data matrix using
#' multi-dimensional scaling.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples. It is strongly recommended that data be
#'   filtered and normalized prior to plotting. Raw counts stored in \code{
#'   \link[edgeR]{DGEList}} or \code{\link[DESeq2]{DESeqDataSet}} objects are
#'   automatically extracted and transformed to the log2-CPM scale, with a
#'   warning.
#' @param group Optional character or factor vector of length equal to sample
#'   size, or up to two such vectors organized into a list or data frame. Supply
#'   legend title(s) by passing a named list or data frame.
#' @param covar Optional continuous covariate. If non-\code{NULL}, then plot can
#'   render at most one \code{group} variable. Supply legend title by passing
#'   a named list or data frame.
#' @param metric Logical. Perform classical (i.e. metric) MDS or nonmetric MDS?
#'   See Details.
#' @param dist Distance measure to be used. Supports all methods available in
#'   \code{\link[stats]{dist}} and \code{\link[vegan]{vegdist}}, as well as
#'   those implemented in the \code{bioDist} package. See Details.
#' @param p Power of the Minkowski distance.
#' @param top Optional number (if > 1) or proportion (if < 1) of top probes to
#'   be used for MDS.
#' @param filter_method String specifying whether to apply a \code{"pairwise"}
#'   or \code{"common"} filter if \code{top} is non-\code{NULL}. See Details.
#' @param pcs Vector specifying which principal coordinates to plot. Must be of
#'   length two unless \code{D3 = TRUE}.
#' @param label Label data points by sample name? Defaults to \code{FALSE}
#'   unless \code{group} and \code{covar} are both \code{NULL}. If \code{TRUE},
#'   then plot can render at most one phenotypic feature.
#' @param pal_group String specifying the color palette to use if \code{group}
#'   is non-\code{NULL}. Options include \code{"ggplot"}, as well as the
#'   complete collection of \code{
#'   \href{https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html}{
#'   ggsci}} palettes, which can be identified by name (e.g., \code{"npg"},
#'   \code{"aaas"}, etc.). Alternatively, any character vector of colors with
#'   length equal to the number of levels in \code{group}.
#' @param pal_covar String specifying the color palette to use if \code{covar}
#'   is non-\code{NULL}. Options include \code{"blues"}, \code{"greens"}, \code{
#'   "purples"}, \code{"greys"}, \code{"oranges"}, and \code{"reds"}.
#'   Alternatively, any character vector of colors representing a smooth
#'   gradient.
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"right"}, \code{
#'   "left"}, \code{"top"}, \code{"bottom"}, \code{"topright"}, \code{
#'   "topleft"}, \code{"bottomright"}, or \code{"bottomleft"}.
#' @param hover Show sample name by hovering mouse over data point? If \code{
#'   TRUE}, the plot is rendered in HTML and will either open in your browser's
#'   graphic display or appear in the RStudio viewer.
#' @param D3 Render plot in three dimensions?
#'
#' @details
#' This function plots the samples of an omic data matrix in a two- or
#' three-dimensional principal coordinate subspace. MDS is an easy and popular
#' method for unsupervised cluster detection. It can also aid in spotting
#' potential outliers, and generally helps to visualize the latent structure of
#' a data set.
#'
#' Classical MDS is implemented by the \code{\link[stats]{cmdscale}} function,
#' which finds the optimal two-dimensional projection of a distance matrix by
#' minimizing the strain of the coordinate mapping (Torgerson, 1958).
#' Nonmetric MDS (NMDS) is implemented by the \code{\link[vegan]{monoMDS}}
#' function, which uses isotonic regression to find the monotonic transformation
#' that minimizes the stress of the embedding (Kruskal, 1964).
#'
#' MDS requires a distance matrix as input. Internal functions median centre the
#' probe vectors and compute samplewise distances using one of the following
#' methods: \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"}, \code{
#' "canberra"}, \code{"minkowski"}, \code{"cosine"}, \code{"pearson"}, \code{
#' "kendall"}, \code{"spearman"}, \code{"bray"}, \code{"kulczynski"}, \code{
#' "jaccard"}, \code{"gower"}, \code{"altGower"}, \code{"morisita"}, \code{
#' "horn"}, \code{"mountford"}, \code{"raup"}, \code{"binomial"}, \code{"chao"},
#' \code{"cao"}, \code{"mahalanobis"}, \code{"MI"}, or \code{"KLD"}. Some
#' distance measures are unsuitable for certain types of data. See \code{
#' \link{dist_mat}} for more details on these methods and links to full
#' documentation on each.
#'
#' If \code{top} is non-\code{NULL}, then data can either be filtered by
#' probewise variance (\code{filter_method = "common"}) or using the leading
#' fold change method of Smyth et al. (\code{filter_method = "pairwise"}). In
#' the latter case, pairwise distances are calculated using only the \code{top}
#' most differentially expressed probes between the two samples. This method is
#' appropriate when different molecular pathways are relevant for distinguishing
#' different pairs of samples. To run MDS on the complete data, set \code{
#' top = NULL}. This is functionally equivalent to running PCA on the full
#' matrix. See \code{\link{plot_pca}}.
#'
#' @references
#' Cox, T.F. & Cox, M.A.A. (2001). \emph{Multidimensional Scaling}. Second
#' edition. Chapman and Hall.
#'
#' Kruskal, J.B. (1964).
#' \href{https://pdfs.semanticscholar.org/e041/08dc293c9cd7cabf32ee1524eaab0d4641b3.pdf}{
#' Multidimensional scaling by optimizing goodness of fit to a nonmetric
#' hypothesis}. \emph{Psychometrika}, \emph{29}(1): 1-27.
#'
#' Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., & Smyth, G.K.
#' (2015).
#' \href{https://www.ncbi.nlm.nih.gov/pubmed/25605792}{limma powers differential
#' expression analyses for RNA-sequencing and microarray studies}. \emph{Nucleic
#' Acids Res.}, \emph{43}(7): e47.
#'
#' Torgerson, W.S. (1958). \emph{Theory and Methods of Scaling}. New York:
#' Wiley.
#'
#' @examples
#' mat <- matrix(rnorm(1000 * 5), nrow = 1000, ncol = 5)
#' plot_mds(mat)
#'
#' library(DESeq2)
#' dds <- makeExampleDESeqDataSet()
#' dds <- rlog(dds)
#' plot_mds(dds, group = colData(dds)$condition)
#'
#' @seealso
#' \code{\link[limma]{plotMDS}}, \code{\link{plot_pca}}
#'
#' @export
#' @importFrom vegan metaMDS
#' @import dplyr
#' @import ggplot2
#'

plot_mds <- function(dat,
                     group = NULL,
                     covar = NULL,
                      dist = 'euclidean',
                         p = 2,
                    metric = TRUE,
                       top = 500,
             filter_method = 'pairwise',
                       pcs = c(1, 2),
                     label = FALSE,
                 pal_group = 'npg',
                 pal_covar = 'blues',
                     title = NULL,
                    legend = 'right',
                     hover = FALSE,
                        D3 = FALSE) {

  # Preliminaries
  if (ncol(dat) < 3L) {
    stop(paste('dat includes only', ncol(dat), 'samples; ',
               'need at least 3 for MDS.'))
  }
  if (!(group %>% is.null)) {
    group <- dat %>% format_features(group, var_type = 'Categorical')
    if (length(group) > 2L) {
      stop('Plot can render at most two categorical features.')
    }
    if (length(group) == 2L && !(covar %>% is.null)) {
      stop('Plot can render at most one categorical feature when a continuous ',
           'covariate is also supplied.')
    }
    group_cols <- colorize(pal = pal_group, var_type = 'Categorical',
                           n = length(levels(group[[1L]])))
  }
  if (!(covar %>% is.null)) {
    covar <- dat %>% format_features(covar, var_type = 'Continuous')
    if (length(covar) != 1L) {
      stop('Plot can render at most one continuous feature.')
    }
    covar_cols <- colorize(pal = pal_covar, var_type = 'Continuous')
  }
  if (!(c(group, covar) %>% is.null)) {
    features <- c(covar, group)
    feature_names <- names(features)
    names(features) <- paste0('Feature', seq_along(features))
  } else {
    features <- NULL
  }
  d <- c('euclidean', 'maximum', 'manhattan', 'canberra', 'minkowski',
         'cosine', 'bray', 'kulczynski', 'jaccard', 'gower', 'altGower',
         'morisita', 'horn', 'mountford', 'raup' , 'binomial', 'chao', 'cao',
         'mahalanobis', 'pearson', 'kendall', 'spearman', 'MI', 'KLD')
  if (!dist %in% d) {
    stop('dist must be one of ', stringify(d, 'or'), '.')
  }
  if (!filter_method %in% c('pairwise', 'common')) {
    stop('filter_method must be either "pairwise" or "common".')
  }
  if (length(pcs) > 2L & !D3) {
    stop('pcs must be of length 2 when D3 = FALSE.')
  } else if (length(pcs) > 3L) {
    stop('pcs must be a vector of length <= 3.')
  }
  if (label && length(features) == 2L) {
    stop('If label is TRUE, then plot can render at most one phenotypic ',
         'feature.')
  }
  if (title %>% is.null) {
    title <- 'MDS'
  }
  loc <- c('right', 'left', 'top', 'bottom',
           'topright', 'topleft', 'bottomright', 'bottomleft')
  if (!legend %in% loc) {
    stop('legend must be one of ', stringify(loc, 'or'), '.')
  }

  # Tidy data
  dat <- matrixize(dat)
  dm <- dist_mat(dat, dist, p, top, filter_method) %>% as.dist(.)
  if (metric) {                                  # MDS
    mds <- suppressWarnings(cmdscale(dm, k = max(pcs)))
  } else {
    converged <- FALSE
    max_iter <- 20L
    while(!converged) {
      mds <- suppressWarnings(
        metaMDS(dm, k = max(pcs), trace = 0L, trymax = max_iter,
                autotransform = FALSE)
      )
      if (mds$converged) {
        converged <- TRUE
      } else {
        max_iter <- max_iter + 20L
      }
    }
    mds <- mds$points
  }
  df <- data_frame(Sample = colnames(dat))       # Melt
  if (length(pcs) == 2L) {
    df <- df %>% mutate(PC1 = mds[, min(pcs)],
                        PC2 = mds[, max(pcs)])
  } else {
    other <- setdiff(pcs, c(min(pcs), max(pcs)))
    df <- df %>% mutate(PC1 = mds[, min(pcs)],
                        PC2 = mds[, other],
                        PC3 = mds[, max(pcs)])
  }
  if (!(features %>% is.null)) {
    df <- df %>% cbind(tbl_df(features))
  }

  # Build plot
  if (top %>% is.null) {
    if (metric) {
      xlab <- paste0('MDS', min(pcs))
      ylab <- paste0('MDS', max(pcs))
    } else {
      xlab <- paste0('NMDS', min(pcs))
      ylab <- paste0('NMDS', max(pcs))
    }
  } else {
    if (metric) {
      xlab <- paste0('Leading logFC, MDS', min(pcs))
      ylab <- paste0('Leading logFC, MDS', max(pcs))
    } else {
      xlab <- paste0('Leading logFC, NMDS', min(pcs))
      ylab <- paste0('Leading logFC, NMDS', max(pcs))
    }
  }
  embed(df, group, covar, group_cols, covar_cols, feature_names,
        label, title, xlab, ylab, legend, hover, D3)

}

# Interactive options:
# 1) filter probes
# 2) filter samples
# 3) change PCs

