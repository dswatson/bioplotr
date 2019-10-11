#' Similarity Matrix Heatmap
#'
#' This function displays pairwise distances between samples as a heatmap.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples. It is strongly recommended that data be
#'   filtered and normalized prior to plotting. Raw counts stored in \code{
#'   \link[edgeR]{DGEList}} or \code{\link[DESeq2]{DESeqDataSet}} objects are
#'   automatically extracted and transformed to the log2-CPM scale, with a
#'   warning. Alternatively, an object of class \code{dist}.
#' @param group Optional character or factor vector of length equal to sample
#'   size. Alternatively, a data frame or list of such vectors, optionally
#'   named.
#' @param covar Optional continuous covariate of length equal to sample size.
#'   Alternatively, a data frame or list of such vectors, optionally named.
#' @param dist Distance measure to be used. Supports all methods available in
#'   \code{\link[stats]{dist}} and \code{\link[vegan]{vegdist}}, as well as
#'   those implemented in the \code{bioDist} package. See Details.
#' @param p Power of the Minkowski distance.
#' @param top Optional number (if > 1) or proportion (if < 1) of top probes to
#'   be used for distance calculations.
#' @param filter_method String specifying whether to apply a \code{"pairwise"}
#'   or \code{"common"} filter if \code{top} is non-\code{NULL}. See Details.
#' @param hclustfun The agglomeration method to be used for hierarchical
#'   clustering. Supports any method available in \code{\link[stats]{hclust}}.
#' @param pal_group String specifying the color palette to use if \code{group}
#'   is non-\code{NULL}, or a vector of such strings with length equal to the
#'   number of vectors passed to \code{group}. Options include \code{"ggplot"}, 
#'   all qualitative color schemes available in \code{\href{
#'   https://bit.ly/2ipuEjn}{RColorBrewer}}, and the complete collection of 
#'   \code{\href{http://bit.ly/2bxnuGB}{ggsci}} palettes. Alternatively, a 
#'   character vector of colors with length equal to the cumulative number of 
#'   levels in \code{group}.
#' @param pal_covar String specifying the color palette to use if \code{covar}
#'   is non-\code{NULL}, or a vector of such strings with length equal to the
#'   number of vectors passed to \code{covar}. Options include the complete 
#'   collection of \code{\href{https://bit.ly/2n7D6tF}{viridis}} palettes, as 
#'   well as all sequential color schemes available in \code{\href{
#'   https://bit.ly/2ipuEjn}{RColorBrewer}}. Alternatively, a character vector 
#'   of colors representing a smooth gradient, or a list of such vectors with 
#'   length equal to the number of continuous variables to visualize.
#' @param pal_tiles String specifying the color palette to use for heatmap
#'   tiles. Options include the complete collection of \code{\href{
#'   https://bit.ly/2n7D6tF}{viridis}} palettes, as well all diverging color 
#'   schemes available in \code{\href{https://bit.ly/2ipuEjn}{RColorBrewer}}. 
#'   Alternatively, a character vector of at least two colors.
#' @param title Optional plot title.
#'
#' @details
#' Similarity matrices are a valuable tool for exploratory data analysis. A
#' hierarchical clustering dendrogram atop the figure helps identify potential
#' clusters and/or outliers in the data. Annotation tracks can help investigate
#' associations with phenotypic features.
#'
#' Different distance metrics and agglomeration methods can lead to different
#' results. The default settings, which perform average linkage hierarchical
#' clustering on a Euclidean distance matrix, are mathematically straightforward
#' and commonly used for omic EDA. Complete linkage is also fairly common.
#'
#' Pearson distance, defined as 1 - the Pearson correlation, is another popular
#' method for evaluating sample similarity. Mutual information and
#' Kullback-Leibler divergence are more complicated distance metrics that
#' require some simplifying assumptions to be efficiently applied to continuous
#' data distributions. See \code{\link[bioDist]{MIdist}} and \code{
#' \link[bioDist]{KLdist.matrix}} for more info.
#'
#' The \code{top} argument optionally filters data using either probewise
#' variance (if \code{filter_method = "common"}) or the leading fold change
#' method of Smyth et al. (if \code{filter_method = "pairwise"}). See \code{
#' \link{plot_mds}} for more details.
#'
#' @examples
#' mat <- matrix(rnorm(5000), nrow = 1000, ncol = 5)
#' plot_similarity(mat, title = "Nothin' Doin'")
#'
#' library(DESeq2)
#' dds <- makeExampleDESeqDataSet()
#' dds <- rlog(dds)
#' plot_similarity(dds, group = colData(dds)$condition,
#'                 title = "Somethin' Cookin'")
#'
#' @export
#' @importFrom RColorBrewer brewer.pal
#'

plot_similarity <- function(dat,
                            group = NULL,
                            covar = NULL,
                             dist = 'euclidean',
                                p = 2L,
                              top = NULL,
                    filter_method = 'pairwise',
                        hclustfun = 'average',
                        pal_group = 'npg',
                        pal_covar = 'Blues',
                        pal_tiles = 'RdBu',
                            title = 'Sample Similarity Matrix') {

  # Preliminaries
  if (!dat %>% is('dist')) {
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
  }
  if (!group %>% is.null) {
    group <- dat %>% format_features(group, var_type = 'Categorical')
    grp_cols <- group %>% track_cols(pal_group, var_type = 'Categorical')
  } else {
    grp_cols <- NULL
  }
  if (!covar %>% is.null) {
    covar <- dat %>% format_features(covar, var_type = 'Continuous')
    cov_cols <- covar %>% track_cols(pal_covar, var_type = 'Continuous')
  } else {
    cov_cols <- NULL
  }
  if (!c(group, covar) %>% is.null) {
    anno <- c(group, covar)
    ann_cols <- c(grp_cols, cov_cols)
  } else {
    anno <- NULL
  }
  hclusts <- c('ward.D', 'ward.D2', 'single', 'complete', 'average',
               'mcquitty', 'median', 'centroid')
  if (!hclustfun %in% hclusts) {
    stop('hclustfun must be one of ', stringify(hclusts, 'or'), '. ',
         'See ?hclust.')
  }
  pal_cols <- colorize(pal_tiles, var_type = 'Continuous')
  if (pal_tiles == 'RdBu') {
    pal_cols <- rev(pal_cols)
  }

  # Tidy data
  if (!dat %>% is('dist')) {
    dat <- matrixize(dat)
    dm <- dist_mat(dat, dist, p, top, filter_method)
  } else {
    dm <- dat
  }

  # Build plot
  suppressPackageStartupMessages(require(NMF))
  if (anno %>% is.null) {
    aheatmap(dm, col = pal_cols, Rowv = FALSE, revC = TRUE, main = title,
             distfun = function(x) as.dist(x), hclustfun = hclustfun,
             border_color = 'grey60')
  } else {
    aheatmap(dm, col = pal_cols, Rowv = FALSE, revC = TRUE, main = title,
             distfun = function(x) as.dist(x), hclustfun = hclustfun,
             annCol = anno, annColors = ann_cols, border_color = 'grey60')
  }

}


### REPLACE AHEATMAP WITH PHEATMAP ###
# pal_covar and pal_group need to adapt to detect number of features 





