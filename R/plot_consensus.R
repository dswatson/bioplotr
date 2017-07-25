#' Consensus Matrix Heatmap
#'
#' This function displays pairwise consensus cluster distances as a heatmap.
#'
#' @param cc A list created by a call to \code{
#'   \link[ConsensusClusterPlus]{ConsensusClusterPlus}}.
#' @param k Integer specifying number of clusters to visualize.
#' @param group Optional character or factor vector of length equal to sample
#'   size. Alternatively, a data frame or list of such vectors, optionally
#'   named.
#' @param covar Optional continuous covariate of length equal to sample size.
#'   Alternatively, a data frame or list of such vectors, optionally named.
#' @param hclustfun The agglomeration method to be used for hierarchical
#'   clustering. Supports any method available in \code{\link[stats]{hclust}}.
#' @param pal_clust String specifying the color palette to use for cluster
#'   assignments. Options include \code{"ggplot"}, all qualitative color schemes
#'   available in \code{RColorBrewer}, and the complete collection of \code{
#'   \href{http://bit.ly/2bxnuGB}{ggsci}} palettes.
#' @param pal_group String specifying the color palette to use if \code{group}
#'   is non-\code{NULL}, or a vector of such strings with length equal to the
#'   number of vectors passed to \code{group}. Options are the same as those for
#'   \code{pal_clust}. Alternatively, a character vector of colors with length
#'   equal to the cumulative number of levels in \code{group}.
#' @param pal_covar String specifying the color palette to use if \code{covar}
#'   is non-\code{NULL}, or a vector of such strings with length equal to the
#'   number of vectors passed to \code{covar}. Options include all sequential
#'   color schemes available in \code{RColorBrewer}. Alternatively, a
#'   character vector of colors representing a smooth gradient, or a list of
#'   such vectors with length equal to the number of continuous variables to
#'   visualize.
#' @param pal_tiles String specifying the color palette to use for heatmap
#'   tiles. Options include all diverging color schemes available in \code{
#'   RColorBrewer}. Alternatively, a character vector of at least two colors.
#' @param title Optional plot title.
#'
#' @details
#' Consensus clustering is a resampling procedure to evaluate cluster stability.
#' A user-specified proportion of samples are held out on each run of the
#' algorithm to test how often the remaining samples do or do not cluster
#' together. The result is a square consensus matrix for each value of cluster
#' numbers \emph{k}. Each cell of the matrix \code{mat[i, j]} represents the
#' proportion of all runs including samples \code{i} and \code{j} in which the
#' two were clustered together.
#'
#' \code{plot_consensus} converts a consensus matrix into a distance matrix by
#' taking the complement of all values and visualising the result as a heatmap.
#' These figures are similar to those presented in the original consensus
#' cluster paper by Monti et al.
#'
#' @references
#' Monti, S., Tamayo, P., Mesirov, J., & Golub, T. (2003).
#' \href{http://link.springer.com/article/10.1023/A:1023949509487}{Consensus
#' Clustering: A Resampling-Based Method for Class Discovery and Visualization
#' of Gene Expression Microarray Data}. \emph{Machine Learning}, \emph{52}:
#' 91-118.
#'
#' @examples
#' # Load libraries
#' library(ConsensusClusterPlus)
#'
#' # Import, filter, transform data
#' data(airway)
#' mat <- assay(airway)
#' mads <- apply(mat, 1, mad)
#' keep <- order(mads, decreasing = TRUE)[seq_len(1000)]
#' mat <- mat[keep, ]
#' mat <- log2((mat + 1L) / 1e6L)
#'
#' # Run consensus clustering, plot results
#' cc <- ConsensusClusterPlus(mat)
#' plot_consensus(cc, k = 2)
#'
#' @export
#' @importFrom fastcluster hclust
#' @importFrom RColorBrewer brewer.pal
#' @import dplyr
#'

plot_consensus <- function(cc,
                           k,
                           group = NULL,
                           covar = NULL,
                       hclustfun = 'average',
                       pal_clust = 'd3',
                       pal_group = 'npg',
                       pal_covar = 'Blues',
                       pal_tiles = 'RdBu',
                           title = NULL) {

  # Preliminaries
  if (cc %>% is.list) {
    if ('consensusMatrix' %in% names(cc[[k]])) {
      dat <- 1L - cc[[k]]$consensusMatrix
    } else {
      dat <- 1L - cc[[k]]
    }
  } else {
    stop('cc must be a list object containing consensus matrices.')
  }
  if (k %>% is.null || !(k %>% is.numeric) || length(k) > 1L) {
    stop('k must be an integer on [2, length(cc)].')
  }
  hclusts <- c('ward.D', 'ward.D2', 'single', 'complete', 'average',
               'mcquitty', 'median', 'centroid')
  if (!hclustfun %in% hclusts) {
    stop('hclustfun must be one of ', stringify(hclusts, 'or'), '. ',
         'See ?hclust.')
  }
  clusters <- as.dist(dat) %>%
    hclust(method = hclustfun) %>%
    cutree(k = k)
  clusters <- list('Cluster' = clusters)
  clust_cols <- clusters %>% track_cols(pal_clust, var_type = 'Categorical')
  if (!(group %>% is.null)) {
    group <- dat %>% format_features(group, var_type = 'Categorical')
    grp_cols <- group %>% track_cols(pal_group, var_type = 'Categorical')
  } else {
    grp_cols <- NULL
  }
  if (!(covar %>% is.null)) {
    covar <- dat %>% format_features(covar, var_type = 'Continuous')
    cov_cols <- covar %>% track_cols(pal_covar, var_type = 'Continuous')
  } else {
    cov_cols <- NULL
  }
  if (!(c(group, covar) %>% is.null)) {
    anno <- c(clusters, group, covar)
    ann_cols <- c(clust_cols, grp_cols, cov_cols)
  } else {
    anno <- clusters
    ann_cols <- clust_cols
  }
  pal_cols <- colorize(pal_tiles, var_type = 'Continuous')
  if (pal_tiles == 'RdBu') {
    pal_cols <- rev(pal_cols)
  }
  if (title %>% is.null) {
    title <- 'Consensus Matrix'
  }

  # Build plot
  suppressPackageStartupMessages(require(NMF))
  aheatmap(dat, col = pal_cols, Rowv = FALSE, revC = TRUE, main = title,
           distfun = function(x) as.dist(x), hclustfun = hclustfun,
           annCol = anno, annColors = ann_cols, border_color = 'grey60')

}


