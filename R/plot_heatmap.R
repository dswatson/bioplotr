#' Omic Heatmap
#'
#' This function visualizes a probe by sample matrix as a heatmap.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples.
#' @param group Optional character or factor vector of length equal to sample
#'   size. Alternatively, a data frame or list of such vectors, optionally
#'   named.
#' @param covar Optional continuous covariate of length equal to sample size.
#'   Alternatively, a data frame or list of such vectors, optionally named.
#' @param dist Distance measure to be used. Supports all methods available in
#'   \code{\link[stats]{dist}}, \code{Rfast::\link[Rfast]{Dist}}, and \code{
#'   \link[vegan]{vegdist}}, as well as those implemented in the \code{bioDist} 
#'   package. See Details.
#' @param p Power of the Minkowski distance.
#' @param top Optional number (if > 1) or proportion (if < 1) of top probes to
#'   be used for t-SNE.
#' @param filter_method String specifying whether to apply a \code{"pairwise"}
#'   or \code{"common"} filter if \code{top} is non-\code{NULL}.
#' @param center Center each probe prior to computing distances?
#' @param hclustfun The agglomeration method to be used for hierarchical
#'   clustering. Supports all methods available in \code{\link[stats]{hclust}}.
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
#'   well as all sequential color schemes available in \code{RColorBrewer}. 
#'   Alternatively, a character vector of colors representing a smooth gradient, 
#'   or a list of such vectors with length equal to the number of continuous 
#'   variables to visualize.
#' @param pal_tiles String specifying the color palette to use for heatmap
#'   tiles. Options include the complete collection of \code{\href{
#'   https://bit.ly/2n7D6tF}{viridis}} palettes, as well all diverging color 
#'   schemes available in \code{\href{https://bit.ly/2ipuEjn}{RColorBrewer}}. 
#'   Alternatively, a character vector of at least two colors.
#' @param title Optional plot title.
#'
#' @details
#' Heatmaps are a common and intuitive way to display the values of an omic data
#' matrix, especially after top probes have been selected for closer
#' investigation. Hierarchical clustering dendrograms cluster both the rows and
#' the columns, revealing latent structure in the data. Annotation tracks atop
#' the figure may be used to investigate associations with phenotypic features.
#' 
#' Available distance measures include: \code{"euclidean"}, \code{"maximum"},
#' \code{"manhattan"}, \code{"canberra"}, \code{"minkowski"}, \code{"cosine"},
#' \code{"pearson"}, \code{"kendall"}, \code{"spearman"}, \code{"bray"}, \code{
#' "kulczynski"}, \code{"jaccard"}, \code{"gower"}, \code{"altGower"}, \code{
#' "morisita"}, \code{"horn"}, \code{"mountford"}, \code{"raup"}, \code{
#' "binomial"}, \code{"chao"}, \code{"cao"}, \code{"mahalanobis"}, \code{"MI"},
#' or \code{"KLD"}. Some distance measures are unsuitable for certain types of
#' data. See \code{\link{dist_mat}} for more details on these methods and links
#' to documentation on each.
#'
#' The \code{top} argument optionally filters data using either probewise
#' variance (if \code{filter_method = "common"}) or the leading fold change
#' method of Smyth et al. (if \code{filter_method = "pairwise"}). See \code{
#' \link{plot_mds}} for more details.
#'
#' @examples
#' mat <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
#' grp <- rep(c("A", "B"), each = 5)
#' plot_heatmap(mat, group = grp)
#'
#' @export
#' @importFrom RColorBrewer brewer.pal
#'

plot_heatmap <- function(dat,
                         group = NULL,
                         covar = NULL,
                          dist = 'pearson',
                             p = 2L,
                           top = NULL,
                 filter_method = 'pairwise',
                        center = FALSE,
                     hclustfun = 'average',
                     pal_group = 'npg',
                     pal_covar = 'Blues',
                     pal_tiles = 'RdBu',
                         title = 'Omic Heatmap') {

  # Preliminaries
  if (!group %>% is.null) {
    group <-  dat %>% format_features(group, var_type = 'Categorical')
    grp_clrs <- group %>% track_cols(pal_group, var_type = 'Categorical')
  } else {
    grp_clrs <- NULL
  }
  if (!covar %>% is.null) {
    covar <- format_features(dat, covar, var_type = 'Continuous')
    cov_clrs <- track_cols(covar, pal_covar, var_type = 'Continuous')
  } else {
    cov_clrs <- NULL
  }
  if (!c(group, covar) %>% is.null) {
    anno <- c(group, covar)
    ann_clrs <- c(grp_clrs, cov_clrs)
  }
  d <- c('euclidean', 'maximum', 'manhattan', 'canberra', 'minkowski',
         'bhattacharyya', 'hellinger', 'kullback_leibler', 'cosine', 
         'bray', 'kulczynski', 'jaccard', 'gower', 'altGower', 'morisita', 
         'horn', 'mountford', 'raup' , 'binomial', 'chao', 'cao',
         'mahalanobis', 'pearson', 'kendall', 'spearman', 'MI')
  if (!dist %in% d) {
    stop('dist must be one of ', stringify(d, 'or'), '.')
  }
  if (!filter_method %in% c('pairwise', 'common')) {
    stop('filter_method must be either "pairwise" or "common".')
  }
  hclusts <- c('ward.D', 'ward.D2', 'single', 'complete', 'average',
               'mcquitty', 'median', 'centroid')
  if (!hclustfun %in% hclusts) {
    stop('hclustfun must be one of ', stringify(hclusts, 'or'), '. ',
         'See ?hclust.')
  }
  pal_tiles <- colorize(pal_tiles, var_type = 'Continuous')

  # Tidy data
  dat <- matrixize(dat)
  dm <- dist_mat(dat, dist, p, top, filter_method, center) %>% as.dist(.)

  # Plot
  suppressPackageStartupMessages(require(NMF))
  if (anno %>% is.null) {
    aheatmap(dat, distfun = dm, hclustfun = hclustfun, scale = 'row',
             col = pal_tiles, main = title, border_color = 'grey60')
  } else {
    aheatmap(dat, distfun = dm, hclustfun = hclustfun, scale = 'row',
             annCol = anno, annColors = ann_clrs, col = pal_tiles,
             main = title, border_color = 'grey60')
  }

}


### Check out d3heatmap?
### The distfun = dm thing is broken. Row dendrogram can't deal.
### Some good way to fix this?
### ann_clrs should be a list

