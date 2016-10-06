#' Visualize a similarity matrix
#'
#' @param dat Omic data matrix with rows corresponding to probes and columns
#'   to samples.
#' @param group Optional factor or character vector of length equal to sample size.
#'   Levels are used to color an annotation track atop the heatmap.
#' @param main Optional plot title.
#'
#' @details
#' This function displays the pairwise scaled Euclidean distance between samples
#' in the form of a heatmap. A hierarchical clustering dendrogram is added
#' atop the figure to help identify potential outliers and/or clusters in
#' the data.
#'
#' @examples
#' mat <- matrix(rnorm(5000), nrow = 1000, ncol = 5)
#' plot_sim_mat(mat, main = "Nothin' Doin'")
#'
#' library(DESeq2)
#' mat <- cbind(matrix(rnbinom(5000, mu = 4, size = 1), nrow = 1000, ncol = 5),
#'              matrix(rnbinom(5000, mu = 4, size = 10), nrow = 1000, ncol = 5))
#' mat <- rlog(mat)
#' grp <- rep(c("A", "B"), each = 5)
#' plot_sim_mat(mat, group = grp, main = "Somethin' Cookin'")
#'
#' @importFrom wordspace dist.matrix
#' @importFrom NMF aheatmap
#' @import RColorBrewer
#'

plot_sim_mat <- function(dat,
	                       group = NULL,
	                       main  = NULL) {

  if (is.null(main)) {
    main <- 'Sample Similarity Matrix'
  }

  dm <- dist.matrix(scale(t(dat)), method = 'euclidean')
  rb <- colorRampPalette(brewer.pal(10, 'RdBu'))(n = 256)

  if (is.null(group)) {
    aheatmap(dm, col = rb, Rowv = FALSE,
      distfun = function(x) as.dist(x), hclustfun = 'average',
      main = main)
  } else {
    aheatmap(dm, col = rb, Rowv = FALSE, annCol = list(Group = group),
      distfun = function(x) as.dist(x), hclustfun = 'average',
      main = main)
  }

}


