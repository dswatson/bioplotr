#' Plot the first two principal components of an omic data matrix
#'
#' @param dat Omic data matrix with rows corresponding to probes and columns 
#'   to samples.
#' @param group Optional factor or character vector of length equal to sample size. 
#'   Levels are used to color and shape points.
#' @param label Label data points by sample name? Defaults to FALSE unless 
#'   \code{group = NULL}.
#' @param main Optional plot title. 
#' @param legend Legend position. Must be one of \code{"outside", 
#'   "bottomleft", "bottomright", "topleft",} or \code{"topright"}.
#'
#' @details 
#' This function plots the first two principal components of an omic data 
#' matrix. The numbers printed in the x- and y-axis labels indicate the 
#' percentage of variance explained by each component. PCA is an easy and 
#' popular method for unsupervised cluster detection. It can also aid in 
#' spotting potential outliers and generally  help to visualize the latent 
#' structure of a data set.
#' 
#' @examples
#' mat <- matrix(rnorm(5000), nrow = 1000, ncol = 5)
#' plot_pca(mat)
#' 
#' library(DESeq2)
#' mat <- cbind(matrix(rnbinom(5000, mu = 4, size = 1), nrow = 1000, ncol = 5),
#'              matrix(rnbinom(5000, mu = 4, size = 10), nrow = 1000, ncol = 5))
#' mat <- rlog(mat)
#' batch <- rep(c("A", "B"), each = 5)
#' plot_pca(mat, group = batch)
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'

plot_PCA <- function(dat, 
                     group  = NULL,
                     label  = FALSE,
                     main   = NULL,
                     legend = 'outside') {

  if (is.null(group)) {
    group <- rep(1, ncol(dat))
  } else {
    if (!is.character(group) &
        !is.factor(group)) {
      stop('group must be a character or factor variable.')
    }
  }
  if (is.null(colnames(dat))) {
    sample <- paste0('Sample', seq_along(1:ncol(dat)))
  } else {
    sample <- colnames(dat)
  }
  if (is.null(main)) {
    main <- 'PCA'
  }

  pca <- prcomp(t(dat), center = TRUE, scale. = TRUE) 
  var1 <- round(pca$sdev[1]^2 / sum(pca$sdev^2) * 100, 2)
  var2 <- round(pca$sdev[2]^2 / sum(pca$sdev^2) * 100, 2)
  df <- data_frame(PC1    = pca$x[, 1], 
                   PC2    = pca$x[, 2],
                   Sample = sample,
                   Group  = group)
  
  p <- ggplot(df, aes(PC1, PC2)) + 
    geom_hline(yintercept = 0, size = 0.2) + 
    geom_vline(xintercept = 0, size = 0.2) +
    labs(title = main, 
         x = paste0('PC1 (', var1, '%)'), 
         y = paste0('PC2 (', var2, '%)')) + 
    theme_bw()

  if (label == TRUE) {
    if (!is.numeric(df$Group)) {
      p <- p + geom_text(aes(label = Sample, color = Group))
    } else {
      p <- p + geom_text(aes(label = Sample))
    }
  } else {
    if (!is.numeric(df$Group)) {
      p <- p + geom_point(aes(color = Group, shape = Group))
    } else {
      p <- p + geom_text(aes(label = Sample))
    }
  }

  if (legend == 'bottomleft') {
    p <- p + theme(legend.justification = c(0, 0), legend.position = c(0, 0))
  } else if (legend == 'bottomright') {
    p <- p + theme(legend.justification = c(1, 0), legend.position = c(1, 0))
  } else if (legend == 'topleft') {
    p <- p + theme(legend.justification = c(0, 1), legend.position = c(0, 1))
  } else if (legend == 'topright') {
    p <- p + theme(legend.justification = c(1, 1), legend.position = c(1, 1))
  }

  print(p)

}


