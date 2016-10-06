#' Create density plots by sample
#'
#' @param dat Omic data matrix with rows corresponding to probes and columns 
#'   to samples.
#' @param group Optional factor or character vector of length equal to sample size. 
#'   Levels are used to color density curves.
#' @param type Optional string specifying omic data type. Currently supports 
#'   \code{"microarray", "RNA-seq",} or \code{"methylation"}.
#' @param xlab Label for x-axis. If left \code{NULL}, this defaults to 
#'   \eqn{log_2} expression for \code{type = "microarray"}, \eqn{log_2} CPM for
#'   \code{type = "RNA-seq"}, and Beta for \code{type = "methylation"}. At least 
#'   one of \code{type} or \code{xlab} must be specified.
#' @param main Optional plot title. 
#' @param legend Legend position. Must be one of \code{"outside", 
#'   "bottomleft", "bottomright", "topleft",} or \code{"topright"}.
#'
#' @details 
#' This function displays each sample's distribution of expression or methylation 
#' values in a single plot. It is especially helpful when contrasting pre- and 
#' post-normalization matrices. It may additionally be used to inspect for batch effects 
#' or other associations with phenotypic categories by using the \code{group} argument.
#' 
#' @examples
#' mat <- matrix(rnorm(5000), nrow = 1000, ncol = 5)
#' plot_density(mat, type = "microarray")
#' 
#' library(DESeq2)
#' mat <- cbind(matrix(rnbinom(5000, mu = 4, size = 1), nrow = 1000, ncol = 5),
#'              matrix(rnbinom(5000, mu = 4, size = 10), nrow = 1000, ncol = 5))
#' mat <- rlog(mat)
#' batch <- rep(c("A", "B"), each = 5)
#' plot_density(mat, group = batch, type = "RNA-seq", 
#'              main = "rlog Transformed Counts")
#'
#' @export
#' @importFrom tidyr gather 
#' @import dplyr
#' @import ggplot2
#'

plot_density <- function(dat, 
                         group  = NULL, 
                         type   = NULL,
                         xlab   = NULL,
                         main   = NULL, 
                         legend = 'outside') {

  densities <- gather(as_tibble(dat), Sample, Expression) %>%
    group_by(Sample) %>%
    do(Expression = density(.$Expression)$x,
       Density    = density(.$Expression)$y) 
  densities <- densities[match(colnames(as_tibble(dat)), densities$Sample), ]

  df <- data_frame(Sample     = rep(densities$Sample, each = 512),
                   Expression = unlist(densities$Expression),
                   Density    = unlist(densities$Density))

  if (!is.null(group)) {
    df <- mutate(df, Group = rep(group, each = 512))
  } 
  if (is.null(main) & is.null(group)) {
    main <- 'Density By Sample'
  }  
  if (is.null(main) & !is.null(group)) {
    main <- 'Density By Group'
  }
  if (is.null(type) & is.null(xlab)) {
    stop('Either data type or xlab must be provided')
  }
  if (type == 'microarray') {
    xlab <- expression('log'[2]*' Expression')
  } 
  if (type == 'RNA-seq') {
    xlab <- expression('log'[2]*' Counts Per Million')
  }
  if (type == 'methylation') {
    xlab <- 'Beta'
  }

  p <- ggplot(df, aes(Expression, Density, group = Sample)) + 
    labs(title = main, x = xlab) + 
    theme_bw() 

  if (!is.null(group)) {
    p <- p + geom_path(aes(colour = Group))
  } else {
    p <- p + geom_path(aes(colour = Sample))
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


