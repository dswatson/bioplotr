#' Plot the mean-variance trend of an omic data matrix
#'
#' @param dat Omic data matrix with rows corresponding to probes and columns
#'   to samples.
#' @param type String specifying data type. Must be one of either
#'   \code{"microarray"} or \code{"RNA-seq"}.
#' @param main Optional plot title.
#' @param hover Show probe name by hovering mouse over data point? If \code{TRUE},
#'   the plot is rendered in HTML and will either open in your browser's graphic
#'   display or appear in the RStudio viewer.
#' @param probes String specifying the name of the column in which to find the probe
#'   identifiers. Only relevant if \code{hover = TRUE}.
#'
#' @details
#' This function plots each probe's mean value against either the logarithm of
#' its standard deviation (if \code{type = "microarray"}) or its square root
#' (if \code{type = "RNA-seq"}). A lowess curve is additionally fit to the data.
#'
#' @examples
#' mat <- matrix(rnorm(5000), nrow = 1000, ncol = 5)
#' plot_mean_var(mat, type = "microarray")
#'
#' library(edgeR)
#' mat <- matrix(rnbinom(5000, mu = 4, size = 1), nrow = 1000, ncol = 5)
#' mat <- cpm(y, log = TRUE)
#' plot_mean_var(mat, type = "RNA-seq")
#'
#' @export
#' @importFrom matrixStats rowSds
#' @import dplyr
#' @import ggplot2
#' @importFrom plotly ggplotly
#'

plot_mean_var <- function(dat,
                          type,
                          ptsize = 0.25,
                          main   = NULL,
                          hover  = FALSE,
                          probes = NULL) {

  # Preliminaries
  if (!type %in% c('microarray', 'RNA-seq')) {
    stop('type must be specified as either "microarray" or "RNA-seq".')
  }
  if (type == 'microarray') {
    vars <- log2(rowSds(dat))
    ylab <- expression('log'[2]*sigma)
    if (is.null(main)) {
      main <- expression('log'[2]*' Expression')
    }
  } else  if (type == 'RNA-seq') {
    vars <- sqrt(rowSds(dat))
    ylab <- expression(sqrt(sigma))
    if (is.null(main)) {
      main <- expression('log'[2]*' Counts Per Million') #MREH
    }
  }

  # Tidy
  df <- data_frame(Probe = rownames(dat),
                   Mean  = rowMeans(dat),
                   Var   = vars)
  lo <- lowess(x = df$Mean, y = df$Var, f = 0.5)
  df <- df %>% mutate(lo.x = lo[['x']],
                      lo.y = lo[['y']])

  # Build plot
  p <- ggplot(df) +
    suppressWarnings(geom_point(aes(Mean, Var, text = Probe),
                                size = ptsize, alpha = 0.25)) +
    geom_smooth(aes(lo.x, lo.y), size = 0.5) +
    labs(title = main, x = expression(mu), y = ylab) +
    theme_bw() +
    theme(plot.title = element_text(hjust = .5))

  # Output
  if (hover == FALSE) {
    print(p)
  } else {
    p <- ggplotly(p, tooltip = 'text', height = 600, width = 600)
    print(p)
  }

}



