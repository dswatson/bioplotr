#' Plot the mean-variance trend of an omic data matrix
#'
#' @param dat Omic data matrix with rows corresponding to probes and columns
#'   to samples. Presumed to be normalized prior to visualization.
#' @param trans Data transformation to be applied to the probewise standard
#'   deviations. Must be either \code{"log"} or \code{"sqrt"}. The former is
#'   recommended for microarrays or any other platform in which data are
#'   approximately log-normally distributed. The latter is recommended for
#'   sequencing and count data.
#' @param main Optional plot title.
#' @param hover Show probe name by hovering mouse over data point? If \code{TRUE},
#'   the plot is rendered in HTML and will either open in your browser's graphic
#'   display or appear in the RStudio viewer.
#' @param probes String specifying the name of the column in which to find the probe
#'   identifiers, assuming they aren't \code{rownames(dat)}. Only relevant if
#'   \code{hover = TRUE}.
#'
#' @details
#' This function plots each probe's mean value against its standard deviation,
#' following the appropriate data transformation. A lowess curve is additionally
#' fit to the data.
#'
#' @examples
#' mat <- matrix(rnorm(5000), nrow = 1000, ncol = 5)
#' plot_mean_var(mat, trans = "log")
#'
#' library(edgeR)
#' mat <- matrix(rnbinom(5000, mu = 4, size = 1), nrow = 1000, ncol = 5)
#' mat <- cpm(y, log = TRUE)
#' plot_mean_var(mat, trans = "sqrt")
#'
#' @export
#' @importFrom matrixStats rowSds
#' @import dplyr
#' @import ggplot2
#' @importFrom plotly ggplotly
#'

plot_mean_var <- function(dat,
                          trans,
                          ptsize = 0.25,
                          main   = NULL,
                          hover  = FALSE,
                          probes = NULL) {

  # Preliminaries
  if (!trans %in% c('log', 'sqrt')) {
    stop('trans must be specified as either "log" or "sqrt"')
  }
  if (trans == 'log') {
    vars <- log2(rowSds(dat))
    ylab <- expression('log'[2]*sigma)
    if (is.null(main)) {
      main <- 'Normalized Expression'
    }
  } else if (trans == 'sqrt') {
    vars <- sqrt(rowSds(dat))
    ylab <- expression(sqrt(sigma))
    if (is.null(main)) {
      main <- 'Normalized Counts'
    }
  }
  if (is.null(probes)) {
    if (is.null(rownames(dat))) {
      probes <- seq_along(nrow(dat))
    } else {
      probes <- rownames(dat)
    }
  } else {
    if (!probes %in% colnames(dat)) {
      stop(paste0('Column "', probes, '" not found'))
    } else {
      probes <- dat[, colnames(dat) == probes]
    }
  }

  # Tidy
  df <- data_frame(Probe = probes,
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
    theme(plot.title = element_text(hjust = 0.5))

  # Output
  if (hover == FALSE) {
    print(p)
  } else {
    p <- ggplotly(p, tooltip = 'text', height = 600, width = 600)
    print(p)
  }

}

# When shiny-ified: add slider for smoother span

