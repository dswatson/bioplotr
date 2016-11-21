#' Plot the mean-variance trend of a gene expression matrix
#'
#' @param dat Matrix of log expression intensities or raw counts.
#' @param type String specifying data type. Must be one of either
#'   \code{"microarray"} or \code{"RNA-seq"}.
#' @param main Optional plot title.
#' @param hover Show probe name by hovering mouse over data point? If \code{TRUE},
#'   the plot is rendered in HTML and will either open in your browser's graphic
#'   display or appear in the RStudio viewer. The plot can also be embedded in an
#'   HTML doc using Rmarkdown so long as \code{knitr = TRUE} and the code chunk
#'   option \code{plotly} is also set to \code{TRUE}.
#' @param knitr Set this to \code{TRUE} if you want to embed a plotly object (viz.,
#'   the \code{plot_mean_var} output when \code{hover = TRUE}) in an HTML doc. Make
#'   sure to set \code{plotly = TRUE} in the corresponding code chunk options.
#'
#' @details
#' This function plots each gene's mean expression against either the
#' logarithm of its standard deviation (if \code{type = "microarray"}) or
#' its square root (if \code{type = "RNA-seq"}). A lowess curve is additionally
#' fit to the data.
#'
#' @examples
#' mat <- matrix(rnorm(5000), nrow = 1000, ncol = 5)
#' plot_mean_var(mat, type = "microarray")
#'
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
                          knitr  = FALSE) {

  if (is.null(type) |
      !type %in% c('microarray', 'RNA-seq')) {
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
      main <- expression('log'[2]*' Counts Per Million')
    }
  }

  df <- data_frame(Gene = rownames(dat),
                   Mean = rowMeans(dat),
                   Var  = vars)
  lo <- lowess(x = df$Mean, y = df$Var, f = 0.5)
  df <- df %>% mutate(lo.x = lo[['x']],
                      lo.y = lo[['y']])

  p <- ggplot(df) +
    suppressWarnings(geom_point(aes(Mean, Var, text = paste('Gene:', Gene)),
                                size = ptsize, alpha = 0.25)) +
    geom_smooth(aes(lo.x, lo.y), size = 0.5) +
    labs(title = main, x = expression(mu), y = ylab) +
    theme_bw() +
    theme(plot.title = element_text(hjust = .5))

  if (hover == FALSE) {
    print(p)
  } else {
    if (knitr == FALSE) {
      p <- ggplotly(p, tooltip = 'text', width = 600, height = 500)
      print(p)
    } else {
      p <- ggplotly(p, tooltip = 'text', width = 600, height = 500,
                    session = 'knitr')
      print(p)
    }
  }

}



