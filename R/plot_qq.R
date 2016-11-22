#' Create QQ plot of expected vs. observed *p*-values
#'
#' @param dat Data frame or matrix representing the results of a test for
#'   differential expression or methylation, such as the output of a call to
#'   limma::topTable, edgeR::topTags, or DESeq2::results. Alternatively, any
#'   object with columns for probes and *p*-values.
#' @param main Optional plot title.
#' @param legend Legend position. Must be one of \code{"outside",
#'   "bottomleft", "bottomright", "topleft",} or \code{"topright"}.
#' @param hover Show probe name by hovering mouse over data point? If \code{TRUE},
#'   the plot is rendered in HTML and will either open in your browser's graphic
#'   display or appear in the RStudio viewer. The plot can also be embedded in an
#'   HTML doc using Rmarkdown so long as \code{knitr = TRUE}.
#' @param knitr Set this to \code{TRUE} if you want to embed a plotly object (viz.,
#'   the \code{plot_md} output when \code{hover = TRUE}) in an HTML doc.
#'
#' @details
#'
#'
#' @examples
#'
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'

plot_qq <- function(dat,
                    main   = NULL,
                    legend = 'outside',
                    hover  = FALSE,
                    knitr  = FALSE) {

  # Preliminaries


  # Tidy
  df <- data_frame(Gene = dat$GeneSymbol,
                   Observed = -log10(sort(dat$p.value, decreasing = FALSE)),
                   Expected = -log10(ppoints(length(dat$p.value))))

  # Basic plot
  p <- ggplot(df, aes(Expected, Observed)) +
    suppressWarnings(geom_point(aes(text = paste('Gene:', Gene)), size = 0.25)) +
    geom_abline(intercept = 0, slope = 1, color = 'red') +
    labs(title = main,
         x = expression('Expected -log'[10](italic(p))),
         y = expression('Observed -log'[10](italic(p)))) +
    theme_bw() +
    theme(plot.title = element_text(hjust = .5))

  # Legend location
  if (legend == 'bottomleft') {
    p <- p + theme(legend.justification = c(0.01, 0.01),
                   legend.position = c(0.01, 0.01))
  } else if (legend == 'bottomright') {
    p <- p + theme(legend.justification = c(0.99, 0.01),
                   legend.position = c(0.99, 0.01))
  } else if (legend == 'topleft') {
    p <- p + theme(legend.justification = c(0.01, 0.99),
                   legend.position = c(0.01, 0.99))
  } else if (legend == 'topright') {
    p <- p + theme(legend.justification = c(0.99, 0.99),
                   legend.position = c(0.99, 0.99))
  }

  # Output
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

# Add: norm, chisq, pois, nbinom, beta



