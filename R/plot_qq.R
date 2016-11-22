#' Create QQ plot of expected vs. observed \emph{p}-values
#'
#' @param dat Data frame or matrix representing the results of a test for
#'   differential expression or methylation, such as the output of a call to
#'   \code{limma::topTable, edgeR::topTags}, or \code{DESeq2::results}.
#'   Alternatively, any object with a column for \emph{p}-values.
#' @param ptsize Size of data points in the plot.
#' @param main Optional plot title.
#' @param legend Legend position. Must be one of \code{"outside",
#'   "bottomleft", "bottomright", "topleft",} or \code{"topright"}.
#' @param hover Show probe name by hovering mouse over data point? If \code{TRUE},
#'   the plot is rendered in HTML and will either open in your browser's graphic
#'   display or appear in the RStudio viewer.
#' @param probes String specifying the colname in which to find the probe names.
#'   Only relevant if \code{hover = TRUE}.
#'
#' @details
#' This function displays a scatterplot of the expected vs. observed quantiles
#' of a \emph{p}-value distribution following a -log10 transform. If the black points
#' deviate too sharply from the red line, especially at low expected values of
#' -log10(\emph{p}), then it suggests that one or several statistical assumptions
#' upon which the test was based have probably been violated.
#'
#' @examples
#' df <- data.frame(p.value = runif(10000))
#' plot_qq(df)
#'
#' @export
#' @import dplyr
#' @import ggplot2
#' @importFrom plotly ggplotly
#'

plot_qq <- function(dat,
                    ptsize = 0.25,
                    main   = NULL,
                    legend = 'outside',
                    hover  = FALSE,
                    probes = NULL) {

  # Preliminaries
  dat <- as_data_frame(dat)
  if (is.null(probes)) {
    dat <- dat %>% mutate(Probe = row_number())
  } else {
    if (!probes %in% colnames(dat)) {
      stop('dat must include a colname that matches the string specified in
  the probes argument.')
    } else {
      colnames(dat)[colnames(dat) == probes] <- 'Probe'
    }
  }
  p <- c('P.Value', 'PValue', 'pvalue', 'p.value')
  for (i in p) {
    if (i %in% colnames(dat)) {
      colnames(dat)[colnames(dat) == i] <- 'p.value'
    }
  }
  if (all(!p %in% colnames(dat))) {
    stop('dat must include a p-value column. Recognized colnames for this
  vector include "p.value", "P.Value", "PValue", and "pvalue". Make sure
  that dat includes exactly one such colname.')
  }

  # Tidy
  df <- dat %>%
    mutate(Observed = -log10(sort(p.value, decreasing = FALSE)),
           Expected = -log10(ppoints(length(p.value)))) %>%
    select(Probe, Observed, Expected)

  # Basic plot
  p <- ggplot(df, aes(Expected, Observed)) +
    suppressWarnings(geom_point(aes(text = Probe), size = ptsize)) +
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



