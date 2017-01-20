#' Create QQ plot of expected vs. observed -log10 \emph{p}-values
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
#' @param probes String specifying the name of the column in which to find the probe
#'   identifiers, assuming they aren't \code{rownames(dat)}. Only relevant if
#'   \code{hover = TRUE}.
#'
#' @details
#' This function displays a scatterplot of the expected vs. observed quantiles
#' of a \emph{p}-value distribution following -log10 transform. If the black points
#' deviate too sharply from the red line, especially at low expected values of
#' -log10(\emph{p}), then it suggests a violation of the statistical assumptions
#' upon which the test was based.
#'
#' @examples
#' df <- data.frame(p.value = runif(10000))
#' plot_qq(df)
#'
#' library(limma)
#' DE_genes <- cbind(matrix(rnorm(250, 5, 1), nrow = 50, ncol = 5),
#'                   matrix(rnorm(250), nrow = 50, ncol = 5))
#' mat <- rbind(DE_genes, matrix(rnorm(45500), nrow = 4550, ncol = 10))
#' treat <- gl(n = 2, k = 5, labels = c("A", "B"))
#' des <- model.matrix(~ treat)
#' fit <- eBayes(lmFit(mat, des))
#' top <- topTable(fit, number = Inf)
#' plot_qq(top)
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
  dat <- as.data.frame(dat)
  p <- c('P.Value', 'PValue', 'pvalue', 'p.value')
  for (i in p) {
    if (i %in% colnames(dat)) {
      colnames(dat)[colnames(dat) == i] <- 'p.value'
    }
  }
  if (all(!p %in% colnames(dat))) {
    stop('dat must include a p-value column. Recognized colnames for this ',
         'vector include "p.value", "P.Value", "PValue", and "pvalue". Make sure ',
         'that dat includes exactly one such colname')
  }
  if (is.null(main)) {
    main <- 'Q-Q Plot'
  }
  if (!legend %in% c('outside', 'bottomleft', 'bottomright', 'topleft', 'topright')) {
    stop('legend must be one of "outside", "bottomleft", "bottomright", ',
         '"topleft", or "topright"')
  }
  if (is.null(probes)) {
    if (is.null(rownames(dat))) {
      dat %>% mutate(Probe = row_number())
    } else {
      dat %>% mutate(Probe = rownames(dat))
    }
  } else {
    if (!probes %in% colnames(dat)) {
      stop(paste0('Column "', probes, '" not found'))
    } else {
      colnames(dat)[colnames(dat) == probes] <- 'Probe'
    }
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
         x = expression('Expected'~-log[10](italic(p))),
         y = expression('Observed'~-log[10](italic(p)))) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

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
    p <- ggplotly(p, tooltip = 'text', height = 600, width = 600)
    print(p)
  }

}


