#' Volcano Plot
#'
#' This function plots the log2 fold changes against the -log10 \emph{p}-values
#' for a given test of differential expression/methylation.
#'
#' @param dat Data frame or matrix representing the results of a test for
#'   differential expression or methylation, such as the output of a call to
#'   \code{limma::topTable, edgeR::topTags}, or \code{DESeq2::results}.
#'   Alternatively, any object with columns for \emph{p}-values, log fold changes,
#'   and FDR.
#' @param fdr Threshold for declaring a probe differentially expressed/methylated.
#' @param ptsize Size of data points in the plot.
#' @param main Optional plot title.
#' @param legend Legend position. Must be one of \code{"outside",
#'   "bottomleft", "bottomright", "topleft",} or \code{"topright"}.
#' @param hover Show probe name by hovering mouse over data point? If \code{TRUE},
#'   the plot is rendered in HTML and will either open in your browser's graphic
#'   display or appear in the RStudio viewer.
#' @param probes String specifying the name of the column in which to find the probe
#'   identifiers. Only relevant if \code{hover = TRUE}.
#'
#' @details
#' Volcano plots visualize the relationship between each probe's log2 fold change and
#' -log10 \emph{p}-value for a given test of differential expression/methylation.
#' Probes are colored to distinguish between those that do and do not meet a
#' user-defined FDR threshold. These figures help to evaluate the symmetry,
#' magnitude, and significance of effects for a given experiment.
#'
#' @examples
#' library(limma)
#' DE_genes <- cbind(matrix(rnorm(50 * 5, mean = 5), nrow = 50, ncol = 5),
#'                   matrix(rnorm(50 * 5), nrow = 50, ncol = 5))
#' eset <- rbind(DE_genes, matrix(rnorm(4950 * 10), nrow = 4950, ncol = 10))
#' treat <- rep(c("A", "B"), each = 5)
#' des <- model.matrix(~ treat)
#' fit <- eBayes(lmFit(eset, des))
#' top <- topTable(fit, number = Inf)
#' plot_volcano(top)
#'
#' @export
#' @importFrom purrr map_lgl
#' @import dplyr
#' @import ggplot2
#' @importFrom plotly ggplotly
#'

plot_volcano <- function(dat,
                         fdr = 0.05,
                      ptsize = 0.25,
                        main = NULL,
                      legend = 'outside',
                       hover = FALSE,
                      probes = NULL) {

  # Preliminaries
  dat <- as.data.frame(dat)
  lfc <- c('log2FoldChange', 'logFC')
  if (sum(lfc %in% colnames(dat)) == 1) {
    colnames(dat)[colnames(dat) %in% lfc] <- 'logFC'
  } else {
    stop('dat must include a log fold change column. Recognized colnames for this ',
         'vector include "logFC" and "log2FoldChange". Make sure that dat includes ',
         'exactly one such colname.')
  }
  p <- c('P.Value', 'PValue', 'pvalue', 'p.value')
  if (sum(p %in% colnames(dat)) == 1) {
    colnames(dat)[colnames(dat) %in% p] <- 'p.value'
  } else {
    stop('dat must include a p-value column. Recognized colnames for this vector ',
         'include "p.value", "P.Value", "PValue", and "pvalue". Make sure that dat',
         'includes exactly one such colname.')
  }
  q <- c('adj.P.Val', 'FDR', 'padj', 'q.value')
  if (sum(q %in% colnames(dat)) == 1) {
    colnames(dat)[colnames(dat) %in% q] <- 'q.value'
  } else {
    stop('dat must include a column for adjusted p-values. Recognized colnames ',
         'for this vector include "q.value", "adj.P.Val", "FDR", "padj", and "FDR". ',
         'Make sure that dat includes exactly one such colname.')
  }
  if (is.null(main)) {
    main <- 'Volcano Plot'
  }
  if (!legend %in% c('outside', 'bottomleft', 'bottomright', 'topleft', 'topright')) {
    stop('legend must be one of "outside", "bottomleft", "bottomright", ',
         '"topleft", or "topright".')
  }
  if (is.null(probes)) {
    if (identical(rownames(dat), as.character(seq_len(nrow(dat)))) ||
        is.null(rownames(dat))) {
      dat <- dat %>% mutate(Probe = row_number())
    } else {
      dat <- dat %>% mutate(Probe = rownames(dat))
    }
  } else {
    if (!probes %in% colnames(dat)) {
      stop(paste0('Column "', probes, '" not found.'))
    } else {
      colnames(dat)[colnames(dat) == probes] <- 'Probe'
    }
  }

  # Tidy data
  df <- dat %>%
    mutate(logP = -log10(p.value),
          is.DE = map_lgl(q.value, function(q) ifelse(q < fdr, TRUE, FALSE))) %>%
    select(Probe, logFC, logP, is.DE) %>%
    na.omit()

  # Build plot
  suppressWarnings(
    p <- ggplot(df, aes(logFC, logP, text = Probe)) +
      labs(title = main,
               x = expression(log[2]*' Fold Change'),
               y = expression(~-log[10](italic(p)))) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
  )
  if (!all(df$is.DE)) {          # Color pts by differential expression?
    warning('No probe meets your fdr threshold. To color data points by differential ',
            'expression/methylation, consider raising your fdr cutoff.')
    p <- p + geom_point(size = ptsize, alpha = 0.25)
  } else {
    p <- p + geom_point(aes(color = is.DE), size = ptsize, alpha = 0.25) +
      scale_colour_manual(name = expression(italic(q)*'-value'),
                        labels = c(paste('\u2265', fdr), paste('<', fdr)),
                        values = c('black', 'red')) +
      guides(col = guide_legend(reverse = TRUE))
  }
  if (legend == 'bottomleft') {  # Locate legend
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
  if (!hover) {
    print(p)
  } else {
    p <- ggplotly(p, tooltip = 'text', height = 600, width = 650)
    print(p)
  }

}


