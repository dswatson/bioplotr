#' Create volcano plot of \eqn{log_2} fold changes against \eqn{-log_10 p}-values
#'
#' @param dat Data frame or matrix representing the results of a test for
#'   differential expression or methylation, such as the output of a call to
#'   limma::topTable, edgeR::topTags, or DESeq2::results. Alternatively, any object
#'   with columns for p-values, log fold changes, and FDR.
#' @param fdr Threshold for declaring a probe differentially expressed/methylated.
#' @param ptsize Size of data points in the plot.
#' @param main Optional plot title.
#' @param legend Legend position. Must be one of \code{"outside",
#'   "bottomleft", "bottomright", "topleft",} or \code{"topright"}.
#'
#' @details
#' This function displays the results of a differential expression or methylation
#' test as a volcano plot, with colors distinguishing probes that meet the
#' user-defined FDR threshold. These figures can be used to evaluate the symmetry,
#' magnitude, and significance of fold changes for a given experiment.
#'
#' @examples
#' library(dplyr)
#' df <- data_frame(logFC   = c(rnorm(50, 0, 10), rnorm(4950)),
#'                  p.value = pnorm(-abs(logFC)),
#'                  FDR     = p.adjust(p.value, method = 'fdr'))
#' plot_volcano(df)
#'
#' library(limma)
#' DE_genes <- cbind(matrix(rnorm(250, 5, 1), nrow = 50, ncol = 5),
#'                   matrix(rnorm(250), nrow = 50, ncol = 5))
#' mat <- rbind(DE_genes, matrix(rnorm(45500), nrow = 4550, ncol = 10))
#' treat <- c(rep("A", 5), rep("B", 5))
#' des <- model.matrix(~ treat)
#' fit <- eBayes(lmFit(mat, des))
#' top <- topTable(fit, number = Inf)
#' plot_volcano(top)
#'
#' @export
#' @importFrom purrr map_lgl
#' @import dplyr
#' @import ggplot2
#'

plot_volcano <- function(dat,
                         fdr    = 0.05,
                         ptsize = 0.25,
                         main   = NULL,
                         legend = 'outside') {

  dat <- as_data_frame(dat)

  for (p in c('P.Value', 'PValue', 'pvalue')) {
    if (p %in% colnames(dat)) {
      colnames(dat)[colnames(dat) == p] <- 'p.value'
    }
  }
  if (!'P.Value' %in% colnames(dat) &
      !'pvalue'  %in% colnames(dat) &
      !'p.value' %in% colnames(dat)) {
    stop('dat must include a p-value column. Recognized colnames for this
  vector include "p.value", "P.Value", "PValue", and "pvalue". Make sure
  that dat includes exactly one such colname.')
  }

  for (q in c('adj.P.Val', 'FDR', 'padj')) {
    if (q %in% colnames(dat)) {
      colnames(dat)[colnames(dat) == q] <- 'q.value'
    }
  }
  if (!'adj.P.Val' %in% colnames(dat) &
      !'FDR' %in% colnames(dat) &
      !'padj' %in% colnames(dat) &
      !'q.value' %in% colnames(dat)) {
    stop('dat must include a column for adjusted p-values. Recognized colnames
  for this vector include "q.value", "adj.P.Val", "FDR", "padj", and "FDR".
  Make sure that dat includes exactly one such colname.')
  }

  if ('log2FoldChange' %in% colnames(dat)) {
    dat <- dat %>% rename(logFC = log2FoldChange)
  }
  if (!'log2FoldChange' %in% colnames(dat) &
      !'logFC' %in% colnames(dat)) {
    stop('dat must include a log fold change column. Recognized colnames for
  this vector include "logFC" and "log2FoldChange". Make sure that dat includes
  exactly one such colname.')
  }

  if (is.null(main)) {
    main <- 'Volcano Plot'
  }

  test <- function(q) ifelse(q < fdr, TRUE, FALSE)
  df <- dat %>%
    mutate(is.DE = map_lgl(q.value, test),
           logP  = -log10(p.value)) %>%
    select(logFC, logP, is.DE) %>%
    na.omit()

  if (sum(df$is.DE == TRUE) == 0) {
    warning('dat returned no differentially expressed/methylated probes at your
  selected fdr threshold. To color points by differential expression/methylation,
  consider raising your fdr cutoff.')
    p <- ggplot(df, aes(logFC, logP)) +
      geom_point(size = ptsize, alpha = 0.25)
  } else {
    p <- ggplot(df, aes(logFC, logP, color = is.DE)) +
      geom_point(size = ptsize, alpha = 0.25) +
      scale_colour_manual(name   = expression(italic(q)*'-value'),
                          labels = c(paste('\u2265', fdr), paste('<', fdr)),
                          values = c('black', 'red')) +
      guides(col = guide_legend(reverse = TRUE))
  }

  p <- p + labs(title = main,
    x = expression('log'[2]*' Fold Change'),
    y = expression('-log'[10]*' '*italic(p))) +
    theme_bw() +
    theme(plot.title = element_text(hjust = .5))

  if (legend == 'bottomleft') {
    p <- p + theme(legend.justification = c(0.01, 0.01),
                   legend.position = c(0.01, 0.01))
  } else if (legend == 'bottomright') {
    p <- p + theme(legend.justification = c(1, 0.01),
                   legend.position = c(0.99, 0.01))
  } else if (legend == 'topleft') {
    p <- p + theme(legend.justification = c(0.01, 0.99),
                   legend.position = c(0.01, 0.99))
  } else if (legend == 'topright') {
    p <- p + theme(legend.justification = c(0.99, 0.99),
                   legend.position = c(0.99, 0.99))
  }

  print(p)

}


