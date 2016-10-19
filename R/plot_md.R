#' Create MD plot of \eqn{log_2} fold changes against mean expression or methylation
#'
#' @param res Data frame or matrix representing the results of a test for
#'   differential expression or methylation, such as the output of a call to
#'   limma::topTable, edgeR::topTags, or DESeq2::results. Alternatively, any
#'   object with columns for log fold changes, and FDR.
#' @param dat A count
#' @param type String specifying data type. Must be one of either
#'   \code{"microarray"}, \code{"RNA-seq"}, or \code{"methylation"}.
#' @param fdr Threshold for declaring a probe differentially expressed or methylated.
#' @param ptsize Size of data points in the plot.
#' @param main Optional plot title.
#' @param legend Legend position. Must be one of \code{"outside",
#'   "bottomleft", "bottomright", "topleft",} or \code{"topright"}.
#'
#' @details
#' This function displays the results of a differential expression or methylation
#' test as an MD plot, with colors distinguishing probes that meet the
#' user-defined FDR threshold. These figures can be used to evaluate the symmetry,
#' magnitude, and significance of fold changes for a given experiment.
#'
#' @examples
#' library(dplyr)
#' df <- data_frame(logFC   = c(rnorm(50, 0, 10), rnorm(4950)),
#'                  AvgExpr = rowMeans(matrix(rnorm(5000), nrow = 5000, ncol = 10)),
#'                  p.value = pnorm(-abs(logFC)),
#'                  FDR     = p.adjust(p.value, method = "fdr"))
#' plot_md(df)
#'
#' library(limma)
#' DE_genes <- cbind(matrix(rnorm(250, 5, 1), nrow = 50, ncol = 5),
#'                   matrix(rnorm(250), nrow = 50, ncol = 5))
#' mat <- rbind(DE_genes, matrix(rnorm(45500), nrow = 4550, ncol = 10))
#' treat <- rep(c("A", "B"), each = 5)
#' des <- model.matrix(~ treat)
#' fit <- eBayes(lmFit(mat, des))
#' top <- topTable(fit, number = Inf)
#' plot_md(top)
#'
#' @export
#' @importFrom purrr map_lgl
#' @import dplyr
#' @import ggplot2
#'

plot_md <- function(dat,
                    fdr    = 0.05,
                    ptsize = 0.25,
                    main   = NULL,
                    legend = 'outside') {

  dat <- as_data_frame(na.omit(dat))
  if ('adj.P.Val' %in% colnames(dat)) {
    dat <- dat %>% rename(q.value = adj.P.Val)
  }
  if ('FDR' %in% colnames(dat)) {
    dat <- dat %>% rename(q.value = FDR)
  }
  if ('padj' %in% colnames(dat)) {
    dat <- dat %>% rename(q.value = padj)
  }
  if (!'adj.P.Val' %in% colnames(dat) &
      !'FDR' %in% colnames(dat) &
      !'padj' %in% colnames(dat) &
      !'q.value' %in% colnames(dat)) {
    stop('dat must include a column for adjusted p-values. Recognized colnames
  for this vector include "adj.P.Val", "padj", "FDR", and "q.value". Make sure
  that dat includes exactly one such colname.')
  }
  if ('AvgMeth' %in% colnames(dat)) {
    dat <- dat %>% rename(AvgExpr = AvgMeth)
  }
  if ('AveExpr' %in% colnames(dat)) {
    dat <- dat %>% rename(AvgExpr = AveExpr)
  }
  if ('logCPM' %in% colnames(dat)) {
    dat <- dat %>% rename(AvgExpr = logCPM)
  }
  if ('baseMean' %in% colnames(dat)) {
    dat <- dat %>% mutate(AvgExpr = log2(baseMean))
  }
  if (!'AveExpr' %in% colnames(dat) &
      !'logCPM' %in% colnames(dat) &
      !'baseMean' %in% colnames(dat) &
      !'AvgExpr' %in% colnames(dat)) {
    stop('dat must include a column for average expression or methylation by
  probe. Recognized colnames for this vector include "AvgExpr", "AvgMeth",
  "AveExpr", "logCPM", and "baseMean". Make sure that dat includes exactly one
  such colname.')
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
    main <- 'MD Plot'
  }

  test <- function(q) ifelse(q < fdr, TRUE, FALSE)
  df <- as_data_frame(dat) %>%
    mutate(is.DE = map_lgl(q.value, test)) %>%
    select(AvgExpr, logFC, is.DE)

  if (sum(df$is.DE == TRUE) == 0) {
    warning('dat returned no differentially expressed/methylated probes at your
  selected fdr threshold. To color points by differential expression/methylation,
  consider raising your fdr cutoff.')
    p <- ggplot(df, aes(AvgExpr, logFC)) +
      geom_point(size = ptsize, alpha = 0.25)
  } else {
    p <- ggplot(df, aes(AvgExpr, logFC, color = is.DE)) +
      geom_point(size = ptsize, alpha = 0.25) +
      scale_colour_manual(name   = expression(italic(q)*'-value'),
                          labels = c(paste('\u2265', fdr), paste('<', fdr)),
                          values = c('black', 'red')) +
      guides(col = guide_legend(reverse = TRUE))
  }

  p <- p + labs(title = main,
                x = expression('Mean log'[2]*' Normalised Expression'),
                y = expression('log'[2]*' Fold Change')) +
    theme_bw()

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


