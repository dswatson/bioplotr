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
#' @param hover Show probe name by hovering mouse over data point? If \code{TRUE},
#'   the plot is rendered in HTML and will either open in your browser's graphic
#'   display or appear in the RStudio viewer. The plot can also be embedded in an
#'   HTML doc using Rmarkdown so long as \code{knitr = TRUE} and the code chunk
#'   option \code{plotly} is also set to \code{TRUE}.
#' @param knitr Set this to \code{TRUE} if you want to embed a plotly object (viz.,
#'   the \code{plot_md} output when \code{hover = TRUE}) in an HTML doc. Make
#'   sure to set \code{plotly = TRUE} in the corresponding code chunk options.
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
                    legend = 'outside',
                    hover  = FALSE,
                    knitr  = FALSE) {

  # Preliminaries
  dat <- as_data_frame(dat)
  q <- c('adj.P.Val', 'FDR', 'padj')
  for (i in q) {
    if (i %in% colnames(dat)) {
      colnames(dat)[colnames(dat) == i] <- 'q.value'
    }
  }
  if (all(!q %in% colnames(dat))) {
    stop('dat must include a column for adjusted p-values. Recognized colnames
         for this vector include "q.value", "adj.P.Val", "FDR", "padj", and "FDR".
         Make sure that dat includes exactly one such colname.')
  }
  avg <- c('AvgMeth', 'AveExpr', 'logCPM', 'baseMean')
  for (i in avg) {
    colnames(dat)[colnames(dat) == i] <- 'AvgExpr'
  }
  if (all(!avg %in% colnames(dat))) {
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
  # Add a prelim for GeneSymbol
  if (is.null(main)) {
    main <- 'MD Plot'
  }

  # Tidy
  test <- function(q) ifelse(q < fdr, TRUE, FALSE)
  df <- as_data_frame(dat) %>%
    mutate(is.DE = map_lgl(q.value, test)) %>%
    select(GeneSymbol, AvgExpr, logFC, is.DE) %>%
    na.omit()

  # Build plot
  if (sum(df$is.DE == TRUE) == 0) {
    warning('dat returned no differentially expressed/methylated probes at your
  selected fdr threshold. To color points by differential expression/methylation,
  consider raising your fdr cutoff.')
    p <- ggplot(df, aes(AvgExpr, logFC))
  } else {
    p <- ggplot(df, aes(AvgExpr, logFC, color = is.DE))
      scale_colour_manual(name   = expression(italic(q)*'-value'),
                          labels = c(paste('\u2265', fdr), paste('<', fdr)),
                          values = c('black', 'red')) +
      guides(col = guide_legend(reverse = TRUE))
  }
  p <- p + suppressWarnings(geom_point(aes(text = paste('Gene:', GeneSymbol)),
                                       size = ptsize, alpha = 0.25)) +
    labs(title = main,
         x = expression('Mean Expression'),
         y = expression('log'[2]*' Fold Change')) +
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


