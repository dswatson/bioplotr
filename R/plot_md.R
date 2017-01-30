#' Create MD plot of log2 fold changes against mean expression/methylation
#'
#' @param dat Data frame or matrix representing the results of a test for
#'   differential expression or methylation, such as the output of a call to
#'   \code{limma::topTable, edgeR::topTags}, or \code{DESeq2::results}.
#'   Alternatively, any object with columns for log fold changes and FDR.
#' @param fdr Threshold for declaring a probe differentially expressed or methylated.
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
#' This function displays the results of a differential expression or methylation
#' test as an MD plot, with colors distinguishing probes that meet the
#' user-defined FDR threshold. These figures can be used to evaluate the symmetry,
#' magnitude, and significance of fold changes for a given experiment.
#'
#' @examples
#' df <- data.frame(logFC   = c(rnorm(50, 0, 10), rnorm(4950)),
#'                  AvgExpr = rowMeans(matrix(rnorm(5000), nrow = 5000, ncol = 10)),
#'                  p.value = pnorm(-abs(logFC)),
#'                  FDR     = p.adjust(p.value, method = "fdr"))
#' plot_md(df)
#'
#' library(limma)
#' DE_genes <- cbind(matrix(rnorm(250, 5, 1), nrow = 50, ncol = 5),
#'                   matrix(rnorm(250), nrow = 50, ncol = 5))
#' mat <- rbind(DE_genes, matrix(rnorm(45500), nrow = 4550, ncol = 10))
#' treat <- gl(n = 2, k = 5, labels = c("A", "B"))
#' des <- model.matrix(~ treat)
#' fit <- eBayes(lmFit(mat, des))
#' top <- topTable(fit, number = Inf)
#' plot_md(top)
#'
#' @export
#' @importFrom purrr map_lgl
#' @import dplyr
#' @import ggplot2
#' @importFrom plotly ggplotly
#'

plot_md <- function(dat,
                    fdr = 0.05,
                 ptsize = 0.25,
                   main = NULL,
                 legend = 'outside',
                  hover = FALSE,
                 probes = NULL) {

  # Preliminaries
  dat <- as.data.frame(dat)
  lfc <- c('log2FoldChange', 'logFC')
  if (any(lfc %in% colnames(dat))) {
    j <- intersect(lfc, colnames(dat))
    colnames(dat)[colnames(dat) == j] <- 'logFC'
  } else {
    stop('dat must include a log fold change column. Recognized colnames for this ',
         'vector include "logFC" and "log2FoldChange". Make sure that dat includes ',
         'exactly one such colname.')
  }
  avg <- c('AvgMeth', 'AveExpr', 'logCPM', 'baseMean', 'AvgExpr')
  if (any(avg %in% colnames(dat))) {
    j <- intersect(avg, colnames(dat))
    colnames(dat)[colnames(dat) == j] <- 'AvgExpr'
  } else {
    stop('dat must include a column for average expression or methylation by ',
         'probe. Recognized colnames for this vector include "AvgExpr", "AvgMeth", ',
         '"AveExpr", "logCPM", and "baseMean". Make sure that dat includes exactly ',
         'one such colname.')
  }
  q <- c('adj.P.Val', 'FDR', 'padj', 'q.value')
  if (any(q %in% colnames(dat))) {
    j <- intersect(q, colnames(dat))
    colnames(dat)[colnames(dat) == j] <- 'q.value'
  } else {
    stop('dat must include a column for adjusted p-values. Recognized colnames ',
         'for this vector include "q.value", "adj.P.Val", "FDR", "padj", and "FDR". ',
         'Make sure that dat includes exactly one such colname.')
  }
  if (is.null(main)) {
    main <- 'MD Plot'
  }
  if (!legend %in% c('outside', 'bottomleft', 'bottomright', 'topleft', 'topright')) {
    stop('legend must be one of "outside", "bottomleft", "bottomright", ',
         '"topleft", or "topright".')
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
    mutate(is.DE = map_lgl(q.value, function(q) ifelse(q < fdr, TRUE, FALSE))) %>%
    select(Probe, AvgExpr, logFC, is.DE) %>%
    na.omit()

  # Build plot
  suppressWarnings(
    p <- ggplot(df, aes(AvgExpr, logFC, text = Probe)) +
      labs(title = main,
               x = expression(mu),
               y = expression('log'[2]*' Fold Change')) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
  )
  if (sum(df$is.DE == TRUE) == 0) {
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
  if (!hover) {
    print(p)
  } else {
    p <- ggplotly(p, tooltip = 'text', height = 600, width = 650)
    print(p)
  }

}


