#' Create taco plot of log2 fold changes against -log10 \emph{p}-values against mean expression/methylation
#'
#' @param dat Data frame or matrix representing the results of a test for
#'   differential expression or methylation, such as the output of a call to
#'   \code{limma::topTable, edgeR::topTags}, or \code{DESeq2::results}.
#'   Alternatively, any object with columns for \emph{p}-values, log fold changes,
#'   and FDR.
#' @param fdr Threshold for declaring a probe differentially expressed/methylated.
#' @param main Optional plot title.
#' @param legend Legend position. Must be one of \code{"outside",
#'   "bottomleft", "bottomright", "topleft",} or \code{"topright"}.
#' @param probes String specifying the name of the column in which to find the probe
#'   identifiers.
#'
#' @details
#' This function displays the results of a differential expression or methylation
#' test as a taco plot, which is essentially a cross between a volcano plot and an MD
#' plot, with colors distinguishing probes that meet the user-defined FDR threshold.
#' These figures can be used to evaluate the symmetry, magnitude, and significance of
#' fold changes for a given experiment.
#'
#' @examples
#' df <- data.frame(logFC   = c(rnorm(50, 0, 10), rnorm(4950)),
#'                  AvgExpr = rowMeans(matrix(rnorm(5000), nrow = 5000, ncol = 10)),
#'                  p.value = pnorm(-abs(logFC)),
#'                  FDR     = p.adjust(p.value, method = 'fdr'))
#' plot_taco(df)
#'
#' library(limma)
#' DE_genes <- cbind(matrix(rnorm(250, 5, 1), nrow = 50, ncol = 5),
#'                   matrix(rnorm(250), nrow = 50, ncol = 5))
#' mat <- rbind(DE_genes, matrix(rnorm(45500), nrow = 4550, ncol = 10))
#' treat <- gl(n = 2, k = 5, labels = c("A", "B"))
#' des <- model.matrix(~ treat)
#' fit <- eBayes(lmFit(mat, des))
#' top <- topTable(fit, number = Inf)
#' plot_taco(top)
#'
#' @export
#' @importFrom purrr map_chr
#' @import dplyr
#' @importFrom plotly plot_ly
#'

plot_taco <- function(dat,
                      fdr    = 0.05,
                      main   = NULL,
                      legend = 'outside',
                      probes = NULL) {

  # Preliminaries
  if (is.null(probes)) {
    if (is.null(rownames(dat))) {
      probes <- 1:nrow(dat)
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
  dat <- as_data_frame(dat)
  p <- c('P.Value', 'PValue', 'pvalue', 'p.value')
  for (i in p) {
    if (i %in% colnames(dat)) {
      colnames(dat)[colnames(dat) == i] <- 'p.value'
    }
  }
  if (all(!p %in% colnames(dat))) {
    stop('dat must include a p-value column. Recognized colnames for this vector ',
         'include "p.value", "P.Value", "PValue", and "pvalue". Make sure that dat',
         'includes exactly one such colname')
  }
  q <- c('adj.P.Val', 'FDR', 'padj', 'q.value')
  for (i in q) {
    if (i %in% colnames(dat)) {
      colnames(dat)[colnames(dat) == i] <- 'q.value'
    }
  }
  if (all(!q %in% colnames(dat))) {
    stop('dat must include a column for adjusted p-values. Recognized colnames ',
         'for this vector include "q.value", "adj.P.Val", "FDR", "padj", and "FDR". ',
         'Make sure that dat includes exactly one such colname')
  }
  avg <- c('AvgMeth', 'AveExpr', 'logCPM', 'baseMean', 'AvgExpr')
  for (i in avg) {
    colnames(dat)[colnames(dat) == i] <- 'AvgExpr'
  }
  if (all(!avg %in% colnames(dat))) {
    stop('dat must include a column for average expression or methylation by ',
         'probe. Recognized colnames for this vector include "AvgExpr", "AvgMeth", ',
         '"AveExpr", "logCPM", and "baseMean". Make sure that dat includes exactly ',
         'one such colname')
  }
  if ('log2FoldChange' %in% colnames(dat)) {
    dat <- dat %>% rename(logFC = log2FoldChange)
  }
  if (!'log2FoldChange' %in% colnames(dat) &
      !'logFC' %in% colnames(dat)) {
    stop('dat must include a log fold change column. Recognized colnames for this ',
         'vector include "logFC" and "log2FoldChange". Make sure that dat includes ',
         'exactly one such colname')
  }
  if (is.null(main)) {
    main <- 'Taco Plot'
  }
  if (!legend %in% c('outside', 'bottomleft', 'bottomright', 'topleft', 'topright')) {
    stop('legend must be one of "outside", "bottomleft", "bottomright", ',
         '"topleft", or "topright"')
  }

  # Tidy
  test <- function(q) ifelse(q < fdr, paste('q <', fdr), paste('q >', fdr))
  df <- dat %>%
    mutate(Probe = probes,
           is.DE = map_chr(q.value, test),
           logP  = -log10(p.value)) %>%
    select(Probe, AvgExpr, logFC, logP, is.DE) %>%
    na.omit()
  if (sum(grepl('<', df$is.DE) == 0)) {
    warning('No probe meets your fdr threshold. To color data points by differential ',
            'expression/methylation, consider raising your fdr cutoff')
  }

  # Plot
  p <- plot_ly(df, x = ~AvgExpr, y = ~logFC, z = ~logP,
               text = ~Probe, color = ~is.DE, colors = c('red', 'black'),
               type = 'scatter3d', mode = 'markers',
               alpha = 0.85, hoverinfo = 'text', marker = list(size = 1)) %>%
    layout(hovermode = 'closest', title = main, scene = list(
      xaxis = list(title = 'Mean Expression'),  # expression() doesn't work here :(
      yaxis = list(title = 'log2 Fold Change'),
      zaxis = list(title = '-log10 p-value')))
  print(p)

}


