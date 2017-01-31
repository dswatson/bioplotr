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
#'   identifiers, assuming they aren't \code{rownames(dat)}.
#'
#' @details
#' This function displays the results of a differential expression or methylation
#' test as a taco plot, which is essentially a combination of a volcano plot and an MD
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
                      fdr = 0.05,
                     main = NULL,
                   legend = 'outside',
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
  p <- c('P.Value', 'PValue', 'pvalue', 'p.value')
  if (any(p %in% colnames(dat))) {
    j <- intersect(p, colnames(dat))
    colnames(dat)[colnames(dat) == j] <- 'p.value'
  } else {
    stop('dat must include a p-value column. Recognized colnames for this vector ',
         'include "p.value", "P.Value", "PValue", and "pvalue". Make sure that dat',
         'includes exactly one such colname.')
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
    main <- 'Taco Plot'
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

  # Tidy data
  test <- function(q) ifelse(q < fdr, paste('q <', fdr), paste('q >', fdr))
  df <- dat %>%
    mutate(is.DE = map_chr(q.value, test),
            logP = -log10(p.value)) %>%
    select(Probe, AvgExpr, logFC, logP, is.DE) %>%
    na.omit()
  if (sum(grepl('<', df$is.DE) == 0)) {
    warning('No probe meets your fdr threshold. To color data points by differential ',
            'expression/methylation, consider raising your fdr cutoff.')
  }

  # Build Plot
  p <- plot_ly(df, x = ~AvgExpr, y = ~logFC, z = ~logP,
               text = ~Probe, color = ~is.DE, colors = c('red', 'black'),
               type = 'scatter3d', mode = 'markers',
               alpha = 0.85, hoverinfo = 'text', marker = list(size = 1)) %>%
    layout(hovermode = 'closest', title = main, scene = list(
      xaxis = list(title = 'Mean Expression'),
      yaxis = list(title = 'log2 Fold Change'),
      zaxis = list(title = '-log10 p-value')))
  print(p)

}



# Plotly doesn't play well with expression() text
# According to their documentation, LaTeX formating should work, e.g.
# xaxis = list(title = '$\\mu\\$')
# Except it doesn't cuz plotly sucks
