#' Taco Plot
#'
#' This function plots log2 fold changes against -log10 \emph{p}-values against
#' probewise means for a given test of differential expression/methylation.
#'
#' @param dat Data frame or matrix representing the results of a test for
#'   differential expression or methylation, such as the output of a call to
#'   \code{\link[limma]{topTable}}, \code{\link[edgeR]{topTags}}, or
#'   \code{\link[DESeq2]{results}}. Alternatively, any object with columns for log
#'   fold changes, probewise means, \emph{p}-values, and FDR.
#' @param fdr Threshold for declaring a probe differentially expressed/methylated.
#' @param main Optional plot title.
#' @param legend Legend position. Must be one of \code{"outside",
#'   "bottomleft", "bottomright", "topleft",} or \code{"topright"}.
#' @param probes String specifying the name of the column in which to find the probe
#'   identifiers, assuming they aren't \code{rownames(dat)}.
#'
#' @details
#' A taco plot combines the elements of a volcano plot and an MD plot into a single
#' three-dimensional figure. Points are colored to distinguish between those that
#' do and do not meet a user-defined FDR threshold. These figures help to evaluate
#' the symmetry, magnitude, and significance of effects for a given experiment.
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
  if (sum(lfc %in% colnames(dat)) == 1) {
    colnames(dat)[colnames(dat) %in% lfc] <- 'logFC'
  } else {
    stop('dat must include a log fold change column. Recognized colnames for this ',
         'vector include "logFC" and "log2FoldChange". Make sure that dat includes ',
         'exactly one such colname.')
  }
  avg <- c('AvgMeth', 'AveExpr', 'logCPM', 'baseMean', 'AvgExpr')
  if (sum(avg %in% colnames(dat)) == 1) {
    colnames(dat)[colnames(dat) %in% avg] <- 'AvgExpr'
  } else {
    stop('dat must include a column for average expression or methylation by ',
         'probe. Recognized colnames for this vector include "AvgExpr", "AvgMeth", ',
         '"AveExpr", "logCPM", and "baseMean". Make sure that dat includes exactly ',
         'one such colname.')
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
    main <- 'Taco Plot'
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
