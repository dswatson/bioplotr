#' Taco Plot
#'
#' This function plots effect size against significance against probewise means
#' for a given test of differential expression.
#'
#' @param dat Data frame representing the results of a test for differential
#'   expression, such as the output of a call to \code{
#'   limma::\link[limma]{topTable}}, \code{edgeR::\link[edgeR]{topTags}}, or
#'   \code{DESeq2::\link[DESeq2]{results}}. Alternatively, any object coercable
#'   to a data frame with columns for log fold changes, probewise means,
#'   \emph{p}-values, and FDR. Missing values are silently removed.
#' @param fdr Significance threshold for declaring a probe differentially
#'   expressed.
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"bottom"}, \code{"left"},
#'   \code{"top"}, \code{"right"}, \code{"bottomright"}, \code{"bottomleft"},
#'   \code{"topleft"}, or \code{"topright"}.
#'
#' @details
#' A taco plot combines the elements of a volcano plot and an MD plot into a
#' single three-dimensional figure. Points are colored to distinguish between
#' those that do and do not meet a user-defined FDR threshold. These
#' figures help to evaluate the symmetry, magnitude, and significance of effects
#' in an omic experiment.
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
#' @import dplyr
#' @importFrom ggsci pal_d3
#' @importFrom purrr map_chr
#'

plot_taco <- function(
  dat,
     fdr = 0.05,
   title = 'Taco Plot',
  legend = 'right'
) {

  # Preliminaries
  dat <- dat %>%
    as.data.frame(.) %>%
    na.omit(.)
  lfc <- c('logFC', 'log2FoldChange')
  if (sum(lfc %in% colnames(dat)) == 1L) {
    colnames(dat)[colnames(dat) %in% lfc] <- 'logFC'
  } else {
    stop('dat must include a log fold change column. Recognized colnames for ',
         'this vector include ', stringify(fc, 'and'), '. Make sure that dat ',
         'includes exactly one such colname.')
  }
  p <- c('P.Value', 'PValue', 'pvalue', 'p.value')
  if (sum(p %in% colnames(dat)) == 1L) {
    colnames(dat)[colnames(dat) %in% p] <- 'p.value'
  } else {
    stop('dat must include a p-value column. Recognized colnames for this ',
         'vector include ', stringify(p, 'and'), '. Make sure that dat ',
         'includes exactly one such colname.')
  }
  if ('baseMean' %in% colnames(dat)) {
    dat$baseMean <- log2(dat$baseMean)
  }
  avg <- c('AveExpr', 'baseMean', 'logCPM', 'AvgExpr', 'AvgMeth')
  if (sum(avg %in% colnames(dat)) == 1L) {
    colnames(dat)[colnames(dat) %in% avg] <- 'AvgExpr'
  } else {
    stop('dat must include a column for average expression by probe. ',
         'Recognized colnames for this vector include ', stringify(avg, 'and'),
         '. Make sure that dat includes exactly one such colname.')
  }
  q <- c('adj.P.Val', 'padj', 'FDR', 'q.value')
  if (sum(q %in% colnames(dat)) == 1L) {
    colnames(dat)[colnames(dat) %in% q] <- 'q.value'
  } else {
    stop('dat must include a column for adjusted p-values. Recognized ',
         'colnames for this vector include ', stringify(q), '. Make sure that ',
         'dat includes exactly one such colname.')
  }
  locations <- c('bottom', 'left', 'top', 'right',
                 'bottomright', 'bottomleft', 'topleft', 'topright')
  legend <- match.arg(legend, locations)

  # Tidy data
  if (rownames(dat) %>% is.null) {
    dat <- dat %>% mutate(Probe = row_number())
  } else {
    dat <- dat %>% mutate(Probe = rownames(dat))
  }
  test <- function(q) ifelse(q < fdr, paste('q <', fdr), paste('q >', fdr))
  df <- dat %>%
    mutate(is.DE = map_chr(q.value, test),
            logP = -log10(p.value)) %>%
    select(Probe, AvgExpr, logFC, logP, is.DE)
  if (sum(grepl('<', df$is.DE) == 0L)) {
    warning('No probe meets your fdr threshold. To color data points by ',
            'differential expression, consider raising your fdr cutoff.')
  }

  # Build Plot
  require(plotly)
  p <- plot_ly(df, x = ~AvgExpr, y = ~logFC, z = ~logP,
               text = ~Probe, color = ~is.DE, 
               colors = c(pal_d3()(4L)[4L], 'black'),
               type = 'scatter3d', mode = 'markers',
               alpha = 0.85, hoverinfo = 'text', marker = list(size = 1)) %>%
    layout(hovermode = 'closest', title = title, scene = list(
      xaxis = list(title = 'Mean Expression'),
      yaxis = list(title = 'log2 Fold Change'),
      zaxis = list(title = '-log10 p-value')))
  print(p)

}

# Plotly doesn't play well with expression() text
# Can't handle plotmath, LaTeX typesetting, etc. in R:
# https://github.com/ropensci/plotly/issues/375
# Altho it used to? https://plot.ly/r/LaTeX/

# Add lfc threshold?

# More elegant solution to probe column so we can get straight to tibble
