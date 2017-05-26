#' Volcano Plot
#'
#' This function plots log2 fold changes against -log10 \emph{p}-values for a
#' given test of differential expression.
#'
#' @param dat Data frame representing the results of a test for differential
#'   expression, such as the output of a call to \code{\link[limma]{topTable}},
#'   \code{\link[edgeR]{topTags}}, or \code{\link[DESeq2]{results}}.
#'   Alternatively, any object coercable to a data frame with columns for \emph{
#'   p}-values, log fold changes, and FDR. Missing values are silently removed.
#' @param fdr Significance threshold for declaring a probe differentially
#'   expressed.
#' @param lfc Optional effect size threshold for declaring a probe
#'   differentially expressed.
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"right"}, \code{
#'   "left"}, \code{"top"}, \code{"bottom"}, \code{"topright"}, \code{
#'   "topleft"}, \code{"bottomright"}, or \code{"bottomleft"}.
#' @param hover Show probe name by hovering mouse over data point? If \code{
#'   TRUE}, the plot is rendered in HTML and will either open in your browser's
#'   graphic display or appear in the RStudio viewer. Probe names are extracted
#'   from \code{dat}.
#'
#' @details
#' Volcano plots visualize the relationship between each probe's log2 fold
#' change and -log10 \emph{p}-value for a given test of differential expression.
#' Points are colored to distinguish between those that do and do not meet a
#' user-defined FDR threshold. Up- and down-regulated genes may also be
#' differentially colored if a minimum absolute fold change is supplied. These
#' figures help to evaluate the symmetry, magnitude, and significance of effects
#' in an omic experiment.
#'
#' @examples
#' # Simulated data
#' library(DESeq2)
#' dds <- makeExampleDESeqDataSet()
#' dds <- DESeq(dds)
#' res <- results(dds)
#' plot_volcano(res)
#'
#' # Real data
#' data(airway)
#' dds <- DESeqDataSet(airway, design = ~ cell + dex)
#' dds <- DESeq(dds)
#' res <- results(dds)
#' plot_volcano(res)
#'
#' @export
#' @import dplyr
#' @import ggplot2
#' @importFrom ggsci pal_d3
#'

plot_volcano <- function(dat,
                         fdr = 0.05,
                         lfc = NULL,
                       title = NULL,
                      legend = 'right',
                       hover = FALSE) {

  # Preliminaries
  dat <- as.data.frame(dat)
  fc <- c('log2FoldChange', 'logFC')
  if (sum(fc %in% colnames(dat)) == 1L) {        # Rename logFC
    colnames(dat)[colnames(dat) %in% fc] <- 'logFC'
  } else {
    stop(paste0('dat must include a log fold change column. Recognized ',
                'colnames for this vector include ', stringify(fc), '. ',
                'Make sure that dat includes exactly one such colname.'))
  }
  p <- c('P.Value', 'pvalue', 'PValue', 'p.value')
  if (sum(p %in% colnames(dat)) == 1L) {         # Rename p.value
    colnames(dat)[colnames(dat) %in% p] <- 'p.value'
  } else {
    stop(paste0('dat must include a p-value column. Recognized colnames for ',
                'this vector include ', stringify(p), '. Make sure that dat ',
                'includes exactly one such colname.'))
  }
  q <- c('adj.P.Val', 'padj', 'FDR', 'q.value')
  if (sum(q %in% colnames(dat)) == 1L) {         # Rename q.value
    colnames(dat)[colnames(dat) %in% q] <- 'q.value'
  } else {
    stop(paste0('dat must include a column for adjusted p-values. Recognized ',
                'colnames for this vector include ', stringify(q), '. Make ',
                'sure that dat includes exactly one such colname.'))
  }
  dat <- dat %>%
    select(logFC, p.value, q.value) %>%
    na.omit()
  if (nrow(dat) == 0L) {
    stop('dat must have at least one row with non-missing values for logFC, ',
         'p.value, and FDR.')
  }
  if (min(dat$p.value) < 0L | max(dat$p.value) > 1L) {
    stop('P-values must be on [0, 1].')
  }
  if (min(dat$q.value) < 0L | max(dat$q.value) > 1L) {
    stop('FDR values must be on [0, 1].')
  }
  if (is.null(title)) {
    title <- 'Volcano Plot'
  }
  if (!legend %in% c('right', 'left', 'top', 'bottom',
                     'topright', 'topleft', 'bottomright', 'bottomleft')) {
    stop('legend must be one of "right", "left", "top", "bottom", ',
         '"topright", "topleft", "bottomright", or "bottomleft".')
  }

  # Tidy data
  df <- dat %>%
    mutate(Probe = ifelse(is.null(rownames(dat)), row_number(), rownames(dat)),
            logP = -log10(p.value))
  if (!is.null(lfc)) {
    df <- df %>%
      mutate(Direction = ifelse(q.value <= fdr & logFC >= lfc, 'Up',
                                ifelse(q.value <= fdr & -logFC >= lfc, 'Down', 'NA')))
  }

  # Build plot
  size <- pt_size(df)
  alpha <- pt_alpha(df)
  p <- ggplot(df, aes(logFC, logP, text = Probe)) +
    labs(title = title,
             x = expression(log[2]~'Fold Change'),
             y = expression(~-log[10](italic(p)))) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (!any(df$q.value <= fdr)) {                 # Color pts by differential expression?
    warning('No probe meets your fdr threshold. To color data points by',
            'differential expression, consider raising your fdr cutoff.')
    p <- p + geom_point(size = size, alpha = alpha, color = 'black')
  } else {                                       # Separate up- and down-regulated probes?
    if (!is.null(lfc) & any(df$Direction != 'NA')) {
      y <- df %>%
        filter(q.value <= fdr) %>%
        filter(logP == min(logP)) %>%
        select(logP) %>%
        as.numeric()
      suppressWarnings(
        p <- p + geom_point(data = df %>% filter(Direction != 'NA'),
                            aes(logFC, logP, color = Direction, text = Probe),
                            size = size, alpha = alpha) +
          geom_point(data = df %>% filter(Direction == 'NA'),
                     aes(logFC, logP, text = Probe),
                     color = '#444444', size = size, alpha = alpha) +
          geom_segment(x = -Inf, xend = -lfc, y = y, yend = y, linetype = 2L) +
          geom_segment(x = lfc, xend = Inf, y = y, yend = y, linetype = 2L) +
          geom_segment(x = -lfc, xend = -lfc, y = y, yend = Inf, linetype = 2L) +
          geom_segment(x = lfc, xend = lfc, y = y, yend = Inf, linetype = 2L) +
          annotate('text', x = min(df$logFC), y, label = paste('FDR =', fdr),
                   hjust = 0L, vjust = -1L) +
          scale_color_manual(guide = FALSE, values = pal_d3()(4)[3:4])
      )
    } else {
      p <- p + geom_point(aes(color = q.value <= fdr), size = size, alpha = alpha) +
        scale_color_manual(name = 'FDR',
                         labels = c(paste('>', fdr), paste('\u2264', fdr)),
                         values = c('black', pal_d3()(4)[4])) +
        guides(color = guide_legend(reverse = TRUE))
      if (!is.null(lfc) & all(df$Direction == 'NA')) {
        warning('No probe meets both your fdr and lfc criteria. To color',
                'probes by the direction of their differential expression,',
                'consider lowering your lfc threshold.')
      }
    }
  }

  # Output
  gg_out(p, hover, legend)

}


