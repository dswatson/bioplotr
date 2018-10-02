#' Volcano Plot
#'
#' This function plots effect size against significance for a given test of
#' differential expression.
#'
#' @param dat Data frame or similar object representing the results of a test
#'   for differential expression, such as the output of a call to \code{
#'   limma::\link[limma]{topTable}}, \code{edgeR::\link[edgeR]{topTags}}, or
#'   \code{DESeq2::\link[DESeq2]{results}}. Alternatively, any object coercable
#'   to a data frame with columns for \emph{p}-values, \emph{q}-values, and log
#'   fold changes. Missing values are silently removed.
#' @param y Plot \emph{p}-values (\code{y = "p"}) or \emph{q}-values (\code{
#'   y = "q"}) on the y-axis?
#' @param fdr Significance threshold for declaring a probe differentially
#'   expressed.
#' @param lfc Optional effect size threshold for declaring a probe
#'   differentially expressed.
#' @param probe Name of column in which to find probe IDs, presuming they are
#'   not stored in \code{rownames(dat)}.
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"bottom"}, \code{"left"},
#'   \code{"top"}, \code{"right"}, \code{"bottomright"}, \code{"bottomleft"},
#'   \code{"topleft"}, or \code{"topright"}.
#' @param hover Show probe name by hovering mouse over data point? If \code{
#'   TRUE}, the plot is rendered in HTML and will either open in your browser's
#'   graphic display or appear in the RStudio viewer. Probe names are extracted
#'   from \code{dat}.
#'
#' @details
#' Volcano plots visualize the relationship between each probe's log2 fold
#' change and -log10 \emph{p}- or \emph{q}-value for a given test of
#' differential expression. Points are colored to distinguish between those that
#' do and do not meet a user-defined FDR threshold. Up- and down-regulated genes
#' may also be differentially colored if a minimum absolute fold change is
#' supplied. These figures help to evaluate the symmetry, magnitude, and
#' significance of effects in an omic experiment.
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
                         y = 'q',
                       fdr = 0.05,
                       lfc = NULL,
                     probe = NULL,
                     title = NULL,
                    legend = 'right',
                     hover = FALSE) {

  # Preliminaries
  dat <- as.data.frame(dat)
  fc <- c('log2FoldChange', 'logFC')
  if (sum(fc %in% colnames(dat)) == 1L) {        # Rename logFC
    colnames(dat)[colnames(dat) %in% fc] <- 'logFC'
  } else {
    stop('dat must include a log fold change column. Recognized colnames for ',
         'this vector include ', stringify(fc, 'and'), '. Make sure that dat ',
         'includes exactly one such colname.')
  }
  p <- c('P.Value', 'pvalue', 'PValue', 'p.value')
  if (sum(p %in% colnames(dat)) == 1L) {         # Rename p.value
    colnames(dat)[colnames(dat) %in% p] <- 'p.value'
  } else {
    stop('dat must include a p-value column. Recognized colnames for this ',
         'vector include ', stringify(p, 'and'), '. Make sure that dat ',
         'includes exactly one such colname.')
  }
  q <- c('adj.P.Val', 'padj', 'FDR', 'qvalue', 'q.value')
  if (sum(q %in% colnames(dat)) == 1L) {         # Rename q.value
    colnames(dat)[colnames(dat) %in% q] <- 'q.value'
  } else {
    stop('dat must include a column for adjusted p-values. Recognized ',
         'colnames for this vector include ', stringify(q, 'and'), '. ',
         'Make sure that dat includes exactly one such colname.')
  }
  if (!y %in% c('p', 'q')) {
    stop('y must be one of "p" or "q".')
  }
  if (!(probe %>% is.null)) {
    if (!probe %in% colnames(dat)) {
      stop('Column "', probe, '" not detected.')
    } else {
      colnames(dat)[which(colnames(dat)) == probe] <- 'Probe'
    }
  } else {
    dat <- dat %>% mutate(Probe = rownames(.))
  }
  dat <- dat %>%
    select(Probe, logFC, p.value, q.value) %>%
    na.omit(.)
  if (nrow(dat) == 0L) {
    stop('dat must have at least one row with non-missing values for logFC, ',
         'p.value, and q.value.')
  }
  if (min(dat$p.value) < 0L || max(dat$p.value) > 1L) {
    stop('P-values must be on [0, 1].')
  }
  if (min(dat$q.value) < 0L || max(dat$q.value) > 1L) {
    stop('Q-values values must be on [0, 1].')
  }
  if (title %>% is.null) {
    title <- 'Volcano Plot'
  }
  loc <- c('bottom', 'left', 'top', 'right',
           'bottomright', 'bottomleft', 'topleft', 'topright')
  if (!legend %in% loc) {
    stop('legend must be one of ', stringify(loc, 'or'), '.')
  }

  # Tidy data
  if (y == 'p') {
    df <- dat %>% mutate(Y = -log10(p.value))
  } else {
    df <- dat %>% mutate(Y = -log10(q.value))
  }
  if (!(lfc %>% is.null)) {
    df <- df %>%
      mutate(Direction = ifelse(q.value <= fdr & logFC >= lfc, 'Up',
                                ifelse(q.value <= fdr & -logFC >= lfc, 'Down',
                                       'None')))
  }

  # Build plot
  size <- pt_size(df)
  alpha <- pt_alpha(df)
  p <- ggplot(df, aes(logFC, Y, text = Probe)) +
    labs(title = title, x = expression(log[2]~'Fold Change')) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (y == 'p') {
    p <- p + ylab(expression(~-log[10](italic(p))))
    y_pt <- df %>%
      filter(q.value <= fdr) %>%
      filter(Y == min(Y)) %>%
      select(Y) %>%
      as.numeric(.)
  } else {
    p <- p + ylab(expression(~-log[10](italic(q))))
    y_pt <- -log10(fdr)
  }
  if (!any(df$q.value <= fdr)) {                 # Color pts by differential expression?
    warning('No probe meets your fdr threshold. To color data points by ',
            'differential expression, consider raising your fdr cutoff.')
    p <- p + geom_point(size = size, alpha = alpha, color = 'black')
  } else {                                       # Separate up- and down-regulated probes?
    if (!(lfc %>% is.null) && any(df$Direction != 'None')) {
      suppressWarnings(
        p <- p + geom_point(data = filter(df, Direction != 'None'),
                            aes(logFC, Y, color = Direction, text = Probe),
                            size = size, alpha = alpha) +
          geom_point(data = filter(df, Direction == 'None'),
                     aes(logFC, Y, text = Probe),
                     color = '#444444', size = size, alpha = alpha) +
          geom_segment(x = -Inf, xend = -lfc, y = y_pt, yend = y_pt, linetype = 2L) +
          geom_segment(x = lfc, xend = Inf, y = y_pt, yend = y_pt, linetype = 2L) +
          geom_segment(x = -lfc, xend = -lfc, y = y_pt, yend = Inf, linetype = 2L) +
          geom_segment(x = lfc, xend = lfc, y = y_pt, yend = Inf, linetype = 2L) +
          # annotate('text', x = min(df$logFC), y, label = paste('FDR =', fdr),
          #          hjust = 0L, vjust = -1L) +
          scale_color_manual(guide = FALSE, values = pal_d3()(4L)[3:4])
      )
    } else {
      p <- p + geom_point(aes(color = q.value <= fdr),
                          size = size, alpha = alpha) +
        scale_color_manual(name = expression(italic(q)*'-value'),
                         labels = c(paste('>', fdr), paste('\u2264', fdr)),
                         values = c('#444444', pal_d3()(4L)[4L])) +
        guides(color = guide_legend(reverse = TRUE))
      if (!(lfc %>% is.null) && all(df$Direction == 'None')) {
        warning('No probe meets both your fdr and lfc criteria. To color ',
                'probes by the direction of their differential expression, ',
                'consider lowering your lfc threshold.')
      }
    }
  }

  # Output
  gg_out(p, hover, legend)

}


