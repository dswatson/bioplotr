#' Mean-Difference Plot
#'
#' This function plots probewise means vs. log2 fold changes for a test of
#' differential expression or between-sample comparison.
#'
#' @param dat Either a data frame representing the results of a test for differential
#'   expression, or an omic data matrix with rows corresponding to probes and columns
#'   to samples. The former will render a study-wide MD plot, the latter a between
#'   between-sample MD plot. See Details.
#' @param fdr Threshold for declaring a probe differentially expressed. Only
#'   relevant for study-wide MD plots.
#' @param sample Column number or name specifying which sample in \code{dat} to
#'   compare with the others. Only relevant for between-sample MD plots.
#' @param ctrls Optional vector of length equal to \code{nrow(dat)} indicating the
#'   control status of each probe. Only relevant for between-sample MD plots.
#' @param ptsize Size of data points in the plot.
#' @param main Optional plot title.
#' @param legend Legend position. Must be one of \code{"outside",
#'   "bottomleft", "bottomright", "topleft",} or \code{"topright"}.
#' @param hover Show probe name by hovering mouse over data point? If \code{TRUE},
#'   the plot is rendered in HTML and will either open in your browser's graphic
#'   display or appear in the RStudio viewer. Probe names are extracted from
#'   \code{dat}.
#'
#' @details
#' MD plots (also known as "Bland-Altman plots" or "MA plots") visualize the
#' relationship between a probe's mean value and its log2 fold change versus some
#' relevant reference group. These figures help to evaluate the symmetry, magnitude,
#' and significance of differential effects across the full omic range.
#'
#' If \code{dat} is a data frame summarizing the results of a test for differential
#' expression, then each point's x-coordinate correponds to its average expression
#' across all samples, while y-coordinates represent the log2 fold change for the
#' given contrast. Points are colored to distinguish between those that do and do not
#' meet a user-defined FDR threshold. \code{plot_md} accepts output from \code{
#' limma::\link[limma]{topTable}}, \code{DESeq2::\link[DESeq2]{results}}, or \code{
#' edgeR::\link[edgeR]{topTags}}. Alternatively, any object with columns for log fold
#' changes, probewise means, and FDR is acceptable.
#'
#' If \code{dat} is probe by sample matrix or matrix-like object, then \code{sample}
#' must be specified. An artificial array is created by averaging probewise values for
#' all other samples in the data. The figure will then represent the mean vs. the
#' difference of expression values for the specified sample vs. the artificial array.
#'
#' @examples
#' library(limma)
#' DE_genes <- cbind(matrix(rnorm(50 * 5, mean = 5), nrow = 50, ncol = 5),
#'                   matrix(rnorm(50 * 5), nrow = 50, ncol = 5))
#' mat <- rbind(DE_genes, matrix(rnorm(4950 * 10), nrow = 4950, ncol = 10))
#' treat <- rep(c("A", "B"), each = 5)
#' des <- model.matrix(~ treat)
#' fit <- eBayes(lmFit(mat, des))
#' top <- topTable(fit, number = Inf)
#' plot_md(top)
#'
#' @seealso
#' \code{\link[limma]{plotMD}}, \code{\link[DESeq2]{plotMA}},
#' \code{\link[Glimma]{glMDPlot}}
#'
#' @export
#' @importFrom limma getEAWP
#' @importFrom purrr map_lgl map_chr
#' @importFrom ggsci pal_d3
#' @import dplyr
#' @import ggplot2
#' @importFrom plotly ggplotly
#'

plot_md <- function(dat,
                    fdr = 0.05,
                 sample = NULL,
                  ctrls = NULL,
                 ptsize = 0.25,
                   main = NULL,
                 legend = 'outside',
                  hover = FALSE) {

  # Preliminaries
  if ((is.data.frame(dat) || is(dat, 'DataFrame') || is(dat, 'TopTags'))) {
    dat <- as.data.frame(dat) %>% na.omit()
    if ('baseMean' %in% colnames(dat)) {
      dat$baseMean <- log2(dat$baseMean)
    }
    avg <- c('AveExpr', 'baseMean', 'logCPM', 'AvgExpr', 'AvgMeth')
    if (sum(avg %in% colnames(dat)) == 1L) {
      colnames(dat)[colnames(dat) %in% avg] <- 'Mean'
    } else {
      stop('dat must include a column for average expression or methylation by ',
           'probe. Recognized colnames for this vector include "AveExpr", "baseMean",
           "logCPM", "AvgExpr", and "AvgMeth".Make sure that dat includes exactly ',
           'one such colname.')
    }
    lfc <- c('logFC', 'log2FoldChange')
    if (sum(lfc %in% colnames(dat)) == 1L) {
      colnames(dat)[colnames(dat) %in% lfc] <- 'Diff'
    } else {
      stop('dat must include a log fold change column. Recognized colnames for this ',
           'vector include "logFC" and "log2FoldChange". Make sure that dat includes ',
           'exactly one such colname.')
    }
    q <- c('adj.P.Val', 'FDR', 'padj', 'q.value')
    if (sum(q %in% colnames(dat)) == 1L) {
      colnames(dat)[colnames(dat) %in% q] <- 'q.value'
    } else {
      stop('dat must include a column for adjusted p-values. Recognized colnames ',
           'for this vector include "q.value", "adj.P.Val", "FDR", "padj", and "FDR". ',
           'Make sure that dat includes exactly one such colname.')
    }
    if (!is.null(ctrls)) {
      warning('ctrls is ignored when dat is a data.frame or object of any class ',
              'inheriting from data.frame.')
    }
    if (!is.null(sample)) {
      warning('sample is ignored when dat is a data.frame or an object of any class ',
              'inheriting from data.frame.')
    }
  } else {
    if (!is.null(control)) {
      if (length(control) != nrow(dat)) {
        stop('control must be NULL or of length equal to nrow(dat).')
      }
    }
    if (is.null(sample)) {
      stop('When dat is an omic matrix, sample must be supplied.')
    } else {
      if (is.numeric(sample)) {
        if (!sample %in% seq_len(ncol(dat))) {
          stop(paste0('dat has no column number ', sample, '.'))
        }
      } else {
        if (!sample %in% colnames(dat)) {
          stop(paste0('dat has no column named "', sample, '".'))
        }
      }
    }
    if (!is.null(fdr)) {
      warning('fdr is ignored when dat is a matrix or matrix-like object.')
    }
    dat <- getEAWP(dat)$expr
    keep <- rowSums(is.finite(dat)) == ncol(dat)
    dat <- dat[keep, , drop = FALSE]
  }
  if (is.null(main)) main <- 'Mean-Difference Plot'
  if (!legend %in% c('outside', 'bottomleft', 'bottomright', 'topleft', 'topright')) {
    stop('legend must be one of "outside", "bottomleft", "bottomright", ',
         '"topleft", or "topright".')
  }

  # Tidy data
  if (!is.null(rownames(dat))) {
    probes <- rownames(dat)
  } else {
    probes <- seq_len(nrow(dat))
  }
  if (!is.matrix(dat)) {
    test <- function(q) ifelse(q < fdr, TRUE, FALSE)
    df <- dat %>%
      mutate(Probe = probes,
             is.DE = map_lgl(q.value, test)) %>%
      select(Probe, Mean, Diff, is.DE)
  } else {
    other <- rowMeans(dat[, -sample])
    df <- data_frame(Probe = probes,
                      Mean = (other + dat[, sample]) / 2L,
                      Diff = dat[, sample] - other)
    if (!is.null(ctrls)) {
      type <- function(p) {   ### THIS IS GONNA NEED WORK ###
        if (probe > 0L) 'Positive'
        else if (probe < 0L) 'Negative'
        else 'Zero'
      }
      df <- df %>% mutate(Ctrl = map_chr(ctrls, type),
                          Size = ifelse(ctrls == 0L, 1L, 2L))
    }
  }

  # Build plot
  suppressWarnings(
    p <- ggplot(df, aes(Mean, Diff, text = Probe)) +
      geom_hline(yintercept = 0L, size = 0.2, color = 'grey') +
      labs(title = main,
               x = expression(mu),
               y = expression('log'[2]*' Fold Change')) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
  )
  if ('is.DE' %in% colnames(df)) {
    if (sum(df$is.DE) == 0L) {                   # Color pts by differential expression?
      warning('No probe meets your fdr threshold. To color data points by differential ',
              'expression/methylation, consider raising your fdr cutoff.')
      p <- p + geom_point(size = ptsize, alpha = 0.25)
    } else {
      p <- p + geom_point(aes(color = is.DE), size = ptsize, alpha = 0.25) +
        scale_color_manual(name = expression(italic(q)*'-value'),
                         labels = c(paste('\u2265', fdr), paste('<', fdr)),
                         values = c('black', pal_d3()(4))) +
        guides(col = guide_legend(reverse = TRUE))
    }
  } else if ('Ctrl' %in% colnames(df)) {
    p <- p + geom_point(aes(color = Ctrl, size = Size), alpha = 0.25) +
      scale_size(range = c(ptsize, 5L * ptsize), guide = FALSE) +
      scale_color_manual(name = 'Control',
                       values = c(pal_d3()(seq_len(2)), 'black'))
  }
  if (legend == 'bottomleft') {                  # Locate legend
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
    if (legend == 'outside') {
      p <- ggplotly(p, tooltip = 'text', height = 525, width = 600)
    } else {
      p <- ggplotly(p, tooltip = 'text', height = 600, width = 600)
    }
  }

}


# Use gganimate, tweenr, and shiny to toggle between contrasts or samples
# Set FDR and/or fold change cutoffs on the fly

