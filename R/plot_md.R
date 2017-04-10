#' Mean-Difference Plot
#'
#' This function plots probewise means vs. log2 fold changes for a test of
#' differential expression or between-sample comparison.
#'
#' @param dat Either a data frame representing the results of a test for differential
#'   expression, or an omic data matrix or matrix-like object with rows corresponding
#'   to probes and columns to samples. The former will render a study-wide MD plot,
#'   the latter a between between-sample MD plot. See Details.
#' @param design Optional design matrix with rows corresponding to samples and columns
#'   to coefficients to be estimated. Only relevant for \code{\link[edgeR]{DGEList}}
#'   objects. See Details.
#' @param fdr Significance threshold for declaring a probe differentially expressed. Only
#'   relevant for study-wide MD plots.
#' @param lfc Optional effect size threshold for declaring a probe differentially
#'   expressed. Only relevant for study-wide MD plots.
#' @param sample Column number or name specifying which sample in \code{dat} to
#'   compare with the others. Only relevant for between-sample MD plots.
#' @param ctrls Optional vector of length equal to \code{nrow(dat)} indicating the
#'   control status of each probe. Only relevant for between-sample MD plots.
#' @param title Optional plot title.
#' @param xlab Optional label for x-axis.
#' @param legend Legend position. Must be one of \code{"outside", "bottomleft",
#'   "bottomright", "topleft",} or \code{"topright"}.
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
#' If \code{dat} summarizes the results of a test for differential expression, then
#' each point's x-coordinate correponds to its average expression across all samples,
#' while y-coordinates represent the log2 fold change for the given contrast. Points
#' are colored to distinguish between those that do and do not meet a user-defined
#' FDR threshold. \code{plot_md} accepts output from
#' \code{limma::\link[limma]{topTable}}, \code{edgeR::\link[edgeR]{topTags}}, or
#' \code{DESeq2::\link[DESeq2]{results}}. Alternatively, any object with columns for
#' log fold changes, probewise means, and FDR is acceptable.
#'
#' If \code{dat} is probe by sample matrix or matrix-like object, then \code{sample}
#' must be specified. An artificial array is created by averaging probewise values for
#' all other samples in the data. The figure will then represent the mean vs. the
#' difference of expression values for the specified sample vs. the artificial array.
#' Acceptable inputs for between-sample MD plots include all \code{limma} expression
#' set objects, as well as \code{\link[edgeR]{DGEList}}, \code{\link[DESeq2]{
#' DESeqDataSet}}, and \code{\link[DESeq2]{DESeqTransform}} objects.
#'
#' @references
#' Bolstad, B.M., Irizarry, R.A., Åstrand, M. & Speed, T.P. (2003).
#' \href{https://www.ncbi.nlm.nih.gov/pubmed/12538238}{A comparison of normalization
#' methods for high density oligonucleotide array data based on variance and bias}.
#' \emph{Bioinformatics}, \emph{19}(2): 185–193.
#'
#' Dudoit, S., Yang, Y.H., Callow, M.J. & Speed, T.P. (2002).
#' \href{https://pdfs.semanticscholar.org/2af3/26eabfc0e81d6f3e687e1283c32cfba25688.pdf}{
#' Statistical methods for identifying differentially expressed genes in replicated cDNA
#' microarray experiments}. \emph{Stat. Sin.}, \strong{12}, 111–140.
#'
#' Martin, B.J. & Altman, D.G. (1986).
#' \href{http://www.thelancet.com/journals/lancet/article/PIIS0140-6736(86)90837-8/abstract}{
#' Statistical methods for assessing agreement between two methods of clinical measurement}.
#' \emph{Lancet}, 327: 307–310.
#'
#' Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., & Smyth, G.K. (2015).
#' \href{https://www.ncbi.nlm.nih.gov/pubmed/25605792}{limma powers differential
#' expression analyses for RNA-sequencing and microarray studies}. \emph{Nucleic
#' Acids Res.}, \emph{43}(7): e47.
#'
#' @examples
#' library(DESeq2)
#' dds <- makeExampleDESeqDataSet()
#'
#' # Between-sample MD plot
#' rld <- rlog(dds)
#' plot_md(rld)
#'
#' # Study-wide MD plot
#' dds <- DESeq(dds)
#' res <- results(dds)
#' plot_md(res)
#'
#' @seealso
#' \code{\link[limma]{plotMD}, \link[DESeq2]{plotMA}, \link[Glimma]{glMDPlot}}
#'
#' @export
#' @import dplyr
#' @import ggplot2
#' @importFrom plotly ggplotly
#'

plot_md <- function(dat,
                    title = NULL,
                   legend = 'outside', ...) {

  # Preliminaries
  if (is.null(title)) {
    title <- 'Mean-Difference Plot'
  }
  if (!legend %in% c('outside', 'bottomleft', 'bottomright', 'topleft', 'topright')) {
    stop('legend must be one of "outside", "bottomleft", "bottomright", ',
         '"topleft", or "topright".')
  }

  # Method
  UseMethod('plot_md')

}


#' @rdname plot_md
#' @method plot_md DGEList
#' @S3method plot_md DGEList
#' @importFrom edgeR calcNormFactors aveLogCPM estimateDisp cpm
#' @importFrom ggsci scale_color_d3

plot_md.DGEList <- function(dat,
                            design = NULL,
                            sample = 1,
                             ctrls = NULL,
                               lfc = NULL,
                             title = NULL,
                              xlab = NULL,
                            legend = 'outside',
                             hover = FALSE) {

  # Preliminaries
  if (is.numeric(sample) & sample > ncol(dat)) {
    stop('Sample number exceeds ncol(dat).')
  }
  if (is.character(sample)) {
    if (!sample %in% colnames(dat)) {
      stop(paste0('Could not detect a sample named "', sample, '" in dat.'))
    } else {
      sample <- which(colnames(dat) == sample)
    }
  }
  if (!is.null(ctrls)) {
    if (length(ctrls) != nrow(dat)) {
      stop('ctrls must be NULL or of length equal to nrow(dat).')
    }
    ctrls <- as.character(ctrls)
    ctrls[ctrls == names(table(ctrls)[which.max(table(ctrls))])] <- '0'
  }
  if (is.null(xlab)) {
    xlab <- expression('Mean'~log[2]*'-CPM')
  }

  # Tidy data
  keep <- rowSums(dat$counts) > 0L               # Minimal count filter
  dat <- dat[keep, ]
  other <- dat[, -sample]
  other <- calcNormFactors(other)
  if (is.null(other$tagwise.dispersion)) {       # Estimate dispersions for aveLogCPM
    if (is.null(dat$design) & !is.null(dat$group)) {
      other$design <- model.matrix(~ dat$group[-sample, ])
    }
    if (!is.null(design)) {
      other$design <- design[-sample, ]
    }
    if (is.null(other$design)) {
      other <- estimateTagwiseDisp(other, dispersion = estimateCommonDisp(other))
    } else {
      other <- estimateDisp(other)
    }
  }
  other <- aveLogCPM(other, prior.count = 1L, dispersion = other$tagwise.dispersion)
  dat <- cpm(dat, log = TRUE, prior.count = 1L)
  df <- data_frame(Probe = rownames(dat),
                    Mean = (other + dat[, sample]) / 2L,
                    Diff = dat[, sample] - other)
  if (!is.null(ctrls)) {
    ctrls <- ctrls[keep]
    df <- df %>% mutate(Control = ctrls)
  }

  # Build plot
  size <- probe_ptsize(df)
  alpha <- probe_alpha(df)
  p <- ggplot(df, aes(Mean, Diff, text = Probe)) +
    geom_hline(yintercept = 0L, color = 'grey') +
    labs(title = title, x = xlab, y = expression(log[2]~'Fold Change')) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (is.null(ctrls)) {
    p <- p + geom_point(size = size, alpha = alpha)
  } else {
    suppressWarnings(
      p <- p + geom_point(data = df %>% filter(Control == '0'),
                          aes(Mean, Diff, text = Probe),
                          size = size, alpha = alpha) +
        geom_point(data = df %>% filter(Control != '0'),
                   aes(Mean, Diff, color = Control, text = Probe),
                   size = 3 * size, alpha = alpha) +
        scale_color_d3()
    )
  }
  if (!is.null(lfc)) {
    p <- p + geom_hline(yintercept = lfc, linetype = 2L) +
      geom_hline(yintercept = -lfc, linetype = 2L)
  }
  p <- locate_legend(p, legend)

  # Output
  gg_out(p, hover, legend)

}


#' @rdname plot_md
#' @method plot_md DESeqDataSet
#' @S3method plot_md DESeqDataSet
#' @importFrom DESeq2 counts sizeFactors normalizationFactors estimateSizeFactors
#'   estimateDispersions dispersions
#' @importFrom edgeR aveLogCPM cpm
#' @importFrom ggsci scale_color_d3

plot_md.DESeqDataSet <- function(dat,
                                 sample = 1,
                                  ctrls = NULL,
                                    lfc = NULL,
                                  title = NULL,
                                   xlab = NULL,
                                 legend = 'outside',
                                  hover = FALSE) {

  # Preliminaries
  if (is.numeric(sample) & sample > ncol(dat)) {
    stop('Sample number exceeds ncol(dat).')
  }
  if (is.character(sample)) {
    if (!sample %in% colnames(dat)) {
      stop(paste0('Could not detect a sample named "', sample, '" in dat.'))
    } else {
      sample <- which(colnames(dat) == sample)
    }
  }
  if (!is.null(ctrls)) {
    if (length(ctrls) != nrow(dat)) {
      stop('ctrls must be NULL or of length equal to nrow(dat).')
    }
    ctrls <- as.character(ctrls)
    ctrls[ctrls == names(table(ctrls)[which.max(table(ctrls))])] <- '0'
  }
  if (is.null(xlab)) {
    xlab <- expression('Mean'~log[2]*'-CPM')
  }
  if (is.null(sizeFactors(dat)) & is.null(normalizationFactors(dat))) {
    dat <- estimateSizeFactors(dat)
  }
  if (is.null(dispersions(dat))) {
    dat <- estimateDispersions(dat, quiet = TRUE)
  }

  # Tidy data
  keep <- mcols(dat)$baseMean > 0L
  dat <- dat[keep, , drop = FALSE]
  cnts <- counts(dat, normalized = TRUE)
  other <- aveLogCPM(cnts[, -sample], prior.count = 1L,
                     dispersion = dispersions(dat))
  dat <- cpm(cnts, log = TRUE, prior.count = 1L)
  df <- data_frame(Probe = rownames(dat),
                    Mean = (other + dat[, sample]) / 2L,
                    Diff = dat[, sample] - other)
  if (!is.null(ctrls)) {
    ctrls <- ctrls[keep]
    df <- df %>% mutate(Control = ctrls)
  }

  # Build plot
  size <- probe_ptsize(df)
  alpha <- probe_alpha(df)
  p <- ggplot(df, aes(Mean, Diff, text = Probe)) +
    geom_hline(yintercept = 0L, color = 'grey') +
    labs(title = title, x = xlab, y = expression(log[2]~'Fold Change')) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (is.null(ctrls)) {
    p <- p + geom_point(size = size, alpha = alpha)
  } else {
    suppressWarnings(
      p <- p + geom_point(data = df %>% filter(Control == '0'),
                          aes(Mean, Diff, text = Probe),
                          size = size, alpha = alpha) +
        geom_point(data = df %>% filter(Control != '0'),
                   aes(Mean, Diff, color = Control, text = Probe),
                   size = 3 * size, alpha = alpha) +
        scale_color_d3()
    )
  }
  if (!is.null(lfc)) {
    p <- p + geom_hline(yintercept = lfc, linetype = 2L) +
      geom_hline(yintercept = -lfc, linetype = 2L)
  }
  p <- locate_legend(p, legend)

  # Output
  gg_out(p, hover, legend)

}


#' @rdname plot_md
#' @method plot_md DESeqTransform
#' @S3method plot_md DESeqTransform
#' @importFrom SummarizedExperiment assay

plot_md.DESeqTransform <- function(dat,
                                   sample = 1,
                                    ctrls = NULL,
                                      lfc = NULL,
                                    title = NULL,
                                     xlab = NULL,
                                   legend = 'outside',
                                    hover = FALSE) {

  # Preliminaries
  if (is.null(xlab)) {
    xlab <- expression('Mean Transformed Counts')
  }

  # Tidy data
  dat <- assay(dat)

  # Export
  plot_md.default(dat = dat, sample = sample, ctrls = ctrls, lfc = lfc,
                  title = title, xlab = xlab, legend = legend, hover = hover)

}


#' @rdname plot_md
#' @method plot_md DESeqResults
#' @S3method plot_md DESeqResults

plot_md.DESeqResults <- function(dat,
                                 fdr = 0.05,
                                 lfc = NULL,
                               title = NULL,
                                xlab = NULL,
                              legend = 'outside',
                               hover = FALSE) {

  # Preliminaries
  if (is.null(xlab)) {
    xlab <- 'Mean of Normalized Counts'
  }
  dat <- na.omit(dat)
  if (nrow(dat) == 0L) {
    stop('dat must have at least one row with non-missing values for baseMean, ',
         'log2FoldChange, and padj.')
  }

  # Tidy data
  df <- as.data.frame(dat) %>%
    na.omit() %>%
    mutate(Probe = rownames(dat)) %>%
    rename(Mean = baseMean,
           Diff = log2FoldChange,
        q.value = padj) %>%
    select(Probe, Mean, Diff, q.value)
  if (!is.null(lfc)) {
    df <- df %>%
      mutate(Direction = ifelse(q.value <= fdr & Diff >= lfc, 'Up',
                                ifelse(q.value <= fdr & -Diff >= lfc, 'Down', 'NA')))
  }

  # Build plot
  size <- probe_ptsize(df)
  alpha <- probe_alpha(df)
  p <- ggplot(df, aes(Mean, Diff, text = Probe)) +
    scale_x_log10() +
    geom_hline(yintercept = 0L, color = 'grey') +
    labs(title = title, x = xlab, y = expression(log[2]~'Fold Change')) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (!any(df$q.value <= fdr)) {               # Color pts by differential expression?
    warning('No probe meets your fdr threshold. To color data points by differential ',
            'expression, consider raising your fdr cutoff.')
    p <- p + geom_point(size = size, alpha = alpha)
  } else {
    if (is.null(lfc)) {
      p <- p + geom_point(aes(color = q.value <= fdr), size = size, alpha = alpha) +
        scale_color_manual(name = 'FDR',
                         labels = c(paste('>', fdr), paste('\u2264', fdr)),
                         values = c('#444444', pal_d3()(4)[4]),
                          guide = guide_legend(reverse = TRUE, override.aes = list(
                           size = rep(1L, 2L), alpha = rep(1L, 2L))))
    } else {
      suppressWarnings(
        p <- p + geom_hline(yintercept = lfc, linetype = 2L) +
          geom_hline(yintercept = -lfc, linetype = 2L) +
          geom_point(data = df %>% filter(Direction != 'NA'),
                     aes(Mean, Diff, color = Direction, text = Probe),
                     size = size, alpha = alpha) +
          geom_point(data = df %>% filter(Direction == 'NA'),
                     aes(Mean, Diff, text = Probe),
                     color = '#444444', size = size, alpha = alpha) +
          scale_color_manual(guide = FALSE, values = pal_d3()(4)[3:4])
      )
    }
  }
  p <- locate_legend(p, legend)

  # Output
  gg_out(p, hover, legend)

}


#' @rdname plot_md
#' @method plot_md TopTags
#' @S3method plot_md TopTags

plot_md.TopTags <- function(dat,
                            fdr = 0.05,
                            lfc = NULL,
                          title = NULL,
                           xlab = NULL,
                         legend = 'outside',
                          hover = FALSE) {

  # Preliminaries
  if (is.null(xlab)) {
    xlab <- expression('Mean'~log[2]*'-CPM')
  }

  # Export
  dat <- as.data.frame(dat)
  plot_md.data.frame(dat = dat, fdr = fdr, lfc = lfc,
                     title = title, xlab = xlab, legend = legend, hover = hover)

}


#' @rdname plot_md
#' @method plot_md data.frame
#' @S3method plot_md data.frame

plot_md.data.frame <- function(dat,
                               probes = NULL,
                                  fdr = 0.05,
                                  lfc = NULL,
                                title = NULL,
                                 xlab = NULL,
                               legend = 'outside',
                                hover = FALSE) {

  # Preliminaries
  if (is.null(probes)) {
    if (is.null(rownames(dat)) |
        identical(rownames(dat), as.character(seq_len(nrow(dat))))) {
      stop('If dat does not have rownames, then the column of unique probe ',
           'identifiers must be specified using the probes argument.')
    } else {
      dat$Probe <- rownames(dat)
    }
  } else {
    if (is.numeric(probes)) {
      if (probes > ncol(dat)) {
        stop('Column number for probes exceeds ncol(dat).')
      } else {
        colnames(dat)[probes] <- 'Probe'
      }
    } else {
      if (!probes %in% colnames(dat)) {
        stop(paste0('Could not detect a column named "', probes, '" in dat.'))
      } else {
        colnames(dat)[colnames(dat) == probes] <- 'Probe'
      }
    }
  }
  if ('baseMean' %in% colnames(dat)) {
    dat$baseMean <- log2(dat$baseMean)
  }
  avg <- c('AveExpr', 'baseMean', 'logCPM', 'AvgExpr', 'AvgMeth')
  if (sum(avg %in% colnames(dat)) == 1L) {       # Rename AvgExpr
    colnames(dat)[colnames(dat) %in% avg] <- 'Mean'
  } else {
    stop('dat must include a column for average expression by probe. Recognized ',
         'colnames for this vector include "AveExpr", "baseMean", and "logCPM". ',
         'Make sure that dat includes exactly one such colname.')
  }
  fc <- c('logFC', 'log2FoldChange')
  if (sum(fc %in% colnames(dat)) == 1L) {        # Rename logFC
    colnames(dat)[colnames(dat) %in% fc] <- 'Diff'
  } else {
    stop('dat must include a log fold change column. Recognized colnames for this ',
         'vector include "logFC" and "log2FoldChange". Make sure that dat includes ',
         'exactly one such colname.')
  }
  q <- c('adj.P.Val', 'FDR', 'padj', 'q.value')
  if (sum(q %in% colnames(dat)) == 1L) {         # Rename FDR
    colnames(dat)[colnames(dat) %in% q] <- 'q.value'
  } else {
    stop('dat must include a column for adjusted p-values. Recognized colnames ',
         'for this vector include "q.value", "adj.P.Val", "FDR", "padj", and "FDR". ',
         'Make sure that dat includes exactly one such colname.')
  }
  if (min(dat$q.value) < 0L | max(dat$q.value) > 1L) {
    stop('FDR values must be on [0, 1].')
  }
  if (is.null(xlab)) {
    xlab <- 'Mean Expression'
  }
  df <- dat %>%
    select(Probe, Mean, Diff, q.value) %>%
    na.omit()
  if (nrow(dat) == 0L) {
    stop('dat must have at least one row with non-missing values for AveExpr, logFC, ',
         'and FDR.')
  }

  # Tidy data
  if (!is.null(lfc)) {
    df <- df %>%
      mutate(Direction = ifelse(q.value <= fdr & Diff >= lfc, 'Up',
                                ifelse(q.value <= fdr & -Diff >= lfc, 'Down', 'NA')))
  }

  # Build plot
  size <- probe_ptsize(df)
  alpha <- probe_alpha(df)
  p <- ggplot(df, aes(Mean, Diff, text = Probe)) +
    geom_hline(yintercept = 0L, color = 'grey') +
    labs(title = title, x = xlab, y = expression(log[2]~'Fold Change')) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (!any(df$q.value <= fdr)) {               # Color pts by differential expression?
    warning('No probe meets your fdr threshold. To color data points by differential ',
            'expression, consider raising your fdr cutoff.')
    p <- p + geom_point(size = size, alpha = alpha)
  } else {
    if (is.null(lfc)) {
      p <- p + geom_point(aes(color = q.value <= fdr), size = size, alpha = alpha) +
        scale_color_manual(name = 'FDR',
                           labels = c(paste('>', fdr), paste('\u2264', fdr)),
                           values = c('#444444', pal_d3()(4)[4]),
                           guide = guide_legend(reverse = TRUE, override.aes = list(
                             size = rep(1L, 2L), alpha = rep(1L, 2L))))
    } else {
      suppressWarnings(
        p <- p + geom_hline(yintercept = lfc, linetype = 2L) +
          geom_hline(yintercept = -lfc, linetype = 2L) +
          geom_point(data = df %>% filter(Direction != 'NA'),
                     aes(Mean, Diff, color = Direction, text = Probe),
                     size = size, alpha = alpha) +
          geom_point(data = df %>% filter(Direction == 'NA'),
                     aes(Mean, Diff, text = Probe),
                     color = '#444444', size = size, alpha = alpha) +
          scale_color_manual(guide = FALSE, values = pal_d3()(4)[3:4])
      )
    }
  }
  p <- locate_legend(p, legend)

  # Output
  gg_out(p, hover, legend)

}


#' @rdname plot_md
#' @method plot_md default
#' @S3method plot_md default
#' @importFrom limma getEAWP
#' @importFrom ggsci scale_color_d3

plot_md.default <- function(dat,
                            sample = 1,
                             ctrls = NULL,
                               lfc = NULL,
                             title = NULL,
                              xlab = NULL,
                            legend = 'outside',
                             hover = FALSE) {

  # Preliminaries
  if (is.numeric(sample) & sample > ncol(dat)) {
    stop('Sample number exceeds ncol(dat).')
  }
  if (is.character(sample)) {
    if (!sample %in% colnames(dat)) {
      stop(paste0('Could not detect a sample named "', sample, '" in dat.'))
    } else {
      sample <- which(colnames(dat) == sample)
    }
  }
  if (!is.null(ctrls)) {
    if (length(ctrls) != nrow(dat)) {
      stop('ctrls must be NULL or of length equal to nrow(dat).')
    }
    ctrls <- as.character(ctrls)
    ctrls[ctrls == names(table(ctrls)[which.max(table(ctrls))])] <- '0'
  }
  if (is.null(xlab)) {
    xlab <- 'Mean Expression'
  }

  # Tidy data
  dat <- getEAWP(dat)$expr
  keep <- rowSums(is.finite(dat)) == ncol(dat)
  dat <- dat[keep, , drop = FALSE]
  other <- rowMeans(dat[, -sample])
  df <- data_frame(Probe = rownames(dat),
                    Mean = (other + dat[, sample]) / 2L,
                    Diff = dat[, sample] - other)
  if (!is.null(ctrls)) {
    ctrls <- ctrls[keep]
    df <- df %>% mutate(Control = ctrls)
  }

  # Build plot
  size <- probe_ptsize(df)
  alpha <- probe_alpha(df)
  p <- ggplot(df, aes(Mean, Diff, text = Probe)) +
    geom_hline(yintercept = 0L, color = 'grey') +
    labs(title = title, x = xlab, y = expression(log[2]~'Fold Change')) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (is.null(ctrls)) {
    p <- p + geom_point(size = size, alpha = alpha)
  } else {
    suppressWarnings(
      p <- p + geom_point(data = df %>% filter(Control == '0'),
                          aes(Mean, Diff, text = Probe),
                          size = size, alpha = alpha) +
        geom_point(data = df %>% filter(Control != '0'),
                   aes(Mean, Diff, color = Control, text = Probe),
                   size = 3 * size, alpha = alpha) +
        scale_color_d3()
    )
  }
  if (!is.null(lfc)) {
    p <- p + geom_hline(yintercept = lfc, linetype = 2L) +
      geom_hline(yintercept = -lfc, linetype = 2L)
  }
  p <- locate_legend(p, legend)

  # Output
  gg_out(p, hover, legend)

}


