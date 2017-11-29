#' Mean-Difference Plot
#'
#' This function plots probewise means vs. log2 fold changes for a test of
#' differential expression or between-sample comparison.
#'
#' @param dat Either a data frame representing the results of a test for
#'   differential expression, or a probe by sample omic data matrix. The former
#'   will render a study-wide MD plot, the latter a between-sample MD plot.
#'   Suitable objects from familiar packages are also acceptable. See Details.
#' @param design Optional design matrix with rows corresponding to samples and
#'   columns to coefficients to be estimated. Only relevant for \code{
#'   \link[edgeR]{DGEList}} objects.
#' @param fdr Optional significance threshold for declaring a probe
#'   differentially expressed. Only relevant for study-wide MD plots.
#' @param lfc Optional effect size threshold for declaring a probe
#'   differentially expressed. Only relevant for study-wide MD plots.
#' @param sample Column number or name specifying which sample in \code{dat} to
#'   compare with the others. Only relevant for between-sample MD plots.
#' @param ctrls Optional vector of length equal to \code{nrow(dat)} indicating
#'   the control status of each probe. Only relevant for between-sample MD
#'   plots.
#' @param probes Optional column number or name specifying where probe names are
#'   stored, presuming they are not stored in \code{rownames(dat)}.
#' @param title Optional plot title.
#' @param xlab Optional label for x-axis.
#' @param legend Legend position. Must be one of \code{"bottom"}, \code{"left"},
#'   \code{"top"}, \code{"right"}, \code{"bottomright"}, \code{"bottomleft"},
#'   \code{"topleft"}, or \code{"topright"}.
#' @param hover Show probe name by hovering mouse over data point? If \code{
#'   TRUE}, the plot is rendered in HTML and will either open in your browser's
#'   graphic display or appear in the RStudio viewer. Probe names are extracted
#'   from \code{dat}.
#'
#' @details
#' MD plots (also known as "Bland-Altman plots" or "MA plots") visualize the
#' relationship between a probe's mean value and its log2 fold change versus
#' some relevant reference group. These figures help to evaluate the symmetry,
#' magnitude, and significance of differential effects across the full omic
#' range.
#'
#' If \code{dat} summarizes the results of a test for differential expression,
#' then each point's x-coordinate correponds to its average expression across
#' all samples, while y-coordinates represent the log2 fold change for the given
#' contrast. Points are colored to distinguish between those that do and do not
#' meet a user-defined FDR threshold. \code{plot_md} accepts output from
#' \code{limma::\link[limma]{topTable}}, \code{edgeR::\link[edgeR]{topTags}}, or
#' \code{DESeq2::\link[DESeq2]{results}}. Alternatively, any object with columns
#' for log fold changes, probewise means, and FDR is acceptable.
#'
#' If \code{dat} is probe by sample matrix or matrix-like object, then \code{
#' sample} must be specified. An artificial array is created by averaging
#' probewise values for all other samples in the data. The figure will then
#' represent the mean vs. the difference of expression values for the specified
#' sample vs. the artificial array. Acceptable inputs for between-sample MD
#' plots include all \code{limma} expression set objects, as well as \code{
#' \link[edgeR]{DGEList}}, \code{\link[DESeq2]{DESeqDataSet}}, and \code{
#' \link[DESeq2]{DESeqTransform}} objects.
#'
#' @references
#' Bolstad, B.M., Irizarry, R.A., Åstrand, M. & Speed, T.P. (2003).
#' \href{https://www.ncbi.nlm.nih.gov/pubmed/12538238}{A comparison of
#' normalization methods for high density oligonucleotide array data based on
#' variance and bias}. \emph{Bioinformatics}, \emph{19}(2): 185–193.
#'
#' Dudoit, S., Yang, Y.H., Callow, M.J. & Speed, T.P. (2002).
#' \href{https://www.jstor.org/stable/24307038?seq=1#page_scan_tab_contents}{
#' Statistical methods for identifying differentially expressed genes in
#' replicated cDNA microarray experiments}. \emph{Stat. Sin.}, \strong{12},
#' 111–140.
#'
#' Martin, B.J. & Altman, D.G. (1986).
#' \href{http://www.sciencedirect.com/science/article/pii/S0140673686908378}{
#' Statistical methods for assessing agreement between two methods of clinical
#' measurement}. \emph{The Lancet}, \emph{327}(8476): 307–310.
#'
#' Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., & Smyth, G.K.
#' (2015).
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
#' \code{\link[limma]{plotMD}}, \code{\link[DESeq2]{plotMA}},
#' \code{\link[Glimma]{glMDPlot}}
#'
#' @export
#' @importFrom ggsci scale_color_d3 pal_d3
#' @import dplyr
#' @import ggplot2
#'

plot_md <- function(dat,
                    title = NULL,
                   legend = 'right', ...) {

  # Preliminaries
  if (title %>% is.null) {
    title <- 'Mean-Difference Plot'
  }
  loc <- c('bottom', 'left', 'top', 'right',
           'bottomright', 'bottomleft', 'topleft', 'topright')
  if (!legend %in% loc) {
    stop('legend must be one of ', stringify(loc, 'or'), '.')
  }

  # Method
  UseMethod('plot_md')

}


#' @rdname plot_md
#' @export
#' @importFrom edgeR calcNormFactors aveLogCPM estimateDisp cpm

plot_md.DGEList <- function(dat,
                            design = NULL,
                            sample = 1L,
                             ctrls = NULL,
                               lfc = NULL,
                             title = NULL,
                              xlab = NULL,
                            legend = 'right',
                             hover = FALSE) {

  # Preliminaries
  if (sample %>% is.numeric && sample > ncol(dat)) {
    stop('Sample number exceeds ncol(dat).')
  }
  if (sample %>% is.character) {
    if (!sample %in% colnames(dat)) {
      stop('Could not detect a sample named "', sample, '" in dat.')
    } else {
      sample <- which(colnames(dat) == sample)
    }
  }
  if (!(ctrls %>% is.null)) {
    if (length(ctrls) != nrow(dat)) {
      stop('ctrls must be NULL or of length equal to nrow(dat).')
    }
    ctrls <- as.character(ctrls)
    ctrls[ctrls == names(which.max(table(ctrls)))] <- '0'
  }
  if (xlab %>% is.null) {
    xlab <- expression('Mean'~log[2]*'-CPM')
  }

  # Tidy data
  keep <- rowSums(dat$counts) > 1L               # Minimal count filter
  dat <- dat[keep, ]
  other <- dat[, -sample]
  other <- calcNormFactors(other)
  if (other$tagwise.dispersion %>% is.null) {    # Estimate dispersions for aveLogCPM
    if (design %>% is.null && !(dat$group %>% is.null)) {
      design <- model.matrix(~ dat$group)
    }
    if (!(design %>% is.null)) {
      design <- design[-sample, ]
    }
    if (design %>% is.null) {
      if (other$common.dispersion %>% is.null) {
        other <- estimateCommonDisp(other)
      }
      other <- estimateTagwiseDisp(other)
    } else {
      other <- estimateDisp(other, design = design)
    }
  }
  other <- aveLogCPM(other, prior.count = 1L,
                     dispersion = other$tagwise.dispersion)
  dat <- cpm(dat, log = TRUE, prior.count = 1L)
  df <- data_frame(Probe = rownames(dat),
                    Mean = (other + dat[, sample]) / 2L,
                    Diff = dat[, sample] - other)
  if (!(ctrls %>% is.null)) {
    ctrls <- ctrls[keep]
    df <- df %>% mutate(Control = ctrls)
  }

  # Build plot
  size <- pt_size(df)
  alpha <- pt_alpha(df)
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
    p <- p + geom_hline(yintercept = lfc, linetype = 'dashed') +
      geom_hline(yintercept = -lfc, linetype = 'dashed')
  }

  # Output
  gg_out(p, hover, legend)

}


#' @rdname plot_md
#' @export
#' @importFrom edgeR aveLogCPM cpm

plot_md.DESeqDataSet <- function(dat,
                                 sample = 1L,
                                  ctrls = NULL,
                                    lfc = NULL,
                                  title = NULL,
                                   xlab = NULL,
                                 legend = 'right',
                                  hover = FALSE) {

  # Preliminaries
  if (sample %>% is.null && sample > ncol(dat)) {
    stop('Sample number exceeds ncol(dat).')
  }
  if (sample %>% is.character) {
    if (!sample %in% colnames(dat)) {
      stop('Could not detect a sample named "', sample, '" in dat.')
    } else {
      sample <- which(colnames(dat) == sample)
    }
  }
  if (!(ctrls %>% is.null)) {
    if (length(ctrls) != nrow(dat)) {
      stop('ctrls must be NULL or of length equal to nrow(dat).')
    }
    ctrls <- as.character(ctrls)
    ctrls[ctrls == names(which.max(table(ctrls)))] <- '0'
  }
  if (xlab %>% is.null) {
    xlab <- expression('Mean'~log[2]*'-CPM')
  }
  suppressPackageStartupMessages(require(DESeq2))
  if (sizeFactors(dat) %>% is.null && normalizationFactors(dat) %>% is.null) {
    dat <- estimateSizeFactors(dat)
  }
  if (dispersions(dat) %>% is.null) {
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
  if (!(ctrls %>% is.null)) {
    ctrls <- ctrls[keep]
    df <- df %>% mutate(Control = ctrls)
  }

  # Build plot
  size <- pt_size(df)
  alpha <- pt_alpha(df)
  p <- ggplot(df, aes(Mean, Diff, text = Probe)) +
    geom_hline(yintercept = 0L, color = 'grey') +
    labs(title = title, x = xlab, y = expression(log[2]~'Fold Change')) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (ctrls %>% is.null) {
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
  if (!(lfc %>% is.null)) {
    p <- p + geom_hline(yintercept = lfc, linetype = 'dashed') +
      geom_hline(yintercept = -lfc, linetype = 'dashed')
  }

  # Output
  gg_out(p, hover, legend)

}


#' @rdname plot_md
#' @export

plot_md.DESeqTransform <- function(dat,
                                   sample = 1L,
                                    ctrls = NULL,
                                      lfc = NULL,
                                    title = NULL,
                                     xlab = NULL,
                                   legend = 'right',
                                    hover = FALSE) {

  # Preliminaries
  if (xlab %>% is.null) {
    xlab <- expression('Mean Transformed Counts')
  }

  # Tidy data
  suppressPackageStartupMessages(require(SummarizedExperiment))
  dat <- assay(dat)

  # Export
  plot_md.default(dat = dat, sample = sample, ctrls = ctrls, lfc = lfc,
                  title = title, xlab = xlab, legend = legend, hover = hover)

}


#' @rdname plot_md
#' @export

plot_md.DESeqResults <- function(dat,
                                 fdr = 0.05,
                                 lfc = NULL,
                               title = NULL,
                                xlab = NULL,
                              legend = 'right',
                               hover = FALSE) {

  # Preliminaries
  if (xlab %>% is.null) {
    xlab <- 'Mean of Normalized Counts'
  }
  dat <- dat %>% na.omit(.)
  if (nrow(dat) == 0L) {
    stop('dat must have at least one row with non-missing values for ',
         'baseMean, log2FoldChange, and padj.')
  }

  # Tidy data
  dat$Probe <- rownames(dat)
  df <- dat %>%
    as_tibble(.) %>%
    rename(Mean = baseMean,
           Diff = log2FoldChange,
        q.value = padj) %>%
    select(Probe, Mean, Diff, q.value)
  if (!(lfc %>% is.null)) {
    df <- df %>%
      mutate(Direction = ifelse(q.value <= fdr && Diff >= lfc, 'Up',
                                ifelse(q.value <= fdr && -Diff >= lfc, 'Down', 'None')))
  }

  # Build plot
  size <- pt_size(df)
  alpha <- pt_alpha(df)
  p <- ggplot(df, aes(Mean, Diff, text = Probe)) +
    scale_x_log10() +
    geom_hline(yintercept = 0L, color = 'grey') +
    labs(title = title, x = xlab, y = expression(log[2]~'Fold Change')) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (!any(df$q.value <= fdr)) {               # Color pts by differential expression?
    warning('No probe meets your fdr threshold. To color data points by ',
            'differential expression, consider raising your fdr cutoff.')
    p <- p + geom_point(size = size, alpha = alpha)
  } else {
    if (lfc %>% is.null) {
      p <- p + geom_point(aes(color = q.value <= fdr), size = size, alpha = alpha) +
        scale_color_manual(name = expression(italic(q)*'-value'),
                         labels = c(paste('>', fdr), paste('\u2264', fdr)),
                         values = c('#444444', pal_d3()(4)[4]),
                          guide = guide_legend(reverse = TRUE, override.aes = list(
                           size = rep(1L, 2L), alpha = rep(1L, 2L))))
    } else {
      suppressWarnings(
        p <- p + geom_hline(yintercept = lfc, linetype = 'dashed') +
          geom_hline(yintercept = -lfc, linetype = 'dashed') +
          geom_point(data = filter(df, Direction != 'None'),
                     aes(Mean, Diff, color = Direction, text = Probe),
                     size = size, alpha = alpha) +
          geom_point(data = filter(df, Direction == 'None'),
                     aes(Mean, Diff, text = Probe),
                     color = '#444444', size = size, alpha = alpha) +
          scale_color_manual(guide = FALSE, values = pal_d3()(4)[3:4])
      )
    }
  }

  # Output
  gg_out(p, hover, legend)

}


#' @rdname plot_md
#' @export

plot_md.TopTags <- function(dat,
                            fdr = 0.05,
                            lfc = NULL,
                          title = NULL,
                           xlab = NULL,
                         legend = 'right',
                          hover = FALSE) {

  # Preliminaries
  if (xlab %>% is.null) {
    xlab <- expression('Mean'~log[2]*'-CPM')
  }

  # Export
  dat <- dat %>% as.data.frame(.)
  plot_md.data.frame(dat = dat, fdr = fdr, lfc = lfc,
                     title = title, xlab = xlab, legend = legend, hover = hover)

}


#' @rdname plot_md
#' @export

plot_md.data.frame <- function(dat,
                               probes = NULL,
                                  fdr = 0.05,
                                  lfc = NULL,
                                title = NULL,
                                 xlab = NULL,
                               legend = 'right',
                                hover = FALSE) {

  # Preliminaries
  if (probes %>% is.null) {
    if (rownames(dat) %>% is.null ||
        rownames(dat) %>% identical(as.character(seq_len(nrow(dat))))) {
      stop('If dat does not have rownames, then the column of unique probe ',
           'identifiers must be specified using the probes argument.')
    } else {
      dat$Probe <- rownames(dat)
    }
  } else {
    if (probes %>% is.numeric) {
      if (probes > ncol(dat)) {
        stop('Column number for probes exceeds ncol(dat).')
      } else {
        colnames(dat)[probes] <- 'Probe'
      }
    } else {
      if (!probes %in% colnames(dat)) {
        stop('Could not detect a column named "', probes, '" in dat.')
      } else {
        colnames(dat)[colnames(dat) == probes] <- 'Probe'
      }
    }
  }
  if ('baseMean' %in% colnames(dat)) {
    dat$baseMean <- log2(dat$baseMean / 1e6L)
  }
  avg <- c('AveExpr', 'baseMean', 'logCPM', 'AvgExpr', 'AvgMeth')
  if (sum(avg %in% colnames(dat)) == 1L) {       # Rename AvgExpr
    colnames(dat)[colnames(dat) %in% avg] <- 'Mean'
  } else {
    stop('dat must include a column for average expression by probe. ',
         'Recognized colnames for this vector include ', stringify(avg, 'and'),
         '. Make sure that dat includes exactly one such colname.')
  }
  fc <- c('logFC', 'log2FoldChange')
  if (sum(fc %in% colnames(dat)) == 1L) {        # Rename logFC
    colnames(dat)[colnames(dat) %in% fc] <- 'Diff'
  } else {
    stop('dat must include a log fold change column. Recognized colnames for ',
         'this vector include ', stringify(fc, 'and'), '. Make sure that dat ',
         'includes exactly one such colname.')
  }
  q <- c('adj.P.Val', 'FDR', 'padj', 'q.value')
  if (sum(q %in% colnames(dat)) == 1L) {         # Rename FDR
    colnames(dat)[colnames(dat) %in% q] <- 'q.value'
  } else {
    stop('dat must include a column for adjusted p-values. Recognized ',
         'colnames for this vector include ', stringify(q, 'and'), '. Make ',
         'sure that dat includes exactly one such colname.')
  }
  if (min(dat$q.value) < 0L || max(dat$q.value) > 1L) {
    stop('Adjusted p-values must be on [0, 1].')
  }
  if (xlab %>% is.null) {
    xlab <- 'Mean Expression'
  }
  df <- dat %>%
    select(Probe, Mean, Diff, q.value) %>%
    na.omit(.)
  if (nrow(dat) == 0L) {
    stop('dat must have at least one row with non-missing values for AveExpr, ',
         'logFC, and FDR.')
  }

  # Tidy data
  if (!(lfc) %>% is.null) {
    df <- df %>%
      mutate(Direction = ifelse(q.value <= fdr & Diff >= lfc, 'Up',
                                ifelse(q.value <= fdr & -Diff >= lfc, 'Down', 'None')))
  }

  # Build plot
  size <- pt_size(df)
  alpha <- pt_alpha(df)
  p <- ggplot(df, aes(Mean, Diff, text = Probe)) +
    geom_hline(yintercept = 0L, color = 'grey') +
    labs(title = title, x = xlab, y = expression(log[2]~'Fold Change')) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (!any(df$q.value <= fdr)) {               # Color pts by differential expression?
    warning('No probe meets your fdr threshold. To color data points by ',
            'differential expression, consider raising your fdr cutoff.')
    p <- p + geom_point(size = size, alpha = alpha)
  } else {
    if (lfc %>% is.null) {
      p <- p + geom_point(aes(color = q.value <= fdr), size = size, alpha = alpha) +
        scale_color_manual(name = expression(italic(q)*'-value'),
                         labels = c(paste('>', fdr), paste('\u2264', fdr)),
                         values = c('#444444', pal_d3()(4L)[4L]),
                          guide = guide_legend(reverse = TRUE, override.aes = list(
                           size = rep(1L, 2L), alpha = rep(1L, 2L))))
    } else {
      suppressWarnings(
        p <- p + geom_hline(yintercept = lfc, linetype = 'dashed') +
          geom_hline(yintercept = -lfc, linetype = 'dashed') +
          geom_point(data = filter(df, Direction != 'None'),
                     aes(Mean, Diff, color = Direction, text = Probe),
                     size = size, alpha = alpha) +
          geom_point(data = filter(df, Direction == 'None'),
                     aes(Mean, Diff, text = Probe),
                     color = '#444444', size = size, alpha = alpha) +
          scale_color_manual(guide = FALSE, values = pal_d3()(4)[3:4])
      )
    }
  }

  # Output
  gg_out(p, hover, legend)

}


#' @rdname plot_md
#' @export
#' @importFrom limma getEAWP

plot_md.default <- function(dat,
                            sample = 1L,
                             ctrls = NULL,
                               lfc = NULL,
                             title = NULL,
                              xlab = NULL,
                            legend = 'right',
                             hover = FALSE) {

  # Preliminaries
  if (sample %>% is.numeric && sample > ncol(dat)) {
    stop('Sample number exceeds ncol(dat).')
  }
  if (sample %>% is.charater) {
    if (!sample %in% colnames(dat)) {
      stop('Could not detect a sample named "', sample, '" in dat.')
    } else {
      sample <- which(colnames(dat) == sample)
    }
  }
  if (!(ctrls %>% is.null)) {
    if (length(ctrls) != nrow(dat)) {
      stop('ctrls must be NULL or of length equal to nrow(dat).')
    }
    ctrls <- ctrls %>% as.character(.)
    ctrls[ctrls == names(which.max(table(ctrls)))] <- '0'
  }
  if (xlab %>% is.null) {
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
  if (!(ctrls %>% is.null)) {
    ctrls <- ctrls[keep]
    df <- df %>% mutate(Control = ctrls)
  }

  # Build plot
  size <- pt_size(df)
  alpha <- pt_alpha(df)
  p <- ggplot(df, aes(Mean, Diff, text = Probe)) +
    geom_hline(yintercept = 0L, color = 'grey') +
    labs(title = title, x = xlab, y = expression(log[2]~'Fold Change')) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (ctrls %>% is.null) {
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
  if (!(lfc %>% is.null)) {
    p <- p + geom_hline(yintercept = lfc, linetype = 'dashed') +
      geom_hline(yintercept = -lfc, linetype = 'dashed')
  }

  # Output
  gg_out(p, hover, legend)

}


