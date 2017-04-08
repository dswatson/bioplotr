#' Plot Covariance Between Omic and Clinical Data
#'
#' This function creates a heatmap visualizing the strength of associations between
#' the principal components of an omic data matrix and a set of technical and/or
#' biological covariates.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples. For best results, data should be filtered and
#'   normalized in preparation for PCA. For count data, this means undergoing some
#'   sort of variance stabilizing transformation, such as\code{\link[edgeR]{cpm}}
#'   (with \code{log = TRUE}), \code{\link[DESeq2]{vst}, \link[DESeq2]{rlog}}, etc.
#'   Count matrices stored in \code{\link[edgeR]{DGEList}} or \code{\link[DESeq2]{
#'   DESeqDataSet}} objects are automatically extracted and transformed to the
#'   log2-CPM scale, with a warning.
#' @param clin Data frame with rows correponding to samples and columns to
#'   technical and/or biological covariates to test for associations with omic data.
#' @param index String specifying the name of the column in \code{clin} containing
#'   sample names. Only samples in \code{index} with matching column names in
#'   \code{dat} will be considered.
#' @param block String specifying the name of the column in which to find the
#'   blocking variable, should one be accounted for. See Details.
#' @param top Optional number (if > 1) or proportion (if < 1) of most variable probes
#'   to be used for PCA.
#' @param n.pc Number of principal components to include in the figure.
#' @param title Optional plot title.
#' @param hover Show \emph{p}-values by hovering mouse over tiles? If \code{TRUE},
#'   the plot is rendered in HTML and will either open in your browser's graphic
#'   display or appear in the RStudio viewer.
#'
#' @details
#' Strength of association is measured using -log10 \emph{p}-values. When
#' \code{block = NULL}, significance is derived from either Pearson correlation
#' tests (for continuous features) or one-way ANOVA \emph{F}-tests (for categorical
#' features).
#'
#' An optional blocking variable may be provided if samples violate the assumption
#' of independence, e.g. for studies in which subjects are observed at multiple
#' time points. If a blocking variable is identified, it will be factored into all
#' subsequent measures of association. Significance is then evaluated using Pearson
#' partial correlation tests (for continuous features) or repeated measures ANOVA
#' \emph{F}-tests (for categorical features). When supplying a blocking variable,
#' be sure to check that it is not confounded with other features. For instance,
#' biological covariates like sex and age are usually nested within subject, while
#' subject may be nested within other features like batch or treatment group.
#'
#' @examples
#' library(edgeR)
#' data(airway)
#' cnts <- assay(airway)
#' keep <- rowSums(cpm(cnts) > 1) >= 4
#' mat <- cpm(cnts[keep, ], log = TRUE)
#' clin <- colData(airway)[, 1:3]        # Only need the first three cols
#' plot_covar(mat, clin)
#'
#' @export
#' @importFrom edgeR calcNormFactors cpm
#' @importFrom DESeq2 estimateSizeFactors counts
#' @importFrom SummarizedExperiment assay
#' @importFrom limma getEAWP is.fullrank
#' @importFrom purrr map_chr
#' @import dplyr
#' @import ggplot2
#' @importFrom plotly ggplotly
#'

plot_covar <- function(dat,
                       clin,
                       index = 'Sample',
                       block = NULL,
                         top = NULL,
                        n.pc = 10,
                       title = NULL,
                       hover = FALSE) {

  # Preliminaries
  if (ncol(dat) < 3L) {
    stop(paste('dat includes only', ncol(dat), 'samples; need at least 3 for PCA.'))
  }
  if (is(dat, 'DGEList')) {
    keep <- rowSums(dat$counts) > 0L             # Minimal count filter
    dat <- dat[keep, , drop = FALSE]
    if (is.null(dat$samples$norm.factors) |      # Calculate size factors
        all(dat$samples$norm.factors == 1L)) {
      dat <- calcNormFactors(dat)
    }
    dat <- cpm(dat, log = TRUE, prior.count = 1L)
    warning('Transforming raw counts to log2-CPM scale.')
  } else if (is(dat, 'DESeqDataSet')) {
    if (is.null(sizeFactors(dat)) & is.null(normalizationFactors(dat))) {
      dat <- estimateSizeFactors(dat)            # Normalize counts
    }
    dat <- counts(dat, normalized = TRUE)
    keep <- rowMeans(dat) > 0L                   # Minimal count filter
    dat <- dat[keep, , drop = FALSE]
    dat <- cpm(dat, log = TRUE, prior.count = 1L)
    warning('Transforming raw counts to log2-CPM scale.')
  } else if (is(dat, 'DESeqTransform')) {
    dat <- assay(dat)
  } else {
    dat <- getEAWP(dat)$expr
    keep <- rowSums(is.finite(dat)) == ncol(dat)
    dat <- dat[keep, , drop = FALSE]
  }
  clin <- tbl_df(clin)
  if (!index %in% colnames(clin)) {
    stop(paste0('Column "', index, '" not found in clin.'))
  }
  if (any(duplicated(clin[[index]]))) {
    stop('Duplicate sample names detected in index.')
  }
  if (!any(clin[[index]] %in% colnames(dat))) {
    stop('None of the samples in index match colnames in dat.')
  }
  overlap <- length(intersect(clin[[index]], colnames(dat)))
  if (!all(clin[[index]] %in% colnames(dat))) {
    warning(paste('The columns in dat correpond to a subset of the samples in index.',
                  'Analysis will proceed with the', overlap, 'samples for which both',
                  'omic and clinical data were detected.'))
  }
  if (!all(colnames(dat) %in% clin[[index]])) {
    warning(paste('The samples in index correspond to a subset of the columns in dat.',
                  'Analysis will proceed with the', overlap, 'samples for which both',
                  'omic and clinical data were detected.'))
  }
  dat <- dat[, match(clin[[index]], colnames(dat))]
  clin <- clin[, -which(colnames(clin) == index)]
  if (!is.null(block)) {
    if (!block %in% colnames(clin)) {
      stop(paste0('Column "', block, '" not found in clin.'))
    } else {
      i <- which(colnames(clin) == block)
      for (j in seq_along(clin)[-i]) {
        mm <- model.matrix(~ clin[[i]] + clin[[j]])
        if (!is.fullrank(mm)) {
          stop(paste(colnames(clin)[i], 'and', colnames(clin)[j], 'are perfectly',
                     'confounded. Nested covariates generate rank deficient models,',
                     'which cannot be meaningfully evaluated. One or both features',
                     'must be removed or revised.'))
        }
      }
    }
  }
  if (!is.null(top)) {
    if (top > 1L) {
      if (top > nrow(dat)) {
        warning(paste('top exceeds nrow(dat), at least after removing probes with
                      missing values and/or applying a minimal expression filter.
                      Proceeding with the complete', nrow(dat), 'x', ncol(dat), 'matrix.'))
      }
      } else {
        top <- round(top * nrow(dat))
      }
    vars <- rowVars(dat)
    keep <- order(vars, decreasing = TRUE)[seq_len(min(top, nrow(dat)))]
    dat <- dat[keep, , drop = FALSE]
  }
  if (n.pc > max(nrow(dat), ncol(dat))) {
    stop('n.pc cannot exceed max(nrow(dat), ncol(dat))')
  }
  data.frame(Feature = colnames(clin),
               Class = map_chr(seq_along(clin), function(j) {
                 ifelse(is.numeric(clin[[j]]), 'numeric', 'factor')
             })) %>%
    print()
  if (is.null(title)) {
    title <- 'Variation By Feature'
  }

  # Tidy data
  pca <- prcomp(t(dat))                          # PCA, % variance explained
  pve <- map_chr(seq_len(n.pc), function(pc) {
    p <- round(pca$sdev[pc]^2L / sum(pca$sdev^2L) * 100L, 2L)
    paste0('\n(', p, '%)')
  })
  sig <- function(var, pc) {                     # p-val fn
    if (is.null(block)) {
      mod <- lm(pca$x[, pc] ~ clin[[var]])
      ifelse(is.numeric(clin[[var]]),
             -log10(summary(mod)$coef[2, 4]), -log10(anova(mod)[1, 5]))
    } else {
      mod <- lm(pca$x[, pc] ~ clin[[var]] + clin[[block]])
      if (identical(clin[[var]], clin[[block]])) {
        -log10(anova(mod)[1, 5])
      } else {
        ifelse(is.numeric(clin[[var]]),
               -log10(summary(mod)$coef[2, 4]), -log10(anova(mod)[1, 5]))
      }
    }
  }
  df <- expand.grid(Feature = colnames(clin),    # Melt
                    PC = paste0('PC', seq_len(n.pc))) %>%
    rowwise() %>%
    mutate(Association = sig(Feature, PC))       # Populate

  # Build plot
  p <- ggplot(df, aes(PC, Feature, fill = Association, text = Association)) +
    geom_tile() +
    coord_equal() +
    scale_fill_gradientn(colors = c('white', 'pink', 'orange', 'red', 'darkred'),
                           name = expression(~-log[10](italic(p)))) +
    scale_x_discrete(labels = paste0(unique(df$PC), pve)) +
    labs(title = title, x = 'Principal Component') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

  # Output
  if (!hover) {
    print(p)
  } else {
    p <- ggplotly(p, tooltip = 'text', height = 500, width = 900)
    print(p)
  }

}

# Plotly probz:
  # no \n for PCs
  # no -log10(p) for legend
  # automate height/width?

# Maybe change -log10(p) to straightup p?


