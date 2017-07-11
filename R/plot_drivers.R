#' Plot Drivers of Omic Variation
#'
#' This function visualizes the strength of associations between the principal
#' components of an omic data matrix and a set of biological and/or technical
#' features.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples. It is strongly recommended that data be
#'   filtered and normalized prior to plotting. Raw counts stored in \code{
#'   \link[edgeR]{DGEList}} or \code{\link[DESeq2]{DESeqDataSet}} objects are
#'   automatically extracted and transformed to the log2-CPM scale, with a
#'   warning.
#' @param clin Data frame or matrix with rows correponding to samples and
#'   columns to technical and/or biological features to test for associations
#'   with omic data.
#' @param index String specifying the name of the column in \code{clin}
#'   containing sample names. Only samples in \code{index} with matching column
#'   names in \code{dat} will be considered.
#' @param block String specifying the name of the column in which to find the
#'   blocking variable, should one be accounted for. See Details.
#' @param top Optional number (if > 1) or proportion (if < 1) of most variable
#'   probes to be used for PCA.
#' @param n.pc Number of principal components to include in the figure.
#' @param label Print association statistics over tiles?
#' @param alpha Optional significance threshold to impose on associations. Those
#'   with \emph{p}-values (optionally adjusted) less than or equal to \code{
#'   alpha} are outlined in black.
#' @param p.adj Optional \emph{p}-value adjustment for multiple testing. Options
#'   include \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{
#'   "bonferroni"}, \code{"BH"}, \code{"BY"}, and \code{"fdr"}. See \code{
#'   \link[stats]{p.adjust}}.
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"bottom"}, \code{"left"},
#'   \code{"top"}, \code{"right"}, \code{"bottomright"}, \code{"bottomleft"},
#'   \code{"topleft"}, or \code{"topright"}.
#' @param hover Show \emph{p}-values by hovering mouse over tiles? If \code{
#'   TRUE}, the plot is rendered in HTML and will either open in your browser's
#'   graphic display or appear in the RStudio viewer.
#'
#' @details
#' Strength of association is measured by -log \emph{p}-values, optionally
#' adjusted for multiple testing. When \code{block = NULL}, significance is
#' derived from either Pearson correlation tests (for continuous features) or
#' one-way ANOVA \emph{F}-tests (for categorical features).
#'
#' An optional blocking variable may be provided if samples violate the
#' assumption of independence, e.g. for studies in which subjects are observed
#' at multiple time points. If a blocking variable is identified, it will be
#' factored into all subsequent measures of association. Significance is then
#' evaluated using Pearson partial correlation tests (for continuous features)
#' or repeated measures ANOVA \emph{F}-tests (for categorical features). When
#' supplying a blocking variable, be sure to check that it is not confounded
#' with other features. For instance, features like sex and age are usually
#' nested within subject, while subject may be nested within other variables
#' like batch or treatment group.
#'
#' @examples
#' library(edgeR)
#' data(airway)
#' cnts <- assay(airway)
#' keep <- rowSums(cpm(cnts) > 1) >= 4
#' mat <- cpm(cnts[keep, ], log = TRUE)
#' clin <- colData(airway)[, 1:3]        # Only need the first three cols
#' plot_drivers(mat, clin)
#'
#' @export
#' @importFrom limma is.fullrank
#' @importFrom purrr map_chr
#' @import dplyr
#' @import ggplot2
#'

plot_drivers <- function(dat,
                         clin,
                         index = 'Sample',
                         block = NULL,
                           top = NULL,
                          n.pc = 10,
                         label = FALSE,
                         alpha = NULL,
                         p.adj = NULL,
                         title = NULL,
                        legend = 'right',
                         hover = FALSE) {

  # Preliminaries
  if (ncol(dat) < 3L) {
    stop('dat includes only ', ncol(dat), ' samples; need at least 3 for PCA.')
  }
  dat <- matrixize(dat)
  clin <- as_tibble(clin)
  if (!index %in% colnames(clin)) {
    stop(paste0('Column "', index, '" not found in clin.'))
  }
  if (any(clin[[index]] %>% duplicated)) {
    stop('Duplicate sample names detected in index.')
  }
  if (!any(clin[[index]] %in% colnames(dat))) {
    stop('None of the samples in index match colnames in dat.')
  }
  overlap <- length(intersect(clin[[index]], colnames(dat)))
  if (!all(clin[[index]] %in% colnames(dat))) {
    warning('The columns in dat correpond to a subset of the samples in ',
            'index. Analysis will proceed with the ', overlap, ' samples  for ',
            'which both omic and clinical data were detected.')
  }
  if (!all(colnames(dat) %in% clin[[index]])) {
    warning('The samples in index correspond to a subset of the columns in ',
            'dat. Analysis will proceed with the ', overlap, ' samples for ',
            'which both omic and clinical data were detected.')
  }
  dat <- dat[, match(clin[[index]], colnames(dat))]
  clin <- clin[, -which(colnames(clin) == index), drop = FALSE]
  if (!is.null(block)) {
    if (!block %in% colnames(clin)) {
      stop(paste0('Column "', block, '" not found in clin.'))
    } else {
      i <- which(colnames(clin) == block)
      for (j in seq_along(clin)[-i]) {
        mm <- model.matrix(~ clin[[i]] + clin[[j]])
        if (!is.fullrank(mm)) {
          stop(colnames(clin)[i],  'and ', colnames(clin)[j], ' are perfectly ',
               'confounded. Nested covariates generate rank deficient models, ',
               'which cannot be meaningfully evaluated. One or both features ',
               'must be removed or revised.')
        }
      }
    }
  }
  if (!(top %>% is.null)) {                      # Filter by variance?
    dat <- var_filt(dat, top, robust = FALSE)
  }
  if (n.pc > max(nrow(dat), ncol(dat))) {
    stop('n.pc cannot exceed max(nrow(dat), ncol(dat))')
  }
  if (!(alpha %>% is.null)) {
    if (alpha <= 0 || alpha >= 1) {
      stop('alpha must be numeric on (0, 1).')
    }
  }
  if (!(p.adj %>% is.null)) {
    p_adj <- c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr')
    if (!p.adj %in% p_adj) {
      stop('p.adj must be one of ', stringify(p_adj, 'or'), '. See ?p.adjust.')
    }
  }
  if (title %>% is.null) {
    title <- 'Variation By Feature'
  }
  loc <- c('bottom', 'left', 'top', 'right',
           'bottomright', 'bottomleft', 'topleft', 'topright')
  if (!legend %in% loc) {
    stop('legend must be one of ', stringify(loc, 'or'), '.')
  }
  data_frame(Feature = colnames(clin),           # Be apprised
               Class = clin %>% map_chr(class)) %>%
    print(n = nrow(.))

  # Tidy data
  pca <- prcomp(t(dat))                          # PCA, % variance explained
  pve <- seq_len(n.pc) %>% map_chr(function(pc) {
    p <- round(pca$sdev[pc]^2L / sum(pca$sdev^2L) * 100L, 2L)
    paste0('\n(', p, '%)')
  })
  sig <- function(var, pc) {                     # p-val fn
    if (block %>% is.null) {
      mod <- lm(pca$x[, pc] ~ clin[[var]])
      ifelse(clin[[var]] %>% is.numeric,
             summary(mod)$coef[2L, 4L], anova(mod)[1L, 5L])
    } else {
      mod <- lm(pca$x[, pc] ~ clin[[var]] + clin[[block]])
      if (clin[[var]] %>% identical(clin[[block]])) {
        anova(mod)[1L, 5L]
      } else {
        ifelse(clin[[var]] %>% is.numeric,
               summary(mod)$coef[2L, 4L], anova(mod)[1L, 5L])
      }
    }
  }
  df <- expand.grid(Feature = colnames(clin),    # Melt
                         PC = paste0('PC', seq_len(n.pc))) %>%
    rowwise() %>%
    mutate(Association = sig(Feature, PC),       # Populate
           Significant = FALSE)
  if (!(alpha %>% is.null)) {
    if (!(p.adj %>% is.null)) {
      df <- df %>% mutate(Association = p.adjust(Association, method = p.adj))
    }
    df <- df %>% mutate(Significant = ifelse(Association <= alpha, TRUE, FALSE))
  }
  df <- df %>% mutate(Association = -log(Association))

  # Build plot
  p <- ggplot(df, aes(PC, Feature, fill = Association, text = Association,
                      color = Significant)) +
    geom_tile(size = 1L, width = 0.9, height = 0.9) +
    coord_equal() +
    scale_fill_gradientn(colors = c('white', 'pink', 'orange', 'red', 'darkred'),
                           name = expression(~-log(italic(p)))) +
    scale_color_manual(values = c('grey90', 'black')) +
    scale_x_discrete(labels = paste0(unique(df$PC), pve)) +
    guides(color = FALSE) +
    labs(title = title, x = 'Principal Component') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (label) {
    p <- p + geom_text(aes(label = round(Association, 2L)))
  }

  # Output
  gg_out(p, hover, legend)

}


# Fit multivariate model?
# Fages & Ferrari, 2014: https://link.springer.com/article/10.1007/s11306-014-0647-9
