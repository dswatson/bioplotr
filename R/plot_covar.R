#' Plot associations between omic and clinical data
#'
#' This function creates a heatmap visualizing the strength of associations
#' between the principal components of an omic data matrix and a set of
#' technical and/or biological covariates.
#'
#' @param dat Omic data matrix with rows corresponding to probes and columns
#'   to samples. For best results, data should be normalized and filtered in
#'   preparation for PCA.
#' @param clin Data frame with rows correponding to samples and columns to
#'   technical and/or biological covariates to test for associations with omic data.
#' @param index String specifying the name of the column in \code{clin} containing
#'   sample names. Only samples in \code{index} with matching column names in
#'   \code{dat} will be considered.
#' @param block String specifying the name of the column in which to find the
#'   blocking variable, should one be accounted for. See Details.
#' @param n.pc Number of principal components to include in the figure.
#' @param main Optional plot title.
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
#' data(Nevins)
#' mat <- exprs(Nevins)
#' clin <- pData(Nevins)
#' clin$Sample <- rownames(clin)
#' plot_covar(mat, clin)
#'
#' @export
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
                        n.pc = 10,
                        main = NULL,
                       hover = FALSE) {

  # Preliminaries
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
        mm <- model.matrix(~ clin[, i] + clin[, j])
        if (!is.fullrank(mm)) {
          stop(paste(colnames(clin)[i], 'and', colnames(clin)[j], 'are perfectly',
                     'confounded. Nested covariates generate rank deficient models,',
                     'which cannot be meaningfully evaluated. One or both features',
                     'must be removed or revised.'))
        }
      }
    }
  }
  if (n.pc > max(nrow(dat), ncol(dat))) {
    stop('n.pc cannot exceed max(nrow(dat), ncol(dat))')
  }
  data.frame(Feature = colnames(clin),
               Class = map_chr(seq_along(clin), function(j) {
                 ifelse(is.numeric(clin[, j]), 'numeric', 'factor')
             })) %>%
    print()
  if (is.null(main)) main <- 'Variation by Feature'

  # Tidy data
  dat <- getEAWP(dat)$expr
  keep <- rowSums(is.finite(dat)) == ncol(dat)
  dat <- dat[keep, , drop = FALSE]
  pca <- prcomp(t(dat))                        # PCA
  pve <- map_chr(seq_len(n.pc), function(pc) {
    p <- round(pca$sdev[pc]^2L / sum(pca$sdev^2L) * 100L, 2L)
    paste0('\n(', p, '%)')
  })
  sig <- function(var, pc) {                   # p-val fn
    if (is.null(block)) {
      mod <- lm(pca$x[, pc] ~ clin[, var])
      ifelse(is.numeric(clin[, var]),
             -log10(summary(mod)$coef[2, 4]), -log10(anova(mod)[1, 5]))
    } else {
      mod <- lm(pca$x[, pc] ~ clin[, var] + clin[[block]])
      if (identical(clin[, var], clin[[block]])) {
        -log10(anova(mod)[1, 5])
      } else {
        ifelse(is.numeric(clin[, var]),
               -log10(summary(mod)$coef[2, 4]), -log10(anova(mod)[1, 5]))
      }
    }
  }
  df <- expand.grid(Feature = colnames(clin),  # Melt
                    PC = paste0('PC', seq_len(n.pc))) %>%
    rowwise() %>%
    mutate(Association = sig(Feature, PC))     # Populate

  # Build plot
  suppressWarnings(
    p <- ggplot(df, aes(PC, Feature, fill = Association, text = Association)) +
      geom_tile() +
      coord_equal() +
      scale_fill_gradientn(colors = c('white', 'pink', 'orange', 'red', 'darkred'),
                           name = expression(~-log[10](italic(p)))) +
      scale_x_discrete(labels = paste0(unique(df$PC), pve)) +
      labs(title = main, x = 'Principal Component') +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
  )

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

# When shiny-ifying: change n.pc?


