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
#' @param unblock Column name(s) of one or more features for which the \code{
#'   block} covariate should not be applied, if one was supplied. See Details.
#' @param kernel The kernel generating function, if using KPCA. Options include
#'   \code{"rbfdot", "polydot", "tanhdot", "vanilladot", "laplacedot",
#'   "besseldot", "anovadot",} and \code{"splinedot"}. To run normal PCA,
#'   set to \code{NULL}. See Details.
#' @param kpar A named list of arguments setting parameters for the kernel
#'   function. Only relevant if \code{kernel} is not \code{NULL}. See Details.
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
#' factored into all subsequent measures of association, except those explicitly
#' exempted by the \code{unblock} argument. Significance is evaluated using
#' Pearson partial correlation tests (for continuous features) or repeated
#' measures ANOVA \emph{F}-tests (for categorical features).
#'
#' When supplying a blocking variable, be careful to consider potential
#' confounding effects. For instance, features like sex and age are usually
#' nested within subject, while subject may be nested within other variables
#' like batch or treatment group. The \code{block} and \code{unblock} arguments
#' are designed to help parse out these relationships.
#'
#' If \code{kernel} is non-\code{NULL}, then KPCA is used instead of PCA. See
#' \code{\link{plot_kpca}} for more info. Details on kernel functions and their
#' input parameters can be found in \code{kernlab::\link[kernlab]{dots}}.
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
#' @seealso
#' \code{\link{plot_pca}}, \code{\link{plot_kpca}}
#'
#' @export
#' @importFrom limma is.fullrank
#' @importFrom purrr map_chr
#' @importFrom kernlab rbfdot
#' @importFrom kernlab polydot
#' @importFrom kernlab tanhdot
#' @importFrom kernlab vanilladot
#' @importFrom kernlab laplacedot
#' @importFrom kernlab besseldot
#' @importFrom kernlab anovadot
#' @importFrom kernlab splinedot
#' @importFrom kernlab kernelMatrix
#' @importFrom kernlab kpca
#' @importFrom kernlab eig
#' @importFrom kernlab rotated
#' @import dplyr
#' @import ggplot2
#'

plot_drivers <- function(dat,
                         clin,
                         index = 'Sample',
                         block = NULL,
                       unblock = NULL,
                        kernel = NULL,
                          kpar = NULL,
                           top = NULL,
                          n.pc = 5L,
                         label = FALSE,
                         alpha = NULL,
                         p.adj = NULL,
                         title = 'Variation By Feature',
                        legend = 'right',
                         hover = FALSE) {

  # Preliminaries
  if (ncol(dat) < 3L) {
    stop('dat includes only ', ncol(dat), ' samples; need at least 3 for PCA.')
  }
  dat <- matrixize(dat)
  clin <- as_tibble(clin)
  if (!index %in% colnames(clin)) {
    stop('Column "', index, '" not found in clin.')
  }
  if (any(clin[[index]] %>% duplicated)) {
    stop('Duplicate sample names detected in index.')
  }
  if (!any(clin[[index]] %in% colnames(dat))) {
    stop('None of the samples in index match colnames in dat.')
  }
  overlap <- intersect(clin[[index]], colnames(dat))
  n <- length(overlap)
  if (!all(clin[[index]] %in% colnames(dat))) {
    warning('The columns in dat correpond to a subset of the samples in ',
            'index. Analysis will proceed with the ', n, ' samples  for ',
            'which both omic and clinical data were detected.')
  }
  if (!all(colnames(dat) %in% clin[[index]])) {
    warning('The samples in index correspond to a subset of the columns in ',
            'dat. Analysis will proceed with the ', n, ' samples for ',
            'which both omic and clinical data were detected.')
  }
  dat <- dat[, overlap]
  clin <- clin[clin[[index]] %in% overlap, ] %>%
    select(-index)
  if (!(block %>% is.null)) {
    if (!block %in% colnames(clin)) {
      stop(paste0('Column "', block, '" not found in clin.'))
    } else {
      if (unblock %>% is.null) {
        j_idx <- which(colnames(clin) != block)
        for (j in j_idx) {
          mm <- model.matrix(~ clin[[block]] + clin[[j]])
          if (!is.fullrank(mm)) {
            stop('"', block,  '"', ' and "', colnames(clin)[j], '" are ', 
                 'perfectly confounded. Consider using the unblock argument.')
          }
        }
      } else {
        if (!all(unblock %in% colnames(clin))) {
          stop('The following column(s) not found in clin: ',
               stringify(unblock[!unblock %in% colnames(clin)], 'and'))
        }
      }
    }
  }
  kernels <- c('rbfdot', 'polydot', 'tanhdot', 'vanilladot', 'laplacedot', 
               'besseldot', 'anovadot', 'splinedot')
  if (!kernel %in% kernels) {
    stop('kernel must be one of ', stringify(kernels, 'or'), '. ', 
         'For more info, see ?plot_kpca or ?kernlab::dots.')
  }
  # SOME WARNING ABOUT KPAR?
  if (!(top %>% is.null)) {                        # Filter by variance?
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
  loc <- c('bottom', 'left', 'top', 'right',
           'bottomright', 'bottomleft', 'topleft', 'topright')
  if (!legend %in% loc) {
    stop('legend must be one of ', stringify(loc, 'or'), '.')
  }
  tibble(Feature = colnames(clin),                 # Be apprised
           Class = clin %>% map_chr(class)) %>%
    print(n = nrow(.))

  # Tidy data
  if (kernel %>% is.null) {
    pca <- prcomp(t(dat))                          # PCA, % variance explained
    pve <- seq_len(n.pc) %>% map_chr(function(pc) {
      p <- round(pca$sdev[pc]^2L / sum(pca$sdev^2L) * 100L, 2L)
      paste0('\n(', p, '%)')
    })
  } else {
    if (kernel == 'rbfdot') {                      # Initialize kernel function
      if (kpar %>% is.null) {
        kpar <- list(sigma = 1L)
      }
      kf <- rbfdot(unlist(kpar))
    } else if (kernel == 'polydot') {
      if (kpar %>% is.null) {
        kpar <- list(degree = 1L, scale = 1L, offset = 1L)
      }
      kf <- polydot(unlist(kpar))
    } else if (kernel == 'tanhdot') {
      if (kpar %>% is.null) {
        kpar <- list(scale = 1L, offset = 1L)
      }
      kf <- tanhdot(unlist(kpar))
    } else if (kernel == 'vanilladot') {
      kf <- vanilladot()
    } else if (kernel == 'laplacedot') {
      if (kpar %>% is.null) {
        kpar <- list(sigma = 1L)
      }
      kf <- laplacedot(unlist(kpar))
    } else if (kernel == 'besseldot') {
      if (kpar %>% is.null) {
        kpar <- list(sigma = 1L, order = 1L, degree = 1L)
      }
      kf <- besseldot(unlist(kpar))
    } else if (kernel == 'anovadot') {
      if (kpar %>% is.null) {
        kpar <- list(sigma = 1L, degree = 1L)
      }
      kf <- anovadot(unlist(kpar))
    } else if (kernel == 'splinedot') {
      kf <- splinedot()
    }
    k_mat <- kernelMatrix(kernel = kf, x = t(dat))
    pca <- kpca(k_mat)                             # PCA, % variance explained
    pve <- seq_len(max(dims)) %>% map_chr(function(pc) {
      p <- as.numeric(eig(pca)[pc] / sum(eig(pca)) * 100L)
      paste0('KPC', pc, ' (', round(p, 2L), '%)')
    })
  }
  sig <- function(var, pc) {                       # p-val fn
    if (block %>% is.null | var %in% unblock | var == block) {
      mod <- lm(pca$x[, pc] ~ clin[[var]])
    } else {
      mod <- lm(pca$x[, pc] ~ clin[[var]] + clin[[block]])
    }
    ifelse(clin[[var]] %>% is.numeric,
           summary(mod)$coef[2L, 4L], anova(mod)[1L, 5L])
  }
  df <- crossing(Feature = colnames(clin),       # Melt
                      PC = paste0('PC', seq_len(n.pc))) %>%
    rowwise(.) %>%
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
# Add limits argument to scale_fill_gradientn to fix number to color mapping
# Some way to facet_grid?




