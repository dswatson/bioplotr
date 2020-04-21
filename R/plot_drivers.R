#' Plot Drivers of Omic Variation
#'
#' This function visualizes the strength of associations between the principal
#' components of an omic data matrix and a set of biological and/or technical
#' features.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples. It is strongly recommended that data be
#'   filtered and normalized prior to plotting. Raw counts stored in 
#'   \code{\link[edgeR]{DGEList}} or \code{\link[DESeq2]{DESeqDataSet}} objects 
#'   are automatically extracted and transformed to the log2-CPM scale, with a
#'   warning.
#' @param clin Data frame or matrix with rows correponding to samples and
#'   columns to technical and/or biological features to test for associations
#'   with omic data.
#' @param parametric Compute \emph{p}-values using parametric association tests?
#'   If \code{FALSE}, rank-based alternatives are used instead. See Details.
#' @param block String specifying the name of the column in which to find the
#'   blocking variable, should one be accounted for. See Details.
#' @param unblock Column name(s) of one or more features for which the 
#'   \code{block} covariate should not be applied, if one was supplied. See 
#'   Details.
#' @param kernel The kernel generating function, if using KPCA. Options include
#'   \code{"rbfdot"}, \code{"polydot"}, \code{"tanhdot"}, \code{"vanilladot"}, 
#'   \code{"laplacedot"}, \code{"besseldot"}, \code{"anovadot"}, and 
#'   \code{"splinedot"}. To run normal PCA, set to \code{NULL}. See Details.
#' @param kpar A named list of arguments setting parameters for the kernel
#'   function. Only relevant if \code{kernel} is not \code{NULL}. See Details.
#' @param top Optional number (if > 1) or proportion (if < 1) of most variable
#'   probes to be used for PCA.
#' @param n_pc Number of principal components to include in the figure.
#' @param label Print association statistics over tiles?
#' @param alpha Optional significance threshold to impose on associations. 
#'   Those with \emph{p}-values (optionally adjusted) less than or equal to 
#'   \code{alpha} are outlined in black.
#' @param p_adj Optional \emph{p}-value adjustment for multiple testing. Options
#'   include \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{
#'   "bonferroni"}, \code{"BH"}, \code{"BY"}, and \code{"fdr"}. See \code{
#'   \link[stats]{p.adjust}}.
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"bottom"}, \code{"left"},
#'   \code{"top"}, \code{"right"}, \code{"bottomright"}, \code{"bottomleft"},
#'   \code{"topleft"}, or \code{"topright"}.
#' @param hover Show \emph{p}-values by hovering mouse over tiles? If 
#'   \code{TRUE}, the plot is rendered in HTML and will either open in your 
#'   browser's graphic display or appear in the RStudio viewer.
#'
#' @details
#' Strength of association is measured by -log \emph{p}-values, optionally
#' adjusted for multiple testing. When \code{parametric = TRUE}, significance
#' is computed from Pearson correlation tests (for continuous features) or 
#' ANOVA \emph{F}-tests (for categorical features). When \code{parametric =
#' FALSE}, significance is computed from rank-based alternatives, i.e. Spearman 
#' correlation tests (for continuous features) or Kruskal-Wallis tests (for 
#' categorical features). 
#'
#' An optional blocking variable may be provided if samples violate the
#' assumption of independence, e.g. for studies in which subjects are observed
#' at multiple time points. If a blocking variable is identified, it will be
#' regressed out prior to testing for all variables except those explicitly
#' exempted by the \code{unblock} argument. Significance is then computed from
#' partial correlation tests for continuous data (Pearson if \code{parametric = 
#' TRUE}, Spearman if \code{parametric = FALSE}) or repeated measures ANOVA 
#' \emph{F}-tests (under rank-transformation if \code{parametric = FALSE}).
#'
#' When supplying a blocking variable, be careful to consider potential
#' confounding effects. For instance, features like sex and age are usually
#' nested within subject, while subject may be nested within other variables
#' like batch or treatment group. The \code{block} and \code{unblock} arguments
#' are designed to help parse out these relationships.
#' 
#' Numeric and categorical features are tested differently. To protect against
#' potential mistakes (e.g., one-hot encoding a Boolean variable), 
#' \code{plot_drivers} automatically prints a data frame listing the class of
#' each feature.
#'
#' If \code{kernel} is non-\code{NULL}, then KPCA is used instead of PCA. See
#' \code{\link{plot_kpca}} for more info. Details on kernel functions and their
#' input parameters can be found in \code{kernlab::\link[kernlab]{dots}}.
#'
#' @examples
#' library(SummarizedExperiment)
#' library(edgeR)
#' library(dplyr)
#' data(airway)
#' cnts <- assay(airway)
#' keep <- rowSums(cpm(cnts) > 1) >= 4
#' mat <- cpm(cnts[keep, ], log = TRUE)
#' clin <- colData(airway) %>%
#'   as_tibble(.) %>%
#'   select(Run, cell, dex)
#' plot_drivers(mat, clin)
#'
#' @seealso
#' \code{\link{plot_pca}}, \code{\link{plot_kpca}}
#'
#' @export
#' @importFrom limma is.fullrank
#' @importFrom purrr map_chr
#' @importFrom kernlab eig
#' @importFrom kernlab rotated
#' @import dplyr
#' @import ggplot2
#'

plot_drivers <- function(dat,
                         clin,
                    parametric = TRUE,
                         block = NULL,
                       unblock = NULL,
                        kernel = NULL,
                          kpar = NULL,
                           top = NULL,
                          n_pc = 5L,
                         label = FALSE,
                         alpha = NULL,
                         p_adj = NULL,
                         title = 'Variation By Feature',
                        legend = 'right',
                         hover = FALSE) {

  # Preliminaries
  if (ncol(dat) < 3L) {
    stop('dat includes only ', ncol(dat), ' samples; need at least 3 for PCA.')
  }
  if (ncol(dat) != nrow(clin)) {
    stop('Number of columns in dat does not match number of rows in clin.')
  }
  dat <- matrixize(dat)
  clin <- as_tibble(clin)
  if (!(block %>% is.null)) {
    if (!block %in% colnames(clin)) {
      stop(paste0('Column "', block, '" not found in clin.'))
    } else {
      if (unblock %>% is.null) {
        j_idx <- which(colnames(clin) != block)
        for (j in j_idx) {
          mm <- model.matrix(~ clin[[block]] + clin[[j]])
          if (!mm %>% is.fullrank) {
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
  clin[[block]] <- as.character(clin[[block]])
  if (!top %>% is.null) {                          # Filter by variance?
    dat <- var_filt(dat, top, robust = FALSE)
  }
  if (n_pc > max(nrow(dat), ncol(dat))) {
    stop('n_pc cannot exceed max(nrow(dat), ncol(dat))')
  }
  if (!alpha %>% is.null) {
    if (alpha <= 0L | alpha >= 1L) {
      stop('alpha must be numeric on (0, 1).')
    }
  } else {
    alpha <- 0L
  }
  if (!p_adj %>% is.null) {
    p_adjes <- c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr')
    if (!p_adj %in% p_adjes) {
      stop('p_adj must be one of ', stringify(p_adjes, 'or'), 
           '. See ?p.adjust.')
    }
  }
  loc <- c('bottom', 'left', 'top', 'right',
           'bottomright', 'bottomleft', 'topleft', 'topright')
  if (!legend %in% loc) {
    stop('legend must be one of ', stringify(loc, 'or'), '.')
  }
  tibble(Feature = colnames(clin),               # Be apprised
           Class = clin %>% map_chr(class)) %>%
    print(n = nrow(.))

  # Compute PCs
  if (kernel %>% is.null) {
    pca <- prcomp(t(dat))                        
    pve <- seq_len(n_pc) %>% map_chr(function(pc) {
      p <- round(pca$sdev[pc]^2L / sum(pca$sdev^2L) * 100L, 2L)
      paste0('PC', pc, '\n(', round(p, 2L), '%)')
    })
    pca <- pca$x
  } else {
    pca <- kpca_fn(dat, kernel, kpar)                 
    pve <- seq_len(max(n_pc)) %>% map_chr(function(pc) {
      p <- as.numeric(eig(pca)[pc] / sum(eig(pca)) * 100L)
      paste0('KPC', pc, '\n(', round(p, 2L), '%)')
    })
    pca <- rotated(pca)
  }
  
  # P-value function
  sig <- function(j, pc) {
    if (clin[[j]] %>% is.numeric) {
      if (block %>% is.null | j %in% unblock | j == block) {
        x <- clin[[j]]
        y <- pca[, pc]
      } else {
        x <- residuals(lm(clin[[j]] ~ clin[[block]], data = df)) 
        y <- residuals(lm(pca[, pc] ~ clin[[block]], data = df))
      }
      if (parametric) {
        p_val <- cor.test(x, y, method = 'pearson')$p.value
      } else {
        p_val <- cor.test(x, y, method = 'spearman')$p.value
      }
    } else {
      if (block %>% is.null | j %in% unblock | j == block) {
        if (parametric) {
          p_val <- anova(lm(pca[, pc] ~ clin[[j]]))[1, 5]
        } else {
          p_val <- kruskal.test(pca[, pc] ~ clin[[j]])$p.value 
        }
      } else {
        if (parametric) {
          f0 <- lm(pca[, pc] ~ clin[[block]])
          f1 <- lm(pca[, pc] ~ clin[[block]] + clin[[j]])
        } else {
          f0 <- lm(rank(pca[, pc]) ~ clin[[block]])
          f1 <- lm(rank(pca[, pc]) ~ clin[[block]] + clin[[j]])
        }
        p_val <- anova(f0, f1)[2, 6]
      }
    }
    return(p_val)
  }
  
  # Tidy data
  df <- expand.grid(Feature = colnames(clin),    # Melt
                         PC = paste0('PC', seq_len(n_pc))) %>%
    rowwise(.) %>%
    mutate(Association = sig(Feature, PC)) %>%   # Populate
    ungroup(.)
  if (!p_adj %>% is.null) {
    df <- df %>% mutate(Association = p.adjust(Association, method = p_adj))
  }
  df <- df %>% 
    mutate(Significant = if_else(Association <= alpha, TRUE, FALSE),
           Association = -log(Association))

  # Build plot
  if (!p_adj %>% is.null & p_adj %in% c('fdr', 'BH', 'BY')) {
    leg_lab <- expression(~-log(italic(q)))
  } else {
    leg_lab <- expression(~-log(italic(p)))
  }
  p <- ggplot(df, aes(PC, Feature, fill = Association, text = Association,
                      color = Significant)) +
    geom_tile(size = 1L, width = 0.9, height = 0.9) +
    coord_equal() +
    scale_fill_gradientn(colors = c('white', 'pink', 'orange', 'red', 'darkred'),
                           name = leg_lab) +
    scale_color_manual(values = c('grey90', 'black')) +
    scale_x_discrete(labels = pve) +
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

# Allow spline fits? 
# Optionally summarise associations with R^2? MSE?
# Allow parametric tests for some assocations and nonparametric for others?
# Fit multivariate models?
# Fages & Ferrari, 2014: https://link.springer.com/article/10.1007/s11306-014-0647-9
# Add limits argument to scale_fill_gradientn to fix number to color mapping
# Some way to facet_grid? A by argument
# Use pData if available


