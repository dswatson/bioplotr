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
#' @param clin Data frame or matrix with rows corresponding to samples and
#'   columns to technical and/or biological features to test for associations
#'   with omic data. 
#' @param bivariate Test associations in isolation, or after adjusting for
#'   all remaining covariates? If \code{FALSE}, then \code{clin} is treated as 
#'   a design matrix against which each PC is sequentially regressed. See 
#'   Details.
#' @param parametric Compute \emph{p}-values using parametric association tests?
#'   If \code{FALSE}, rank-based alternatives are used instead. Must either be
#'   a single Boolean, in which case it applies to all tests, or Boolean vector
#'   of length equal to \code{ncol(clin)}. See Details.
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
#' @param stat Association statistic of choice. Currently supports \code{"p"} 
#'   (-log \emph{p}-values) and \code{"r2"} (\eqn{R^2}). Interpretations vary 
#'   depending on whether covariates are included. See Details.
#' @param alpha Optional significance threshold to impose on associations. 
#'   Those with \emph{p}-values (optionally adjusted) less than or equal to 
#'   \code{alpha} are outlined in black.
#' @param p_adj Optional \emph{p}-value adjustment for multiple testing. Options
#'   include \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{
#'   "bonferroni"}, \code{"BH"}, \code{"BY"}, and \code{"fdr"}. See \code{
#'   \link[stats]{p.adjust}}.
#' @param r_adj Return adjusted partial R-squared? Only relevant if \code{stat = 
#'   "r2"}. 
#' @param label Print association statistics over tiles?
#' @param lim Optional vector of length two defining lower and upper bounds for 
#'   the scale range. Default is observed extrema.
#' @param coord_equal Plot tiles of equal width and height?
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"bottom"}, \code{"left"},
#'   \code{"top"}, \code{"right"}, \code{"bottomright"}, \code{"bottomleft"},
#'   \code{"topleft"}, or \code{"topright"}.
#' @param hover Show \emph{p}-values by hovering mouse over tiles? If 
#'   \code{TRUE}, the plot is rendered in HTML and will either open in your 
#'   browser's graphic display or appear in the RStudio viewer.
#'
#' @details
#' If \code{bivariate = TRUE}, then association tests are performed between
#' each PC and each clinical covariate, optionally adjusting for a blocking 
#' variable (if \code{block} is non-\code{NULL}). If \code{bivariate = FALSE},
#' then all tests are partial association tests, in the sense that they adjust 
#' for all remaining covariates. 
#' 
#' When \code{bivariate = TRUE, block = NULL} and \code{parametric = TRUE}, 
#' significance is computed from Pearson correlation tests (for continuous 
#' features) or ANOVA \emph{F}-tests (for categorical features). When 
#' \code{parametric = FALSE}, significance is computed from rank-based 
#' alternatives, i.e. Spearman correlation tests (for continuous features) or 
#' Kruskal-Wallis tests (for categorical features). 
#' 
#' When \code{bivariate = FALSE} or \code{block} is non-\code{NULL}, 
#' significance is computed from partial correlation tests for continuous data 
#' (Pearson if \code{parametric = TRUE}, Spearman if \code{parametric = FALSE}) 
#' or repeated measures ANOVA \emph{F}-tests (under rank-transformation if 
#' \code{parametric = FALSE}). In all cases, the alternative hypothesis assumes
#' a monotonic relationship between variables.
#' 
#' Strength of association may be measured either by -log \emph{p}-values (if
#' \code{stat = "p"}) or R-squared (if \code{stat = "r2"}). The former may be
#' adjusted for multiple testing, while the latter can be adjusted for 
#' covariates. 
#'
#' A blocking variable may be provided if samples violate the assumption of 
#' independence, e.g. for studies in which subjects are observed at multiple 
#' time points. If a blocking variable is identified, it will be regressed out 
#' prior to testing for all variables except those explicitly exempted by the 
#' \code{unblock} argument. When supplying a blocking variable, be careful to 
#' consider potential collinearities in the data. For instance, clinical 
#' features may be invariant with respect to subject, while subject may be 
#' nested within other variables like batch or treatment group. The \code{block} 
#' and \code{unblock} arguments are intended to help parse out these 
#' relationships.
#' 
#' Numeric and categorical features are tested differently. To protect against
#' potential mistakes (e.g., one-hot encoding a Boolean variable), \code{
#' plot_drivers} automatically prints a data frame listing the class of each 
#' feature.
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
#'   select(cell, dex)
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
#' @importFrom rsq rsq.partial
#' @import foreach
#' @import dplyr
#' @import ggplot2
#'

plot_drivers <- function(
  dat,
  clin,
    bivariate = TRUE,
   parametric = TRUE,
        block = NULL,
      unblock = NULL,
       kernel = NULL,
         kpar = NULL,
          top = NULL,
         n_pc = 5L,
         stat = 'p',
        alpha = NULL,
        p_adj = NULL,
        r_adj = FALSE,
        label = FALSE,
          lim = NULL,
  coord_equal = TRUE,
        title = 'Variation By Feature',
       legend = 'right',
        hover = FALSE
) {

  # Preliminaries
  if (ncol(dat) < 3L) {
    stop('dat includes only ', ncol(dat), ' samples; need at least 3 for PCA.')
  }
  if (ncol(dat) != nrow(clin)) {
    stop('Number of columns in dat does not match number of rows in clin.')
  }
  if (length(parametric) > 1 & length(parametric) != ncol(clin)) {
    stop('parametric must be of length 1 or length ncol(clin).')
  }
  dat <- matrixize(dat)
  clin <- as_tibble(clin)
  if (!block %>% is.null) {
    if (!block %in% colnames(clin)) {
      stop(paste0('Column "', block, '" not found in clin.'))
    } else {
      if (unblock %>% is.null) {
        j_idx <- which(colnames(clin) != block)
        for (j in j_idx) {
          mm <- model.matrix(~ clin[[block]] + clin[[j]])
          if (!mm %>% is.fullrank) {
            stop('"', block,  '"', ' and "', colnames(clin)[j], '" are ', 
                 'perfectly collinear. Consider using the unblock argument.')
          }
        }
      } else {
        if (!all(unblock %in% colnames(clin))) {
          stop('The following column(s) not found in clin: ',
               stringify(unblock[!unblock %in% colnames(clin)], 'and'))
        }
      }
    }
    clin[[block]] <- as.character(clin[[block]])
  }
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
    p_adj <- match.arg(p_adj, p_adjes)
  }
  locations <- c('bottom', 'left', 'top', 'right',
                 'bottomright', 'bottomleft', 'topleft', 'topright')
  legend <- match.arg(legend, locations)
  tibble(Feature = colnames(clin),               # Be apprised
           Class = clin %>% map_chr(class)) %>%
    print(n = nrow(.))
  if (sum(is.na(clin)) > 0) {
    warning(sum(is.na(clin)), ' NA values detected in clin. Association tests ',
            'will proceed with pairwise deletion.')
  }

  # Compute PCs
  if (kernel %>% is.null) {
    pca <- prcomp(t(dat), rank. = n_pc)                        
    pve <- seq_len(n_pc) %>% map_chr(function(pc) {
      p <- round(pca$sdev[pc]^2L / sum(pca$sdev^2L) * 100L, 2L)
      paste0('PC', pc, '\n(', round(p, 2L), '%)')
    })
    pca <- pca$x
  } else {
    pca <- kpca_fn(dat, kernel, kpar, n_pc)                 
    pve <- seq_len(max(n_pc)) %>% map_chr(function(pc) {
      p <- as.numeric(eig(pca)[pc] / sum(eig(pca)) * 100L)
      paste0('KPC', pc, '\n(', round(p, 2L), '%)')
    })
    pca <- rotated(pca)
  }
  
  # List full models if bivariate = FALSE
  if (!bivariate) {
    f1_list <- lapply(seq_len(n_pc), function(pc) {
      if_else(parametric, lm(pca[, pc] ~ clin), lm(rank(pca[, pc]) ~ clin))
    })
  }
  
  # Association testing function
  association_test <- function(j, pc) {
    if (length(parametric) > 1) {
      parametric <- parametric[which(colnames(clin) == j)]
    }
    if (bivariate) {
      # The tmp data frame allows pairwise NA deletion
      tmp <- data.frame(x = clin[[j]], y = pca[, pc])
      if (!(block %>% is.null || j %in% unblock || j == block)) {
        tmp <- tmp %>% mutate(z = clin[[block]])
      }
      tmp <- na.omit(tmp)
      if (tmp$x %>% is.numeric) {
        # Regress out blocking effects if necessary
        if (!(block %>% is.null || j %in% unblock || j == block)) {
          x <- residuals(lm(x ~ z, data = tmp)) 
          y <- residuals(lm(y ~ z, data = tmp))
          tmp <- tmp %>% mutate(x = x, y = y)
        }
        # For continuous data, tests are Pearson or Spearman correlations
        # (optionally partial) depending on whether parametric = TRUE
        tst <- if_else(parametric, 
                       cor.test(tmp$x, tmp$y, method = 'pearson'),
                       cor.test(tmp$x, tmp$y, method = 'spearman'))
        p_value <- tst$p.value
        est <- if_else(stat == 'r', tst$estimate^2, p_value)
      } else {
        if (block %>% is.null || j %in% unblock || j == block) {
          # For categorical data with no blocking variable, options are 
          # ANOVA or Kruskal-Wallis test
          if (parametric) {
            f <- lm(y ~ x, data = tmp)
            p_value <- est <- anova(f)[1, 5]
          } else {
            f <- lm(rank(y) ~ x, data = tmp)
            p_value <- est <- kruskal.test(y ~ x, data = tmp)$p.value
          }
          if (stat == 'r') {
            est <- if_else(r_adj, summary(f)$adj.r.squared, summary(f)$r.squared)
          }
        } else {
          # When blocking variable is present, we perform repeated measures 
          # ANOVA, optionally on ranks (if parametric = FALSE)
          if (parametric) {
            f0 <- lm(y ~ z, data = tmp)
            f1 <- lm(y ~ z + x, data = tmp)
          } else {
            f0 <- lm(rank(y) ~ z, data = tmp)
            f1 <- lm(rank(y) ~ z + x, data = tmp)
          }
          p_value <- est <- anova(f0, f1)[2, 6]
          if (stat == 'r') {
            est <- rsq.partial(f1, f0, adj = r_adj, type = 'v')$partial.rsq
          }
        }
      }
    } else {
      # For multivariate tests, we run simple ANOVA, optionally on ranks 
      # (if parametric = FALSE)
      f1 <- f1_list[[pc]]
      f0 <- if_else(parametric, 
                    lm(pca[, pc] ~ select(clin, -j)),
                    lm(rank(pca[, pc]) ~ select(clin, -j)))
      p_value <- est <- anova(f0, f1)[2, 6]
      if (stat == 'r') {
        est <- rsq.partial(f1, f0, adj = r_adj, type = 'v')$partial.rsq
      }
    }
    out <- data.frame('Association' = est, 'p_value' = p_value)
    return(out)
  }
  
  # Tidy data
  df <- expand.grid(Feature = colnames(clin),    # Melt
                         PC = paste0('PC', seq_len(n_pc))) 
  tst_out <- foreach(i = seq_len(nrow(df)), .combine = rbind) %do% 
    association_test(df$Feature[i], df$PC[i])
  df <- cbind(df, tst_out)
  if (!p_adj %>% is.null) {
    df <- df %>% mutate(p_value = p.adjust(p_value, method = p_adj))
  }
  if (stat == 'p') {
    df <- df %>% mutate(Association = -log(p_value))
  }
  df <- df %>% mutate(Significant = if_else(p_value <= alpha, TRUE, FALSE))

  # Build plot
  if (stat == 'p') {
    leg_lab <- if_else(p_adj %in% c('fdr', 'BH', 'BY'),
                       expression(~-log(italic(q))), 
                       expression(~-log(italic(p))))
  } else {
    lim <- c(0, 1)
    leg_lab <- if_else(bivariate, expression(italic(R)^2), 
                       expression('Partial'~italic(R)^2))
  }

  p <- ggplot(df, aes(PC, Feature, fill = Association, text = Association,
                      color = Significant)) +
    geom_tile(size = 1L, width = 0.9, height = 0.9) +
    scale_color_manual(values = c('grey90', 'black')) +
    scale_x_discrete(labels = pve) +
    guides(color = FALSE) +
    labs(title = title, x = 'Principal Component') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (lim %>% is.null) {
    p <- p + scale_fill_gradientn(colors = c('white', 'pink', 'orange', 
                                             'red', 'darkred'),
                                  name = leg_lab)
  } else {
    p <- p + scale_fill_gradientn(colors = c('white', 'pink', 'orange', 
                                             'red', 'darkred'),
                                  name = leg_lab, limits = lim)
  }
  if (label) {
    p <- p + geom_text(aes(label = round(Association, 2L)))
  }
  if (coord_equal) {
    p <- p + coord_equal()
  }

  # Output
  gg_out(p, hover, legend)

}

 

# Some way to facet_grid? A "by" argument
# Use pData if available


