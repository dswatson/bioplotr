#' Concordance Plot
#'
#' This function plots the pairwise concordance between categorical features.
#'
#' @param dat A sample by feature data frame or matrix, e.g. of clinical
#'   variables or patient cluster assignments. All columns are converted to
#'   factors with a warning, if possible.
#' @param method String specifying which measure of association to
#'   compute. Currently supports \code{"chisq"}, \code{USP}, \code{"fisher"}, 
#'   and \code{"MI"}. See Details.
#' @param alpha Optional significance threshold to impose on association
#'   statistics. Those with \emph{p}-values (optionally adjusted) less than or
#'   equal to \code{alpha} are outlined in black.
#' @param p_adj Optional \emph{p}-value adjustment for multiple testing. Options
#'   include \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{
#'   "bonferroni"}, \code{"BH"}, \code{"BY"}, and \code{"fdr"}. See \code{
#'   \link[stats]{p.adjust}}.
#' @param sim_p Calculate \emph{p}-values via Monte Carlo simulation? Only 
#'   relevant if \code{method = "chisq"} or \code{method = "fisher"}.
#' @param B Number of replicates or permutations to sample when computing
#'   \emph{p}-values.
#' @param label Print association statistic over tiles?
#' @param diag Include principal diagonal of the concordance matrix? Only
#'   advisable if \code{method = "MI"}.
#' @param pal_tiles String specifying the color palette to use for heatmap
#'   tiles. Options include the complete collection of \code{\href{
#'   https://bit.ly/2n7D6tF}{viridis}} palettes, as well as all sequential and
#'   divergent color schemes available in \code{\href{
#'   https://bit.ly/2ipuEjn}{RColorBrewer}}. Alternatively, a character vector 
#'   of at least two colors.
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"bottom"}, \code{"left"},
#'   \code{"top"}, \code{"right"}, \code{"bottomright"}, \code{"bottomleft"},
#'   \code{"topleft"}, or \code{"topright"}.
#' @param hover Show association statistic by hovering mouse over the
#'   corresponding tile or circle? If \code{TRUE}, the plot is rendered in HTML
#'   and will either open in your browser's graphic display or appear in the
#'   RStudio viewer.
#' @param export Export concordance matrix? If \code{TRUE} and \code{alpha} is
#'   non-\code{NULL}, then the \emph{p}-value matrix will also be returned.
#'
#' @details
#' Concordance plots visualize associations between categorical features. They
#' are useful when evaluating the dependencies between clinical factors and/or
#' patient clusters.
#'
#' When \code{method = "chisq"}, concordance is measured by the Pearson
#' chi-squared statistic. This test is based on several assumptions that may not
#' be met in practice (see \href{http://bit.ly/2te2u5h}{Wikipedia} for a quick
#' overview). When one or several of these assumptions are violated, more
#' accurate \emph{p}-values can be estimated via Monte Carlo simulation with
#' \code{B} replicates.
#' 
#' When \code{method = "USP"}, concordance is measured by the negative logarithm
#' of the \emph{p}-value of a \emph{U}-statistic permutation test (Berrett et 
#' al., 2021), which is minimax optimal under mild conditions. The test uses 
#' \code{B} permutations.
#' 
#' When \code{method = "fisher"}, concordance is measured by the negative
#' logarithm of the test's \emph{p}-value. If \code{sim_p = FALSE}, then the
#' function will attempt to calculate an exact \emph{p}-value. If this cannot be
#' executed in the available workspace, or if \code{sim_p = TRUE}, then \emph{
#' p}-values are estimated via Monte Carlo simulation with \code{B} replicates.
#'
#' When \code{method = "MI"}, concordance is measured by the mutual information
#' statistic. If \code{alpha} is non-\code{NULL}, then \emph{p}-values are
#' estimated via permutation testing with \code{B} permutations.
#'
#' @return
#' If \code{export = TRUE}, a list with up to two elements:
#' \itemize{
#'   \item The concordance matrix, computed via the chosen \code{method}.
#'   \item The matrix of \emph{p}-values (optionally adjusted), if \code{alpha}
#'   is non-\code{NULL}.
#' }
#'
#' @references
#' Berrett, T.B., Kontoyiannis, I. & Samworth, R. (2021).
#' \href{https://arxiv.org/pdf/2001.05513.pdf}{Optimal rates for independence
#' testing via \emph{U}-statistic permutation tests}.
#' \emph{Ann. Statist.}.
#' 
#' @examples
#' df <- data.frame(A = sample.int(2, 20, replace = TRUE),
#'                  B = sample.int(3, 20, replace = TRUE),
#'                  C = sample.int(4, 20, replace = TRUE))
#' plot_concordance(df)
#'
#' @export
#' @importFrom USP USP.test
#' @importFrom purrr some rerun map_dbl keep
#' @importFrom tidyr pivot_longer
#' @importFrom infotheo mutinformation natstobits
#' @import dplyr
#' @import ggplot2
#'

plot_concordance <- function(
  dat,
     method = 'chisq',
      alpha = 0.05,
      p_adj = NULL,
      sim_p = FALSE,
          B = 1999L,
      label = FALSE,
       diag = FALSE,
  pal_tiles = 'PiRdBr',
      title = 'Concordance Plot',
     legend = 'right',
      hover = FALSE,
     export = FALSE
) {

  # Preliminaries
  p <- ncol(dat)
  if (p < 2L) {
    stop('dat must have at least two columns to generate a concordance matrix.')
  }
  if (some(dat, is.numeric)) {
    for (j in seq_along(dat)) {
      if (dat[[j]] %>% is.numeric) {
        if (!all.equal(dat[[j]], as.integer(dat[[j]]))) {
          dat <- dat[[-j]]
        }
      }
    }
    if (ncol(dat) < 2L) {
      stop('dat must have at least two non-numeric columns to generate a ',
           'concordance matrix.')
    }
    if (p != ncol(dat)) {
      warning('Continuous features have been detected and removed.')
    }
  }
  if (colnames(dat) %>% is.null) {
    colnames(dat) <- paste0('V', seq_len(ncol(dat)))
  }
  methods <- c('chisq', 'USP', 'fisher', 'MI')
  method <- match.arg(method, methods)
  if (!alpha %>% is.null) {
    if (alpha <= 0 | alpha >= 1) {
      stop('alpha must be numeric on (0, 1).')
    }
  }
  if (!p_adj %>% is.null) {
    p_adjes <- c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr')
    p_adj <- match.arg(p_adj, p_adjes)
  }
  pal_cols <- colorize(pal_tiles, var_type = 'Continuous')
  locations <- c('bottom', 'left', 'top', 'right',
                 'bottomright', 'bottomleft', 'topleft', 'topright')
  legend <- match.arg(legend, locations)


  # Tidy Data
  mat <- p_mat <- matrix(nrow = ncol(dat), ncol = ncol(dat),
                         dimnames = list(colnames(dat), colnames(dat)))
  for (i in 2:ncol(dat)) {
    for (j in 1:(i - 1L)) {
      tmp <- dat[, c(i, j)] %>%
        as_tibble(.) %>%
        na.omit(.)
      if (method == 'chisq') {
        test <- chisq.test(tmp[[1]], tmp[[2]])
        stat <- test$statistic
        p <- test$p.value
      } else if (method == 'USP') {
        freq <- table(tmp[[1]], tmp[[2]])
        test <- USP.test(freq)
        p <- test$p.value
        stat <- -log10(p)
      } else if (method == 'fisher') {
        if (sim_p) {
          p <- fisher.test(tmp[[1]], tmp[[2]],
                           simulate.p.value = TRUE, B = B)$p.value
        } else {
          p <- try(fisher.test(tmp[[1]], tmp[[2]], workspace = 2e8L)$p.value,
                   silent = TRUE)
          if (p %>% is('try-error')) {
            p <- fisher.test(tmp[[1]], tmp[[2]],
                             simulate.p.value = TRUE, B = B)$p.value
          } 
          stat <- -log10(p)
        }
      } else if (method == 'MI') {
        stat <- mutinformation(tmp[[1]], tmp[[2]]) %>% natstobits(.)
        null <- B %>%
          rerun(x = tmp[[1]][sample.int(nrow(tmp))],
                y = tmp[[2]]) %>%
          map_dbl(~ mutinformation(.x$x, .x$y)) %>%
          natstobits(.)
        p <- (sum(null >= stat) + 1L) / (B + 1L)
      }
      mat[i, j] <- stat
      p_mat[i, j] <- p
    }
  }
  df <- mat %>%                                  # Melt concordance matrix
    as_tibble(.) %>%
    pivot_longer(everything(), names_to = 'x', values_to = 'Association') %>%
    mutate(y = rep(rownames(mat), nrow(mat))) %>%
    mutate(x = factor(x, levels = unique(x)),
           y = factor(y, levels = rev(unique(x))),
           Significant = FALSE) %>%
    select(x, y, Association, Significant) %>%
    na.omit(.)
  if (!alpha %>% is.null) {
    p_val <- p_mat %>% keep(lower.tri(.))
  }
  if (!p_adj %>% is.null) {
    p_val <- p.adjust(p, method = p_adj)
  }
  df <- df %>% mutate(Significant = if_else(p_val <= alpha, TRUE, FALSE))
  if (export) {
    out <- list(
      'Concordance' = mat,
         'p.values' = p_mat,
            'p_adj' = p_adj
    )
  }

  # Build Plot
  if (method == 'chisq') {
    leg.txt <- expression(chi^2)
  } else if (method %in% c('USP', 'fisher')) {
    leg.txt <- expression(~-log[10](italic(p)))
  } else if (method == 'MI') {
    leg.txt <- 'Mutual Information\n(Bits)'
  }
  p <- ggplot(df, aes(x, y, fill = Association, text = Association,
                      color = Significant)) +
    geom_tile(size = 1L, width = 0.9, height = 0.9) +
    coord_equal() +
    scale_fill_gradientn(colors = pal_cols, name = leg.txt) +
    scale_color_manual(values = c('grey90', 'black')) +
    guides(color = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
         axis.text.x = element_text(angle = 45L, hjust = 1L))
  if (label) {
    p <- p + geom_text(aes(label = round(Association, 2L)))
  }

  # Output
  gg_out(p, hover, legend)
  if (export) {
    return(out)
  }

}


# NAs?
# diag?
# Rand index, adjusted Rand index, other stats from Dudoit & Fridlyand, 2002?


