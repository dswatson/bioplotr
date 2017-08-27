#' Concordance Plot
#'
#' This function plots the pairwise concordance between categorical features.
#'
#' @param dat A sample by feature data frame or matrix, e.g. of clinical
#'   variables or patient cluster assignments. All columns are converted to
#'   factors with a warning, if possible.
#' @param method String specifying which measure of association to
#'   compute. Currently supports \code{"fisher"}, \code{"chisq"}, and \code{
#'   "MI"}. See Details.
#' @param alpha Optional significance threshold to impose on association
#'   statistics. Those with \emph{p}-values (optionally adjusted) less than or
#'   equal to \code{alpha} are outlined in black.
#' @param p.adj Optional \emph{p}-value adjustment for multiple testing. Options
#'   include \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{
#'   "bonferroni"}, \code{"BH"}, \code{"BY"}, and \code{"fdr"}. See \code{
#'   \link[stats]{p.adjust}}.
#' @param sim.p Calculate \emph{p}-values via Monte Carlo simulation?
#' @param B Number of replicates or permutations to generate when computing
#'   \emph{p}-values.
#' @param label Print association statistic over tiles?
#' @param diag Include principal diagonal of the concordance matrix? Only
#'   advisable if \code{method = "MI"}.
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
#' When \code{method = "fisher"}, concordance is measured by the negative
#' logarithm of the test's \emph{p}-value. If \code{sim.p = FALSE}, then the
#' function will attempt to calculate an exact \emph{p}-value. If this cannot be
#' executed in the available workspace, or if \code{sim.p = TRUE}, then \emph{
#' p}-values are estimated via Monte Carlo simulation with \code{B} replicates.
#'
#' When \code{method = "chisq"}, concordance is measured by the Pearson
#' chi-squared statistic. This test is based on several assumptions that may not
#' be met in practice (see \href{http://bit.ly/2te2u5h}{Wikipedia} for a quick
#' overview). When one or several of these assumptions are violated, more
#' accurate \emph{p}-values can be estimated via Monte Carlo simulation with
#' \code{B} replicates.
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
#' @examples
#' df <- data.frame(A = sample.int(2, 20, replace = TRUE),
#'                  B = sample.int(3, 20, replace = TRUE),
#'                  C = sample.int(4, 20, replace = TRUE))
#' plot_concordance(df)
#'
#' @export
#' @importFrom purrr some rerun map_dbl keep
#' @importFrom tidyr gather
#' @importFrom infotheo mutinformation natstobits
#' @import dplyr
#' @import ggplot2
#'

plot_concordance <- function(dat,
                             method = 'fisher',
                              alpha = NULL,
                              p.adj = NULL,
                              sim.p = FALSE,
                                  B = 2000,
                              label = FALSE,
                               diag = FALSE,
                              title = NULL,
                             legend = 'right',
                              hover = FALSE,
                             export = FALSE) {

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
  if (!method %in% c('fisher', 'chisq', 'MI')) {
    stop('method must be one of "fisher", "chisq", or "MI".')
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
    title <- 'Concordance Plot'
  }
  loc <- c('bottom', 'left', 'top', 'right',
           'bottomright', 'bottomleft', 'topleft', 'topright')
  if (!legend %in% loc) {
    stop('legend must be one of ', stringify(loc, 'or'), '.')
  }


  # Tidy Data
  mat <- matrix(nrow = ncol(dat), ncol = ncol(dat),
                dimnames = list(colnames(dat), colnames(dat)))
  for (i in 2:ncol(dat)) {
    for (j in 1:(i - 1L)) {
      tmp <- dat[, c(i, j)] %>%
        as_tibble(.) %>%
        na.omit(.)
      if (method == 'MI') {
        mat[i, j] <- mutinformation(tmp[[1]], tmp[[2]]) %>% natstobits(.)
      } else if (method == 'fisher') {
        if (sim.p) {
          p <- fisher.test(tmp[[1]], tmp[[2]],
                           simulate.p.value = TRUE, B = B)$p.value
        } else {
          p <- try(fisher.test(tmp[[1]], tmp[[2]], workspace = 2e8L)$p.value,
                   silent = TRUE)
          if (p %>% is('try-error')) {
            mat[i, j] <- fisher.test(tmp[[1]], tmp[[2]],
                                     simulate.p.value = TRUE, B = B)$p.value
          } else {
            mat[i, j] <- p
          }
        }
        p_mat <- mat
        mat[lower.tri(mat)] <- -log(mat[lower.tri(mat)])
      } else if (method == 'chisq') {
        mat[i, j] <- chisq.test(tmp[[1]], tmp[[2]])$statistic
      }
    }
  }
  df <- mat %>%                                  # Melt concordance matrix
    as_tibble(.) %>%
    gather('x', 'Association') %>%
    mutate(y = rep(rownames(mat), nrow(mat))) %>%
    mutate(x = factor(x, levels = unique(x)),
           y = factor(y, levels = rev(unique(x))),
           Significant = FALSE) %>%
    select(x, y, Association, Significant) %>%
    na.omit(.)
  if (!(alpha %>% is.null)) {                    # Calculate p-value matrix?
    if (method != 'fisher') {
      p_mat <- matrix(nrow = nrow(mat), ncol = ncol(mat))
      for (i in 2:ncol(p_mat)) {
        for (j in 1:(i - 1L)) {
          tmp <- dat[, c(i, j)] %>% na.omit(.)
          if (method == 'MI') {
            null <- B %>%
              rerun(x = tmp[[1]][sample.int(nrow(tmp))],
                    y = tmp[[2]]) %>%
              map_dbl(~ mutinformation(.x$x, .x$y)) %>%
              natstobits(.)
            p_mat[i, j] <- (sum(null >= mat[i, j]) + 1L) / (B + 1L)
          } else if (method == 'chisq') {
            if (sim.p) {
              p_mat[i, j] <- chisq.test(tmp[[1]], tmp[[2]],
                                        simulate.p.value = TRUE, B = B)$p.value
            } else {
              p_mat[i, j] <- chisq.test(tmp[[1]], tmp[[2]])$p.value
            }
          }
        }
      }
    }
    p <- p_mat %>% keep(lower.tri(.))
    if (!(p.adjust %>% is.null)) {
      p <- p.adjust(p, method = p.adj)
    }
    df <- df %>% mutate(Significant = ifelse(p <= alpha, TRUE, FALSE))
  }
  if (export) {
    out <- list(Concordance = mat)
    if (!(alpha %>% is.null)) {
      out$p.values <- p_mat
    }
  }

  # Build Plot
  if (method == 'MI') {
    leg.txt <- 'Mutual Information\n (Bits)'
  } else if (method == 'fisher') {
    leg.txt <- expression(~-log(italic(p)))
  } else if (method == 'chisq') {
    leg.txt <- expression(chi^2)
  }
  p <- ggplot(df, aes(x, y, fill = Association, text = Association,
                      color = Significant)) +
    geom_tile(size = 1L, width = 0.9, height = 0.9) +
    coord_equal() +
    scale_fill_gradientn(colors = c('white', 'pink', 'orange', 'red', 'darkred'),
                           name = leg.txt) +
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
# Rand index, other stats from Dudoit & Fridlyand, 2002?


