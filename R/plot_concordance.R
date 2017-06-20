#' Concordance Plot
#'
#' This function plots the pairwise concordance between categorical features.
#'
#' @param dat A sample by feature data frame or matrix, e.g. of clinical
#'   variables or patient cluster assignments. All columns are converted to
#'   factors with a warning, if possible.
#' @param label Print association statistic over tiles?
#' @param diag Include principal diagonal of the concordance matrix? Only
#'   advisable if \code{method = "MI"}.
#' @param method String specifying which measure of association to
#'   compute. Currently supports \code{"MI"}, \code{"fisher"}, and \code{
#'   "chisq"}. See Details.
#' @param alpha Optional significance threshold to impose on association
#'   statistics. Those with \emph{p}-values (optionally adjusted) less than or
#'   equal to \code{alpha} are outlined in black.
#' @param p.adj Optional \emph{p}-value adjustment for multiple testing. Options
#'   include \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{
#'   "bonferroni"}, \code{"BH"}, \code{"BY"}, and \code{"fdr"}. See \code{
#'   \link[stats]{p.adjust}}.
#' @param B Number of permutations to generate when computing \emph{p}-values
#'   via Monte Carlo simulation. See Details.
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"right"}, \code{
#'   "left"}, \code{"top"}, \code{"bottom"}, \code{"topright"}, \code{
#'   "topleft"}, \code{"bottomright"}, or \code{"bottomleft"}.
#' @param hover Show correlation coefficient by hovering mouse over the
#'   corresponding tile or circle? If \code{TRUE}, the plot is rendered in HTML
#'   and will either open in your browser's graphic display or appear in the
#'   RStudio viewer.
#'
#' @details
#' Concordance plots visualize the
#'
#' When \code{method = "fisher"}, concordance is measured by the negative
#' logarithm of the test's \emph{p}-value. If this cannot be computed in the
#' available workspace, then it is estimated via simulation using \code{B}
#' Monte Carlo permutations.
#'
#' When \code{method = "chisq"}, concordance is measured by the Pearson
#' chi-squared statistic. This is easier to compute than Fisher exact tests, but
#' may be inappropriate if BLAH
#'
#' When \code{method = "MI"}, concordance is measured by the mutual information
#' statistic. If \code{alpha} is non-\code{NULL}, then \code{p}-values are
#' estimated via simulation using \code{B} Monte Carlo permutations.
#'
#' @examples
#' df <- data.frame(A = sample(1:2, 20, replace = TRUE),
#'                  B = sample(1:3, 20, replace = TRUE),
#'                  C = sample(1:4, 20, replace = TRUE))
#' plot_concordance(df)
#'
#' @export
#' @importFrom purrr map map_lgl
#' @importFrom tidyr gather
#' @importFrom infotheo mutinformation natstobits
#' @import dplyr
#' @import ggplot2
#'

plot_concordance <- function(dat,
                             label = FALSE,
                              diag = FALSE,
                            method = 'MI',
                             alpha = NULL,
                             p.adj = NULL,
                                 B = 2000,
                             title = NULL,
                            legend = 'right',
                             hover = FALSE) {

  # Preliminaries
  if (ncol(dat) < 2) {
    stop('dat must have at least two columns to generate a concordance matrix.')
  }
  if (!all(map_lgl(dat, is.numeric))) {
    dat <- dat[, map_lgl(dat, is.numeric)]
    if (ncol(dat) < 2) {
      stop('dat must have at least two numeric columns to generate a ',
           'concordance matrix.')
    } else {
      warning('Non-numeric variables have been detected and removed.')
    }
  }
  if (is.null(colnames(dat))) {
    colnames(dat) <- paste0('V', seq_len(ncol(dat)))
  }
  if (!geom %in% c('tile', 'circle')) {
    stop('geom must be either "tile" or "circle".')
  }
  if (!method %in% c('MI', 'fisher', 'chisq')) {
    stop('method must be one of "MI", "fisher", or "chisq".')
  }
  if (!is.null(alpha)) {
    if (alpha <= 0 | alpha >= 1) {
      stop('alpha must be numeric on (0, 1).')
    }
  }
  if (!is.null(p.adj)) {
    if (!p.adj %in% c('holm', 'hochberg', 'hommel',
                      'bonferroni', 'BH', 'BY', 'fdr')) {
      stop('p.adj must be one of "holm", "hochberg", "hommel", "bonferroni", ',
           '"BH", "BY", or "fdr". See ?p.adjust.')
    }
  }
  if (is.null(title)) {
    title <- 'Concordance Plot'
  }

  # Tidy data
  mat <- matrix(nrow = ncol(dat), ncol = ncol(dat),
                dimnames = list(colnames(dat), colnames(dat)))
  for (i in 2L:ncol(dat)) {
    for (j in 1L:(i - 1L)) {
      tmp <- na.omit(dat[, c(i, j)])
      if (method == 'MI') {
        mat[i, j] <- natstobits(mutinformation(tmp[[1L]], tmp[[2L]]))
      } else if (method == 'fisher') {
        p <- try(fisher.test(tmp[[1L]], tmp[[2L]], workspace = 2e8L)$p.value,
                 silent = TRUE)
        if (is(p, 'try-error')) {
          mat[i, j] <- fisher.test(tmp[[1L]], tmp[[2L]], simulate.p.value = TRUE,
                                   B = B)$p.value
        } else {
          mat[i, j] <- p
        }
        p_mat <- mat
        mat[lower.tri(mat)] <- -log(mat[lower.tri(mat)])
      } else if (method == 'chisq') {
        mat[i, j] <- as.numeric(chisq.test(tmp[[1L]], tmp[[2L]])$statistic)
      }
    }
  }
  df <- data.frame(mat) %>%
    gather(x, Association) %>%
    mutate(y = rep(rownames(mat), nrow(mat))) %>%
    mutate(x = factor(x, levels = unique(x)),
           y = factor(y, levels = rev(unique(x))),
           Significant = FALSE) %>%
    select(x, y, Association) %>%
    na.omit()
  if (!is.null(alpha)) {                         # p-value matrix?
    if (!method == 'fisher') {
      p_mat <- matrix(nrow = nrow(mat), ncol = ncol(mat))
      for (i in 2L:ncol(p_mat)) {
        for (j in 1L:(i - 1L)) {
          tmp <- na.omit(dat[, c(i, j)])
          if (method == 'MI') {
            null <- numeric(length = B)
            for (i in seq_len(B)) {
              x <- sample(tmp[[1L]], nrow(tmp))
              y <- sample(tmp[[2L]], nrow(tmp))
              null[i] <- natstobits(mutinformation(x, y))
            }
            p_mat[i, j] <- sum(null >= mat[i, j]) / B
          } else if (method == 'chisq') {
            p_mat[i, j] <- chisq.test(tmp[[1L]], tmp[[2L]])$p.value
          }
        }
      }
    }
    p <- p_mat[lower.tri(p_mat)]
    if (!is.null(p.adjust)) {
      p <- p.adjust(p, method = p.adj)
    }
    df <- df %>% mutate(Significant = ifelse(p <= alpha, TRUE, FALSE))
  }

  # Build plot
  if (method == 'MI') {
    leg.txt <- 'Mutual Information (Bits)'
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
    labs(x = NULL, y = NULL, title = title) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
         axis.text.x = element_text(angle = 45L, hjust = 1L))
  if (label) {
    p <- p + geom_text(aes(label = round(Association, 2L)))
  }

  # Output
  gg_out(p, hover, legend)

}


# NAs?
# diag?

