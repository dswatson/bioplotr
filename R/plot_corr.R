#' Correlation Plot
#'
#' This function plots pairwise correlations between variables.
#'
#' @param dat A sample by feature data frame or matrix, e.g. of clinical
#'   variables. Non-numeric features are dropped with a warning.
#' @param geom String specifying whether to visualize correlation coefficients
#'   as \code{"tile"} or \code{"circle"}.
#' @param label Print correlation coefficient over \code{geom}?
#' @param diag Include principal diagonal of the correlation matrix?
#' @param method String specifying which correlation coefficient to compute.
#'   Must be one of \code{"pearson"}, \code{"kendall"}, or \code{"spearman"}.
#'   See \code{\link[stats]{cor}}.
#' @param alpha Optional significance threshold to impose on correlations.
#'   Those with \emph{p}-values (optionally adjusted) less than or equal to
#'   \code{alpha} are outlined in black.
#' @param p.adj Optional \emph{p}-value adjustment for multiple testing. Options
#'   include \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{
#'   "bonferroni"}, \code{"BH"}, \code{"BY"}, and \code{"fdr"}. See \code{
#'   \link[stats]{p.adjust}}.
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"bottom"}, \code{"left"},
#'   \code{"top"}, \code{"right"}, \code{"bottomright"}, \code{"bottomleft"},
#'   \code{"topleft"}, or \code{"topright"}.
#' @param hover Show correlation coefficient by hovering mouse over the
#'   corresponding tile or circle? If \code{TRUE}, the plot is rendered in HTML
#'   and will either open in your browser's graphic display or appear in the
#'   RStudio viewer.
#'
#' @details
#' Correlation plots visualize the associations between numeric features. They
#' are a valuable tool in exploratory data analysis for biological experiments,
#' where they may help identify dependencies among clinical covariates, leading
#' to better omic models.
#'
#' @examples
#' mat <- matrix(rnorm(100), 10, 10)
#' plot_corr(mat)
#'
#' @export
#' @importFrom purrr some keep
#' @importFrom tidyr gather
#' @importFrom RColorBrewer brewer.pal
#' @import dplyr
#' @import ggplot2
#'

plot_corr <- function(dat,
                      geom = 'tile',
                     label = FALSE,
                      diag = FALSE,
                    method = 'pearson',
                     alpha = NULL,
                     p.adj = NULL,
                     title = NULL,
                    legend = 'right',
                     hover = FALSE) {

  # Preliminaries
  if (ncol(dat) < 2L) {
    stop('dat must have at least two columns to generate a correlation matrix.')
  }
  if (dat %>% some(!is.numeric)) {
    dat <- dat %>% keep(is.numeric)
    if (ncol(dat) < 2L) {
      stop('dat must have at least two numeric columns to generate a ',
           'correlation matrix.')
    } else {
      warning('Non-numeric variables have been detected and removed.')
    }
  }
  if (colnames(dat) %>% is.null) {
    colnames(dat) <- paste0('V', seq_len(ncol(dat)))
  }
  if (!geom %in% c('tile', 'circle')) {
    stop('geom must be either "tile" or "circle".')
  }
  if (!method %in% c('pearson', 'kendall', 'spearman')) {
    stop('method must be one of "pearson", "kendall", or "spearman". See ?cor.')
  }
  if (!(alpha %>% is.null)) {
    if (alpha <= 0 | alpha >= 1) {
      stop('alpha must be numeric on (0, 1).')
    }
  }
  if (!(p.adj %>% is.null)) {
    if (!p.adj %in% c('holm', 'hochberg', 'hommel',
                      'bonferroni', 'BH', 'BY', 'fdr')) {
      stop('p.adj must be one of "holm", "hochberg", "hommel", "bonferroni", ',
           '"BH", "BY", or "fdr". See ?p.adjust.')
    }
  }
  if (title %>% is.null) {
    title <- 'Correlation Plot'
  }
  loc <- c('bottom', 'left', 'top', 'right',
           'bottomright', 'bottomleft', 'topleft', 'topright')
  if (!legend %in% loc) {
    stop('legend must be one of ', stringify(loc, 'or'), '.')
  }

  # Tidy data
  mat <- cor(dat, method = method)
  mat[!lower.tri(mat)] <- NA
  if (diag) {
    diag(mat) <- 1L
  }
  df <- mat %>%                                  # Melt correlation matrix
    tbl_df(.) %>%
    gather('x', 'Correlation') %>%
    mutate(y = rep(rownames(mat), nrow(mat))) %>%
    mutate(x = factor(x, levels = unique(x)),
           y = factor(y, levels = rev(unique(x))),
           Significant = FALSE) %>%
    select(x, y, Correlation) %>%
    na.omit(.)
  if (!(alpha %>% is.null)) {                    # Calculate p-value matrix?
    p_mat <- matrix(nrow = nrow(mat), ncol = ncol(mat))
    for (i in 2L:ncol(p_mat)) {
      for (j in 1L:(i - 1L)) {
        p_mat[i, j] <- cor.test(dat[, i], dat[, j], method = method)$p.value
      }
    }
    p <- p_mat[lower.tri(p_mat)]
    if (!(p.adj %>% is.null)) {
      p <- p.adjust(p, method = p.adj)
    }
    if (diag) {
      diag(p_mat) <- 1L
    }
    df <- df %>% mutate(Significant = ifelse(p <= alpha, TRUE, FALSE))
  }

  # Build plot
  p <- ggplot(df, aes(x, y)) +
    coord_equal() +
    labs(x = NULL, y = NULL, title = title) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
         axis.text.x = element_text(angle = 45L, hjust = 1L))
  if (geom == 'tile') {
    p <- p + geom_tile(aes(fill = Correlation, color = Significant),
                       size = 1L, width = 0.9, height = 0.9) +
      scale_fill_gradientn(colors = brewer.pal(10L, 'RdBu')) +
      scale_color_manual(values = c('grey90', 'black')) +
      guides(color = FALSE)
  } else if (geom == 'circle') {
    p <- p + geom_point(data = df %>% filter(Significant),
                        aes(x, y, size = 1.25 * abs(Correlation)),
                        color = 'black', show.legend = FALSE) +
      geom_point(data = df %>% filter(!Significant),
                 aes(x, y, size = 1.25 * abs(Correlation)),
                 color = 'grey60', show.legend = FALSE) +
      geom_point(aes(color = Correlation, size = abs(Correlation))) +
      scale_color_gradientn(colors = brewer.pal(10L, 'RdBu')) +
      guides(size = FALSE)
  }
  if (label) {
    p <- p + geom_text(aes(label = round(Correlation, 2)))
  }

  # Output
  gg_out(p, hover, legend)

}

# NAs?
