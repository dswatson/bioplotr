#' Correlation Plot
#'
#' This function plots pairwise correlations between variables.
#'
#' @param dat A sample by feature data frame or matrix, e.g. of clinical
#'   variables. Non-numeric features are dropped with a warning.
#' @param method String specifying which correlation coefficient to compute.
#'   Must be one of \code{"pearson"}, \code{"kendall"}, or \code{"spearman"}.
#'   See \code{\link[stats]{cor}}.
#' @param use Optional character string giving a method for computing
#'   covariances in the presence of missing values. Must be one of \code{
#'   "everything"}, \code{"all.obs"}, \code{"complete.obs"}, \code{
#'   "na.or.complete"}, or \code{"pairwise.complete.obs"}.
#' @param alpha Optional significance threshold to impose on correlations.
#'   Those with \emph{p}-values (optionally adjusted) less than or equal to
#'   \code{alpha} are outlined in black.
#' @param p_adj Optional \emph{p}-value adjustment for multiple testing. Options
#'   include \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{
#'   "bonferroni"}, \code{"BH"}, \code{"BY"}, and \code{"fdr"}. See \code{
#'   \link[stats]{p.adjust}}.
#' @param lim Optional vector of length two defining lower and upper bounds for 
#'   the scale range. Default is observed extrema.
#' @param geom String specifying whether to visualize correlation coefficients
#'   as \code{"tile"} or \code{"circle"}.
#' @param label Print correlation coefficient over \code{geom}?
#' @param diag Include principal diagonal of the correlation matrix?
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"bottom"}, \code{"left"},
#'   \code{"top"}, \code{"right"}, \code{"bottomright"}, \code{"bottomleft"},
#'   \code{"topleft"}, or \code{"topright"}.
#' @param hover Show correlation coefficient by hovering mouse over the
#'   corresponding tile or circle? If \code{TRUE}, the plot is rendered in HTML
#'   and will either open in your browser's graphic display or appear in the
#'   RStudio viewer.
#' @param export Export correlation matrix? If \code{TRUE} and \code{alpha} is
#'   non-\code{NULL}, then the \emph{p}-value matrix will also be returned.
#'
#' @details
#' Correlation plots visualize the associations between numeric features. They
#' are a valuable tool in exploratory data analysis for biological experiments,
#' where they may help identify dependencies among clinical covariates, leading
#' to better omic models.
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
                      method = 'pearson',
                         use = 'everything',
                       alpha = NULL,
                       p_adj = NULL,
                         lim = NULL,
                        geom = 'tile',
                       label = FALSE,
                        diag = FALSE,
                       title = 'Correlation Plot',
                      legend = 'right',
                       hover = FALSE,
                      export = FALSE) {

  # Preliminaries
  if (ncol(dat) < 2L) {
    stop('dat must have at least two columns to generate a correlation matrix.')
  }
  if (!every(dat, is.numeric)) {
    dat <- keep(dat, is.numeric)
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
  dat <- as_tibble(dat)
  if (!method %in% c('pearson', 'kendall', 'spearman')) {
    stop('method must be one of "pearson", "kendall", or "spearman". See ?cor.')
  }
  uses <- c('everything', 'all.obs', 'complete.obs', 'na.or.complete',
            'pairwise.complete.obs')
  if (!use %in% uses) {
    stop('use must be one of ', stringify(uses, 'or'), '.')
  }
  if (!alpha %>% is.null) {
    if (alpha <= 0 | alpha >= 1) {
      stop('alpha must be numeric on (0, 1).')
    }
  }
  if (!p_adj %>% is.null) {
    p_adjes <- c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr')
    if (!p_adj %in% p_adjes) {
      stop('p_adj must be one of ', stringify(p_adjes, 'or'), 
           '. See ?p.adjust.')
    }
  }
  if (!geom %in% c('tile', 'circle')) {
    stop('geom must be either "tile" or "circle".')
  }
  loc <- c('bottom', 'left', 'top', 'right',
           'bottomright', 'bottomleft', 'topleft', 'topright')
  if (!legend %in% loc) {
    stop('legend must be one of ', stringify(loc, 'or'), '.')
  }

  # Tidy data
  mat <- cor(dat, method = method, use = use)
  mat[!lower.tri(mat)] <- NA_real_
  if (diag) {
    diag(mat) <- 1L
  }
  df <- mat %>%                                  # Melt correlation matrix
    as_tibble(.) %>%
    pivot_longer(everything(), names_to = 'x', values_to = 'Correlation') %>%
    mutate(y = rep(rownames(mat), each = nrow(mat))) %>%
    mutate(x = factor(x, levels = unique(x)),
           y = factor(y, levels = rev(unique(x))),
           Significant = FALSE) %>%
    select(x, y, Correlation, Significant) %>%
    na.omit(.)
  if (!alpha %>% is.null) {                      # Calculate p-value matrix?
    p_mat <- matrix(nrow = nrow(mat), ncol = ncol(mat))
    for (i in 2:ncol(p_mat)) {
      for (j in 1:(i - 1L)) {
        p_mat[i, j] <- cor.test(dat[[i]], dat[[j]],
                                method = method, use = use)$p.value
      }
    }
    p_val <- p_mat %>% keep(lower.tri(.))
    if (!p_adj %>% is.null) {
      p_val <- p.adjust(p_val, method = p_adj)
    }
    if (diag) {
      diag(p_mat) <- 1L
    }
    df <- df %>% mutate(Significant = if_else(p_val <= alpha, TRUE, FALSE))
  }
  if (export) {
    out <- list(Correlation = mat)
    if (!(alpha %>% is.null)) {
      out$p.value <- p_mat
    }
  }

  # Build plot
  p <- ggplot(df, aes(x, y)) +
    coord_equal() +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
         axis.text.x = element_text(angle = 45L, hjust = 1L))
  if (geom == 'tile') {
    p <- p + geom_tile(aes(fill = Correlation, color = Significant),
                       size = 1L, width = 0.9, height = 0.9) +
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
      guides(size = FALSE)
  }
  if (lim %>% is.null) {
    scale_color_gradientn(colors = brewer.pal(10L, 'RdBu'))
  } else {
    scale_color_gradientn(limits = lim, colors = brewer.pal(10L, 'RdBu'))
  }
  if (label) {
    p <- p + geom_text(aes(label = round(Correlation, 2)))
  }

  # Output
  gg_out(p, hover, legend)
  if (export) {
    return(out)
  }

}


