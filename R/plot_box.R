#' Box Plots by Sample
#'
#' This function displays each sample's omic data distribution as a box and
#' whisker plot.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples. Raw counts stored in \code{
#'   \link[edgeR]{DGEList}} or \code{\link[DESeq2]{DESeqDataSet}} objects are
#'   automatically extracted and transformed to the log2-CPM scale, with a
#'   warning.
#' @param group Optional character or factor vector of length equal to sample
#'   size. Numeric or logical vectors are silently coerced to factor. Levels are
#'   used to color box plots. If supplied, legend title defaults to "Group".
#'   Override this feature by passing a named list instead.
#' @param pal String specifying the color palette to use if \code{group} is not
#'   \code{NULL}. Options include \code{"ggplot"}, as well as the complete
#'   collection of \code{
#'   \href{https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html}{
#'   ggsci}} palettes, which can be identified by name (e.g., \code{"npg"},
#'   \code{"aaas"}, etc.). Alternatively, a character vector of colors with
#'   length equal to the number of levels in \code{group}.
#' @param ylab Optional label for y-axis.
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"right"}, \code{
#'   "left"}, \code{"top"}, \code{"bottom"}, \code{"bottomright"},
#'   \code{"bottomleft"}, \code{"topright"}, or \code{"topleft"}.
#' @param hover Show sample name by hovering mouse over data point? If \code{
#'   TRUE}, the plot is rendered in HTML and will either open in your browser's
#'   graphic display or appear in the RStudio viewer.
#'
#' @details
#' Box plots offer an intuitive way to visualize an omic data distribution. They
#' are common in quality control pipelines, and may be especially helpful when
#' contrasting pre- and post-normalization matrices. \code{plot_box} can
#' additionally be used to inspect for batch effects or associations with
#' phenotypic factors by using the \code{group} argument.
#'
#' @examples
#' # Box plots by sample
#' mat <- matrix(rnorm(500 * 10), nrow = 500, ncol = 10)
#' plot_box(mat)
#'
#' # Box plots by group
#' grp <- rep(c("A", "B"), each = 5)
#' plot_box(mat, group = grp)
#'
#' @export
#' @importFrom tidyr gather
#' @import dplyr
#' @import ggplot2
#'

plot_box <- function(dat,
                     group = NULL,
                       pal = 'd3',
                      type = NULL,
                      ylab = NULL,
                     title = NULL,
                    legend = 'outside',
                     hover = FALSE) {

  # Preliminaries
  if (!is.null(group)) {
    if (!is.list(group)) {
      group <- list(group)
    }
    if (length(group) > 1L) {
      stop('group cannot be a list of length > 1.')
    }
    if (length(group[[1]]) != ncol(dat)) {
      stop('group length must match number of samples in dat.')
    }
    if (length(unique(group[[1]])) == 1L) {
      warning('group is invariant.')
    }
    group[[1]] <- as.factor(group[[1]])
    if (is.null(names(group))) {
      names(group) <- 'Group'
    }
  }
  if (length(pal) == 1L & !is.color(pal)) {
    if (!pal %in% c('ggplot', 'npg', 'aaas', 'nejm', 'lancet', 'jco', 'ucscgb',
                    'd3', 'locuszoom', 'igv', 'uchicago', 'startrek',
                    'futurama', 'rickandmorty', 'simpsons', 'gsea')) {
      stop('pal not recognized.')
    }
  } else {
    if (!all(is.color(pal))) {
      stop('When passing multiple strings to pal, each must denote a valid ',
           'color in R.')
    }
    if (length(levels(group[[1]])) != length(pal)) {
      stop('When passing individual colors to pal, length(pal) must equal the ',
           'number of unique groups.')
    }
  }
  if (is.null(title)) {
    if (is.null(group)) {
      title <- 'Box Plots By Sample'
    } else {
      title <- paste('Box Plots By', names(group))
    }
  }
  if (!legend %in% c('right', 'left', 'top', 'bottom', 'bottomright',
                     'bottomleft', 'topright', 'topleft')) {
    stop('legend must be one of "right", "left", "top", "bottom", ',
         '"bottomright", "bottomleft", "topright", or "topleft".')
  }

  # Tidy data
  if (is.null(ylab)) {
    if (is(dat, 'DGEList') | is(dat, 'DESeqDataSet')) {
      ylab <- expression(log[2]*'-CPM Counts')
    } else if (is(dat, 'DESeqTransform')) {
      ylab <- 'Transformed Counts'
    } else {
      ylab <- 'Value'
    }
  }
  dat <- matrixize(dat)
  df <- gather(tbl_df(dat), Sample, Value)
  if (!is.null(group)) {
    df <- df %>%
      mutate(Group = rep(group[[1]], each = nrow(dat))) %>%
      arrange(Group) %>%
      mutate(Sample = factor(Sample, levels = unique(Sample)))
  }

  # Build plot
  p <- ggplot(df, aes(Sample, Value, text = Sample)) +
    labs(title = title, x = 'Sample', y = ylab) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
         axis.text.x = element_text(angle = 45L, hjust = 1L))
  if (!is.null(group)) {                         # Fill by group?
    p <- p + geom_boxplot(aes(fill = Group)) +
      scale_fill_manual(name = names(group),
                      values = colorize(pal, length((levels(group[[1]]))),
                                        var_type = 'Categorical'))
  } else {
    p <- p + geom_boxplot()
  }

  # Output
  gg_out(p, hover, legend)

}

# Fun fact: plotly won't display text for boxplots:
# https://community.plot.ly/t/boxplot-hoverinfo-text-not-display/1959


