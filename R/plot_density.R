#' Density Plots by Sample
#'
#' This function displays each sample's omic data distribution as a density
#' curve.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples. Raw counts stored in \code{
#'   \link[edgeR]{DGEList}} or \code{\link[DESeq2]{DESeqDataSet}} objects are
#'   automatically extracted and transformed to the log2-CPM scale, with a
#'   warning.
#' @param group Optional character or factor vector of length equal to sample
#'   size. Numeric or logical vectors are silently coerced to factor. Levels are
#'   used to color density curves. If supplied, legend title defaults to
#'   "Group". Override this feature by passing a named list instead.
#' @param pal String specifying the color palette to use if \code{group} is
#'   non-\code{NULL}. Options include \code{"ggplot"}, as well as the complete
#'   collection of \code{
#'   \href{https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html}{
#'   ggsci}} palettes, which can be identified by name (e.g., \code{"npg"},
#'   \code{"aaas"}, etc.). Alternatively, a character vector of colors with
#'   length equal to the number of levels in \code{group}.
#' @param xlab Optional label for x-axis.
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"right"}, \code{
#'   "left"}, \code{"top"}, \code{"bottom"}, \code{"bottomright"},
#'   \code{"bottomleft"}, \code{"topright"}, or \code{"topleft"}.
#' @param hover Show sample name by hovering mouse over data point? If \code{
#'   TRUE}, the plot is rendered in HTML and will either open in your browser's
#'   graphic display or appear in the RStudio viewer.
#'
#' @details
#' Density curves offer an intuitive way to visualize an omic data distribution.
#' They are common in quality control pipelines, and may be especially helpful
#' when contrasting pre- and post-normalization matrices. \code{plot_density}
#' can additionally be used to inspect for batch effects or associations with
#' phenotypic factors by using the \code{group} argument.
#'
#' @examples
#' # Density plots by sample
#' mat <- matrix(rnorm(500 * 10), nrow = 500, ncol = 10)
#' plot_density(mat)
#'
#' # Density plots by group
#' grp <- rep(c("A", "B"), each = 5)
#' plot_density(mat, group = grp)
#'
#' @seealso
#' \code{\link[limma]{plotDensities}}
#'
#' @export
#' @importFrom tidyr gather
#' @import dplyr
#' @import ggplot2
#'

plot_density <- function(dat,
                         group = NULL,
                           pal = 'd3',
                          xlab = NULL,
                         title = NULL,
                        legend = 'right',
                         hover = FALSE) {

  # Preliminaries
  if (!is.null(group)) {
    group <- format_features(dat, group, 'Categorical')
    cols <- colorize(pal = pal_group, var_type = 'Categorical',
                     n = length(levels(group[[1L]])))
  }
  if (is.null(title)) {
    if (is.null(group)) {
      title <- 'Density By Sample'
    } else {
      title <- paste('Density By', names(group))
    }
  }
  if (!legend %in% c('right', 'left', 'top', 'bottom',
                     'bottomright', 'bottomleft', 'topright', 'topleft')) {
    stop('legend must be one of "right", "left", "top", "bottom", ',
         '"bottomright", "bottomleft", "topright", or "topleft".')
  }

  # Tidy data
  if (is.null(xlab)) {
    if (is(dat, 'DGEList') | is(dat, 'DESeqDataSet')) {
      xlab <- expression(log[2]*'-CPM Counts')
    } else if (is(dat, 'DESeqTransform')) {
      xlab <- 'Transformed Counts'
    } else {
      xlab <- 'Value'
    }
  }
  dat <- matrixize(dat)
  df <- gather(tbl_df(dat), Sample, Value)
  if (!is.null(group)) {
    df <- df %>% mutate(Group = rep(group[[1L]], each = nrow(dat)))
  }

  # Build plot
  p <- ggplot(df, aes(Value, group = Sample, text = Sample)) +
    labs(title = title, x = xlab, y = 'Density') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (!is.null(group)) {                         # Color by group?
    p <- p + geom_path(stat = 'density', aes(color = Group)) +
      scale_color_manual(name = names(group), values = cols)
  } else {
    p <- p + geom_path(stat = 'density')
  }

  # Output
  gg_out(p, hover, legend)

}


