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
#' @param pal_group String specifying the color palette to use if \code{group}
#'   is non-\code{NULL}, or a vector of such strings with length equal to the
#'   number of vectors passed to \code{group}. Options include \code{"ggplot"}, 
#'   all qualitative color schemes available in \code{\href{
#'   https://bit.ly/2ipuEjn}{RColorBrewer}}, and the complete collection of 
#'   \code{\href{http://bit.ly/2bxnuGB}{ggsci}} palettes. Alternatively, a 
#'   character vector of colors with length equal to the cumulative number of 
#'   levels in \code{group}.
#' @param ylab Optional label for y-axis.
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"bottom"}, \code{"left"},
#'   \code{"top"}, \code{"right"}, \code{"bottomright"}, \code{"bottomleft"},
#'   \code{"topleft"}, or \code{"topright"}.
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

plot_box <- function(
  dat,
      group = NULL,
  pal_group = 'npg',
       type = NULL,
       ylab = NULL,
      title = 'Box Plot',
     legend = 'right',
      hover = FALSE
) {

  # Preliminaries
  if (!group %>% is.null) {
    group <- dat %>% format_features(group, 'Categorical')
    cols <- colorize(pal_group, var_type = 'Categorical',
                     n = length(levels(group[[1]])))
  }
  locations <- c('bottom', 'left', 'top', 'right',
                 'bottomright', 'bottomleft', 'topleft', 'topright')
  legend <- match.arg(legend, locations)

  # Tidy data
  if (ylab %>% is.null) {
    if (dat %>% is('DGEList') || dat %>% is('DESeqDataSet')) {
      ylab <- expression(log[2]*'-CPM Counts')
    } else if (dat %>% is('DESeqTransform')) {
      ylab <- 'Transformed Counts'
    } else {
      ylab <- 'Value'
    }
  }
  dat <- matrixize(dat)
  df <- as_tibble(dat) %>% gather('Sample', 'Value')
  if (!group %>% is.null) {
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
  if (!group %>% is.null) {                      # Fill by group?
    p <- p + geom_boxplot(aes(fill = Group)) +
      scale_fill_manual(name = names(group), values = cols)
  } else {
    p <- p + geom_boxplot()
  }

  # Output
  gg_out(p, hover, legend)

}

# Fun fact: plotly won't display text for boxplots:
# https://community.plot.ly/t/boxplot-hoverinfo-text-not-display/1959


