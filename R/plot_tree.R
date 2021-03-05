#' Tree Plots
#'
#' This function plots tree structures such as dendrograms and regression trees.
#'
#' @param dat An object representing a dendrogram (e.g. of class \code{
#'   \link[stats]{hclust}}), regression tree (e.g. of class \code{\link[rpart]{
#'   rpart}}), or other tree-like structure.
#' @param labels Plot labels? This refers to items for dendrograms, and break
#'   points for regression or classification trees.
#' @param leaf_labels Plot leaf labels? Only relevant for regression and
#'   classification trees.
#' @param rotate Rotate plot 90 degrees?
#' @param title Optional plot title.
#'
#' @details
#' This function is essentially a wrapper for tools in the \code{\href{
#' https://github.com/andrie/ggdendro}{ggdendro}} package, with a few amendments
#' and extra options. See the package documentation for more info.
#'
#' @examples
#' hc <- hclust(dist(USArrests))
#' plot_tree(hc)
#'
#' @export
#' @importFrom ggdendro is.dendro dendro_data
#' @import dplyr
#' @import ggplot2
#'

plot_tree <- function(
  dat,
       labels = TRUE,
  leaf_labels = TRUE,
       rotate = FALSE,
        title = NULL
) {

  # Preliminaries
  dat_class <- if_else(dat %>% inherits('dendro'), dat$class, class(dat))
  angle <- if (dat_class %in% c('dendrogram', 'hclust')) {
    if_else(rotate, 0, 90)
  } else {
    if_else(rotate, 90, 0)
  }
  hjust <- if (dat_class %in% c('dendrogram', 'hclust')) {
    if_else(rotate, 0, 1)
  } else {
    0.5
  }

  # Tidy
  if (!dat %>% is.dendro) {
    dat <- dendro_data(dat, type = 'rectangle')
  }

  # Build plot
  p <- ggplot() +
    geom_blank() +
    theme_bw() +
    labs(title = title, x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = angle, hjust = 1L)) +
    theme(axis.text.y = element_text(angle = angle, hjust = 1L))
  if (dat$class == 'tree') {
    p <- p + geom_segment(data = segment(dat),
                          aes_string(size = 'n', x = 'x', y = 'y',
                                     xend = 'xend', yend = 'yend'))
  } else {
    p <- p + geom_segment(data = segment(dat),
                          aes_string(x = 'x', y = 'y',
                                     xend = 'xend', yend = 'yend'))
  }
  if (leaf_labels && !(dat$leaf_labels %>% is.null)) {
    p <- p + geom_text(data = leaf_label(dat),
                       aes_string(x = 'x', y = 'y', label = 'label'),
                       hjust = hjust, angle = angle)
  }
  if (labels && dat$class %in% c('tree', 'rpart')) {
    p <- p + geom_text(data = dat$labels,
                       aes_string(x = 'x', y = 'y', label = 'label'),
                       size = 3L, vjust = -0.5)
  }
  if (rotate) {
    p <- p + coord_flip() +
      scale_y_reverse()
    if (labels) {
      p <- p + scale_x_continuous(breaks = seq_along(dat$labels$label),
                                  labels = dat$labels$label,
                                position = 'top')
    }
  } else {
    p <- p + scale_y_continuous()
    if (labels) {
      p <- p + scale_x_continuous(breaks = seq_along(dat$labels$label),
                                  labels = dat$labels$label)
    }
  }

  # Output
  print(p)

}


# Optionally color clades and/or items?
# Add auto-text sizing
