#' Density Plots by Sample
#'
#' This function displays each sample's omic data distribution as a density curve.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples.
#' @param group Optional character or factor vector of length equal to sample size.
#'   Numeric or logical vectors will be silently coerced to factor. Levels are used
#'   to color density curves. If supplied, legend title defaults to "Group". Override
#'   this feature by passing a named list instead.
#' @param xlab Optional label for x-axis.
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"outside", "bottomleft",
#'   "bottomright", "topleft",} or \code{"topright"}.
#' @param hover Show sample name by hovering mouse over data point? If \code{TRUE},
#'   the plot is rendered in HTML and will either open in your browser's graphic
#'   display or appear in the RStudio viewer.
#'
#' @details
#' Density curves are an intuitive way to visualize an omic data distribution. They
#' are especially helpful when contrasting pre- and post-normalization matrices.
#' \code{plot_density} may additionally be used to inspect for batch effects
#' or associations with phenotypic factors by using the \code{group} argument.
#'
#' @examples
#' # Simulated data
#' mat <- matrix(rnorm(1000 * 5), nrow = 1000, ncol = 5)
#' plot_density(mat)
#'
#' # Real data: raw counts
#' data(airway)
#' library(edgeR)
#' cnts <- assay(airway)
#' keep <- rowSums(cpm(cnts) > 1) >= 4           # Filter out underexpressed genes
#' y <- DGEList(cnts[keep, ])                    # Create DGEList object
#' plot_density(y, group = colData(airway)$dex)
#'
#' # Real data: log2-CPM transformed counts
#' y <- calcNormFactors(y)
#' y <- cpm(y, log = TRUE)                       # Apply log2-CPM transformation
#' plot_density(y, group = colData(airway)$dex)
#'
#' @seealso
#' \code{\link[limma]{plotDensities}}
#'
#' @export
#' @importFrom DESeq2 counts
#' @importFrom SummarizedExperiment assay
#' @importFrom limma getEAWP
#' @importFrom tidyr gather
#' @importFrom ggsci scale_color_d3
#' @import dplyr
#' @import ggplot2
#' @importFrom plotly ggplotly
#'

plot_density <- function(dat,
                         group = NULL,
                          xlab = NULL,
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
  if (is.null(title)) {
    if (is.null(group)) {
      title <- 'Density By Sample'
    } else {
      title <- paste('Density By', names(group))
    }
  }
  if (!legend %in% c('outside', 'bottomleft', 'bottomright', 'topleft', 'topright')) {
    stop('legend must be one of "outside", "bottomleft", "bottomright" ',
         '"topleft", or "topright".')
  }

  # Tidy data
  if (is(dat, 'DGEList')) {
    dat <- dat$counts
    if (is.null(ylab)) {
      xlab <- 'Raw Counts'
    }
  } else if (is(dat, 'DESeqDataSet')) {
    dat <- counts(dat)
    if (is.null(ylab)) {
      xlab <- 'Raw Counts'
    }
  } else if (is(dat, 'DESeqTransform')) {
    dat <- assay(dat)
    if (is.null(ylab)) {
      xlab <- 'Transformed Counts'
    }
  } else {
    dat <- getEAWP(dat)$expr
    keep <- rowSums(is.finite(dat)) == ncol(dat)
    dat <- dat[keep, , drop = FALSE]
    if (is.null(ylab)) {
      xlab <- 'Value'
    }
  }
  df <- gather(tbl_df(dat), Sample, Value)
  if (!is.null(group)) {
    df <- df %>% mutate(Group = rep(group[[1]], each = nrow(dat)))
  }

  # Build plot
  p <- ggplot(df, aes(Value, group = Sample, text = Sample)) +
    labs(title = title, x = xlab, y = 'Density') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (!is.null(group)) {                         # Color by group?
    p <- p + geom_path(stat = 'density', aes(color = Group)) +
      guides(color = guide_legend(title = names(group))) +
      scale_color_d3()
  } else {
    p <- p + geom_path(stat = 'density')
  }
  p <- locate_legend(p, legend)

  # Output
  gg_out(p, hover, legend)

}


