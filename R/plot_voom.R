#' Voom Plot
#'
#' This function visualizes the mean-variance relationship of count data after
#' applying the voom transformation.
#'
#' @param dat Raw count matrix, or an \code{\link[BioBase]{ExpressionSet}}
#'   object containing raw counts, or a \code{\link[edgeR]{DGEList}}. Data
#'   should be filtered prior to applying the \code{\link[limma]{voom}}
#'   transformation.
#' @param design Optional design matrix with rows corresponding to samples and
#'   columns to coefficients to be estimated.
#' @param lib.size Numeric vector containing total library sizes for each
#'   sample. If \code{NULL} and \code{dat} is a \code{DGEList}, then normalized
#'   library sizes are taken from counts. Otherwise library sizes are calculated
#'   from the columnwise count totals.
#' @param normalize.method Normalization method to be applied to the transformed
#'   counts. Choices are the same as for the \code{method} argument of
#'   \code{\link[limma]{normalizeBetweenArrays}} when the data is
#'   single-channel.
#' @param span Width of the LOWESS smoothing window as a proportion.
#' @param size Point size. 
#' @param alpha Point transparency.
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"bottom"}, \code{"left"},
#'   \code{"top"}, \code{"right"}, \code{"bottomright"}, \code{"bottomleft"},
#'   \code{"topleft"}, or \code{"topright"}.
#' @param hover Show probe name by hovering mouse over data point? If \code{
#'   TRUE}, the plot is rendered in HTML and will either open in your browser's
#'   graphic display or appear in the RStudio viewer. Probe names are extracted
#'   from \code{dat}.
#' @param ... Additional arguments to be passed to \code{\link[limma]{lmFit}}.
#'
#' @details
#' The \code{voom} function from the \code{limma} package offers a unique
#' approach to modeling count data. Rather than fitting negative binomial
#' regressions directly to genewise counts, \code{voom} applies a log2-CPM
#' transformation that renders the distribution approximately normal. A
#' (preliminary) linear model is fit to the transformed counts, from which a
#' mean-variance trend is estimated using LOWESS. Observation weights are then
#' computed as inverse predicted residual variance. These can be applied during
#' a final linear model fitting stage to counteract the heteroskedasticity
#' inherent to count data.
#'
#' The \code{voom} function optionally plots mean log2-CPM counts against
#' quarter-root residual variance. This plot is similar in principle to the
#' output of \code{\link{plot_mv}}, although the y-axis in that case represents
#' raw, not residual variance.
#'
#' @references
#' Law, C.W., Chen, Y., Shi, W., & Smyth, G.K. (2014).
#' \href{http://bit.ly/2sJxAg0}{voom: precision weights unlock linear model
#' analysis tools for RNA-seq read counts}. \emph{Genome Biology}, \strong{15}:
#' R29.
#'
#' @examples
#' library(limma)
#' mat <- matrix(rnbinom(1000 * 10, mu = 5, size = 5))
#' grp <- rep(c("ctl", "trt"), each = 5)
#' des <- model.matrix(~ grp)
#' plot_voom(mat, design = des)
#'
#' @seealso
#' \code{\link[limma]{voom}}
#'
#' @export
#' @importFrom limma getEAWP normalizeBetweenArrays lmFit
#' @importFrom purrr discard
#' @import dplyr
#' @import ggplot2
#' @importFrom ggsci pal_d3
#'

plot_voom <- function(
  dat,
            design = NULL,
          lib.size = NULL,
  normalize.method = 'none',
              span = 0.5,
              size = NULL,
             alpha = NULL,
             title = 'Voom Plot',
            legend = 'right',
             hover = FALSE, ...
) {

  # Preliminaries
  if (nrow(dat) < 2L) {
    stop('At least 2 probes required to fit a mean-variance trend.')
  }
  locations <- c('bottom', 'left', 'top', 'right',
                 'bottomright', 'bottomleft', 'topleft', 'topright')
  legend <- match.arg(legend, locations)

  # Tidy data
  if (dat %>% is('DGEList')) {
    if (design %>% is.null && length(dat$sample$group %>% levels) > 1L) {
      design <- model.matrix(~ group, data = dat$samples)
    }
    if (lib.size %>% is.null) {
      lib.size <- with(dat$samples, lib.size * norm.factors)
      dat <- dat$counts
    }
  } else {
    dat <- getEAWP(dat)$expr
  }
  if (design %>% is.null) {
    design <- matrix(1L, ncol(dat), 1L)
    rownames(design) <- colnames(dat)
    colnames(design) <- 'GrandMean'
  }
  if (lib.size %>% is.null) {
    lib.size <- colSums(dat)
  }
  y <- t(log2(t(dat + 0.5) / (lib.size + 1L) * 1e6L))      # Fit linear model on log2-CPM scale
  y <- normalizeBetweenArrays(y, method = normalize.method)
  fit <- lmFit(y, design, ...)
  if (fit$Amean %>% is.null) {
    fit$Amean <- rowMeans(y, na.rm = TRUE)
  }
  mu <- fit$Amean + mean(log2(lib.size + 1L)) - log2(1e6L)
  sigma <- sqrt(fit$sigma)
  allzero <- rowSums(dat) == 0L
  if (any(allzero)) {
    mu <- mu %>% discard(allzero)
    sigma <- sigma %>% discard(allzero)
  }
  lo <- lowess(mu, sigma, f = span)                        # Fit LOWESS curve
  df <- tibble(Probe = rownames(fit),
                  Mu = mu,
               Sigma = sigma) %>%
    arrange(Mu) %>%
    mutate(lfit = lo[['y']])

  # Build plot
  if (size %>% is.null) {
    size <- pt_size(df)
  }
  if (alpha %>% is.null) {
    alpha <- pt_alpha(df)
  }
  suppressWarnings(
    p <- ggplot(df) +
      geom_point(aes(Mu, Sigma, text = Probe), size = size, alpha = alpha) +
      geom_path(aes(Mu, lfit, color = 'LOWESS'), size = 0.5) +
      scale_color_manual(name = 'Curve', values = pal_d3()(1)) +
      labs(title = title,
               x = expression('Mean'~log[2]*'-CPM'),
               y = expression(sqrt(sigma))) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
  )

  # Output
  gg_out(p, hover, legend)

}


