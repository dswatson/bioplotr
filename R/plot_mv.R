#' Mean-Variance Plot
#'
#' This function visualizes the mean-variance relationship of omic data before
#' or after modeling.
#'
#' @param dat Either an omic data matrix with rows corresponding to probes and
#'   columns to samples, or an object of class \code{\link[limma]{MArrayLM}} as
#'   created by a call to \code{limma}'s \code{\link[limma]{lmFit}} or
#'   \code{\link[limma]{eBayes}} functions. If \code{dat} is a matrix, then data
#'   are presumed to be normalized prior to visualization. Any matrix-like object
#'   that can be processed by \code{\link[limma]{getEAWP}} is also acceptable.
#' @param trans Data transformation to be applied to probewise standard deviations
#'   if \code{dat} is a matrix. Must be one of either \code{"log"} or \code{"sqrt"}.
#'   See Details.
#' @param zero.weights Should spots with zero or negative weights be plotted?
#' @param ptsize Size of data points in the plot.
#' @param main Optional plot title.
#' @param legend Legend position. Must be one of \code{"outside", "bottomleft",
#'   "bottomright", "topleft",} or \code{"topright"}. Only relevant if \code{dat}
#'   is an \code{MArrayLM} object created by a call to \code{eBayes} with
#'   \code{trend = TRUE}.
#' @param hover Show probe name by hovering mouse over data point? If \code{TRUE},
#'   the plot is rendered in HTML and will either open in your browser's graphic
#'   display or appear in the RStudio viewer. Probe names are extracted from
#'   \code{dat}.
#'
#' @details
#' Mean-variance plots are a quick and easy way to visualize the relationship
#' between the first two moments of probewise data distributions. When used prior to
#' modeling, they may help better understand the internal structure of the data
#' and inspect for potential outliers. When used after modeling, they can be
#' useful in evaluating the assumptions of the regression. \code{plot_mv} fits
#' a smooth curve to the points using a generalized additive model with a cubic
#' regression spline. See \code{mgcv::\link[mgcv]{gam}} for more details.
#'
#' If \code{dat} is a matrix, then the appropriate data transformation for probewise
#' standard deviations must be specified. \code{trans = "log"} is recommended for
#' microarrays or any other platform in which data are approximately log-normally
#' distributed. \code{trans = "sqrt"} is recommended for sequencing and count data.
#'
#' If \code{dat} is an \code{MArrayLM} object, then the y-axis will be log2
#' transformed probewise residual standard deviations. If standard errors were
#' moderated with \code{eBayes}, then prior variance will also be plotted as either
#' a horizontal line or a smooth curve, depending on whether a global or
#' intensity-dependent prior was used. See \code{\link[limma]{squeezeVar}} for more
#' details.
#'
#' @examples
#' mat <- matrix(rnorm(1000 * 10), nrow = 1000, ncol = 10)
#' plot_mv(mat, trans = "log")
#'
#' library(limma)
#' grp <- rep(c("ctl", "trt"), each = 5)
#' des <- model.matrix(~ grp)
#' fit <- eBayes(lmFit(mat, des))
#' plot_mv(fit)
#'
#' @seealso
#' \code{\link[limma]{plotSA}} \code{\link[vsn]{meanSdPlot}}
#'
#' @export
#' @importFrom limma getEAWP
#' @importFrom matrixStats rowSds
#' @import dplyr
#' @importFrom purrr map_lgl
#' @import ggplot2
#' @importFrom plotly ggplotly
#'

plot_mv <- function(dat,
                    trans = NULL,
             zero.weights = FALSE,
                   ptsize = 0.25,
                     main = NULL,
                   legend = 'outside',
                    hover = FALSE) {

  # Preliminaries
  if (!is(dat, 'MArrayLM')) {
    dat <- getEAWP(dat)
    wts <- dat$weights
    dat <- dat$expr
    if (is.null(trans) || !trans %in% c('log', 'sqrt')) {
      stop('trans must be specified as either "log" or "sqrt". The former is ',
           'recommended for microarrays or any other platform in which data are ',
           'approximately log-normally distributed. The latter is recommended for ',
           'sequencing and count data.')
    }
  } else {
    wts <- dat$weights
    if (!is.null(trans)) {
      warning('trans is not evaluated when dat is of class MArrayLM. Residual standard ',
              'deviations are always plotted under log2 transform.')
    }
  }
  if (!is.null(rownames(dat))) {
    probes <- rownames(dat)
  } else {
    probes <- seq_len(nrow(dat))
  }
  if (!zero.weights && !is.null(wts)) {
    allzero <- rowSums(wts > 0, na.rm = TRUE) == 0
    dat <- dat[!allzero, ]
  }
  if (is.null(main)) {
    main <- 'Mean-Variance Plot'
  }
  if (!legend %in% c('outside', 'bottomleft', 'bottomright', 'topleft', 'topright')) {
    stop('legend must be one of "outside", "bottomleft", "bottomright" ',
         '"topleft", or "topright".')
  }

  # Tidy data
  if (is.matrix(dat)) {             # Unmodeled data
    if (trans == 'log') {
      Sigma <- log2(rowSds(dat))
      ylab <- expression('log'[2]*(sigma))
    } else {
      Sigma <- sqrt(rowSds(dat))
      ylab <- expression(sqrt(sigma))
    }
    df <- data_frame(Probe = probes,
                        Mu = rowMeans(dat),
                     Sigma = Sigma)
  } else {                          # Modeled data
    df <- data_frame(Probe = probes,
                        Mu = dat$Amean,
                     Sigma = log2(dat$sigma))
    if ('s2.prior' %in% names(dat)) {
      df <- df %>% mutate(Prior = log2(dat$s2.prior))
    }
    if (length(dat$s2.prior) > 1) {
      df2 <- max(dat$df.prior)
      s2 <- dat$sigma^2 / dat$s2.prior
      pdn <- pf(s2, df1 = dat$df.residual, df2 = df2)
      pup <- pf(s2, df1 = dat$df.residual, df2 = df2, lower.tail = FALSE)
      FDR <- p.adjust(2 * pmin(pdn, pup), method = 'BH')
      df <- df %>% mutate(Outlier = map_lgl(seq_len(nrow(dat)), function(i) {
        ifelse(FDR[i] <= 0.05, TRUE, FALSE)
      }))
    }
    ylab <- expression('log'[2]*(sigma))
  }

  # Build plot
  p <- ggplot(df) +
    labs(title = main, x = expression(mu), y = ylab) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if ('Outlier' %in% colnames(df) && any(df$Outlier)) {
    suppressWarnings(
      p <- p + geom_point(aes(Mu, Sigma, text = Probe, color = Outlier),
                          size = ptsize, alpha = 0.25)
    )
  } else {
    suppressWarnings(
      p <- p + geom_point(aes(Mu, Sigma, text = Probe),
                          size = ptsize, alpha = 0.25)
    )
  }
  if ('Prior' %in% colnames(df)) {  # Plot prior
    p <- p + geom_smooth(aes(Mu, Sigma, color = 'GAM fit'),
                         method = 'gam', formula = y ~ s(x, bs = 'cs'),
                         size = 0.5, se = FALSE)
    if (length(dat$s2.prior) == 1) {
      p <- p + geom_abline(aes(color = 'Prior'), slope = 0, intercept = dat$s2.prior)
    } else {
      p <- p + geom_smooth(aes(Mu, Prior, color = 'Prior'),
                           method = 'gam', formula = y ~ s(x, bs = 'cs'),
                           size = 0.5, se = FALSE)
    }
    p <- p + scale_color_manual(name = 'Curves', values = c('red', 'blue')) +
      guides(col = guide_legend(reverse = TRUE))
  } else {
    p <- p + geom_smooth(aes(Mu, Sigma), method = 'gam', formula = y ~ s(x, bs = 'cs'),
                         size = 0.5, se = FALSE)
  }
  if (legend == 'bottomleft') {     # Locate legend
    p <- p + theme(legend.justification = c(0.01, 0.01),
                   legend.position = c(0.01, 0.01))
  } else if (legend == 'bottomright') {
    p <- p + theme(legend.justification = c(0.99, 0.01),
                   legend.position = c(0.99, 0.01))
  } else if (legend == 'topleft') {
    p <- p + theme(legend.justification = c(0.01, 0.99),
                   legend.position = c(0.01, 0.99))
  } else if (legend == 'topright') {
    p <- p + theme(legend.justification = c(0.99, 0.99),
                   legend.position = c(0.99, 0.99))
  }

  # Output
  if (!hover) {
    print(p)
  } else {
    p <- ggplotly(p, tooltip = 'text', height = 600, width = 600)
    print(p)
  }

}


# Extend to DESeq2 model objects
# Add option for sigma ~ rank(mu)

