#' Q-Q Plot
#'
#' This function plots expected vs. observed \emph{p}-values following -log10
#' transform.
#'
#' @param dat Either a vector of \emph{p}-values, optionally named, or any
#'   object with a column for \emph{p}-values coercable to a data frame. Missing
#'   values are silently removed.
#' @param lambda Calculate genomic inflation factor? See Details.
#' @param title Optional plot title.
#' @param hover Show probe name by hovering mouse over data point? If \code{
#'   TRUE}, the plot is rendered in HTML and will either open in your browser's
#'   graphic display or appear in the RStudio viewer. Probe names are extracted
#'   from \code{dat}.
#'
#' @details
#' Q-Q plots are a common way to visually assess the applicability of a
#' statistical test to a given data set. If the black points deviate too sharply
#' from the red line, especially at low expected values of -log10(\emph{p}),
#' then it suggests a violation of the assumptions upon which the test was based.
#'
#' In addition, \code{plot_qq} optionally calculates the genomic inflation
#' factor \eqn{lambda}, defined as the ratio of the median of the observed
#' distribution of the test statistic to the expected median. Inflated \eqn{
#' lambda}-values (i.e., \eqn{lambda > 1}) are indicative of a high false
#' positive rate, possibly due to some systematic and unaccounted for bias in
#' the data.
#'
#' @examples
#' df <- data.frame(p.value = runif(10000))
#' plot_qq(df, lambda = TRUE)
#'
#' library(DESeq2)
#' dds <- makeExampleDESeqDataSet()
#' dds <- DESeq(dds)
#' res <- results(dds)
#' plot_qq(res)
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'

plot_qq <- function(dat,
                    lambda = FALSE,
                     title = NULL,
                     hover = FALSE) {

  # Preliminaries
  if (is.numeric(dat)) {
    dat <- data.frame(p.value = dat)
  } else {
    dat <- as.data.frame(dat)
  }
  p <- c('P.Value', 'PValue', 'pvalue', 'p.value')
  if (sum(p %in% colnames(dat)) == 1) {
    colnames(dat)[colnames(dat) %in% p] <- 'p.value'
  } else {
    stop('dat must include a p-value column. Recognized colnames for this vector ',
         'include "p.value", "P.Value", "PValue", and "pvalue". Make sure that dat',
         'includes exactly one such colname.')
  }
  dat <- na.omit(dat)
  if (nrow(dat) == 0L) {
    stop('No non-missing p-values.')
  }
  if (min(dat$p.value < 0L) || max(dat$p.value > 1L)) {
    stop('P-values must be on [0, 1].')
  }
  if (is.null(title)) {
    title <- 'Q-Q Plot'
  }

  # Tidy
  if (is.null(rownames(dat))) {
    dat <- dat %>% mutate(Probe = seq_len(nrow(dat)))
  } else {
    dat <- dat %>% mutate(Probe = rownames(dat))
  }
  df <- dat %>%
    mutate(Observed = -log10(sort(p.value)),
           Expected = -log10(ppoints(length(p.value)))) %>%
    select(Probe, Observed, Expected)
  if (lambda) {
    chisq <- qchisq(p = 1L - dat$p.value, df = 1L)
    lambda_val <- median(chisq) / qchisq(p = 0.5, df = 1L)
    lambda_lbl <- paste('lambda ==',  round(lambda_val, 2))
  }

  # Build plot
  size <- pt_size(df)
  p <- ggplot(df, aes(Expected, Observed, text = Probe)) +
    geom_point(size = size) +
    geom_abline(intercept = 0L, slope = 1L, color = 'red') +
    labs(title = title,
             x = expression('Expected'~-log[10](italic(p))),
             y = expression('Observed'~-log[10](italic(p)))) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (lambda) {
    p <- p + annotate('text', x = max(df$Expected), y = 0L, size = 5L,
                      hjust = 1L, label = lambda_lbl, parse = TRUE)
  }

  # Output
  gg_out(p, hover)

}

# lambda not eqn-ing right

