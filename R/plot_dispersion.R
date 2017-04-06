#' Plot Dispersion Estimates
#'
#' This function plots the mean-dispersion relationship for a given count matrix.
#'
#' @param dat A \code{\link[edgeR]{DGEList}} object or \code{\link[DESeq2]{
#'   DESeqDataSet}}. If normalization factors or dispersions have not already been
#'   estimated for \code{dat}, they will be internally computed using the appropriate
#'   \code{edgeR} or \code{DESeq2} functions. For the former, a design matrix may be
#'   required to calculate the genewise adjusted profile log-likelihoods. See
#'   Details.
#' @param design Optional study design matrix. This will be extracted from \code{dat}
#'   if available and need not be set explicitly in the call. If provided, however,
#'   \code{design} will override the relevant slots or components of \code{dat}.
#' @param trans Data transformation to be applied to genewise dispersion estimates.
#'  Must be one of \code{"log"} or \code{"sqrt"}.
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"outside", "bottomleft",
#'   "bottomright", "topleft",} or \code{"topright"}.
#' @param hover Show probe name by hovering mouse over data point? If \code{TRUE},
#'   the plot is rendered in HTML and will either open in your browser's graphic
#'   display or appear in the RStudio viewer. Gene names are extracted from
#'   \code{dat}.
#'
#' @details
#' Count data in omic studies are often presumed to follow a negative binomial
#' distribution, which may be uniquely identified by its mean and dispersion
#' parameters. For RNA-seq pipelines that rely on negative binomial generalized linear
#' models (GLMs), such as \code{edgeR} and \code{DESeq2}, estimating genewise
#' dispersions is an essential step in the model fitting process. There are rarely
#' sufficient replicates to reliably infer these values independently for each gene,
#' however, which is why both packages use empirical Bayes methods to pool information
#' across genes.
#'
#' Details vary across the two pipelines, which is why output will differ for \code{
#' DGEList} objects and \code{DESeqDataSet}s. \code{edgeR} begins by computing a common
#' dispersion parameter for the entire dataset, rendered by \code{plot_dispersion} as
#' a blue horizontal line; then fits a trend line to the maximized negative binomial
#' likelihoods, represented by an orange curve; and finally calculates posterior
#' estimates, depicted by black points, using a weighted empirical Bayes likelihood
#' procedure. Be aware that likelihood maximization methods for \code{DGEList}
#' objects vary depending on whether or not a model matrix is supplied. See
#' \code{\link[edgeR]{estimateDisp}} for more details. A thorough explication of the
#' statistical theory behind this pipeline can be found in the original papers by
#' the packages authors: Robinson & Smyth (2007), Robinson & Smyth (2008), and
#' McCarthy et al. (2012).
#'
#' \code{DESeq2} also fits a trend line through likelihood estimates of genewise
#' dispersions, depicted by orange points in the \code{plot_dispersion} output.
#' Posterior values are calculated following regression toward a log-normal prior with
#' mean equal to the predicted value from the trended fit and variance equal to the
#' difference between the observed variance of the log dispersion estimates and the
#' expected sampling variance. These maximum a posteriori values are colored blue,
#' while outliers, defined as genes with log dispersion values more than two median
#' absolute deviations away from the trend line, are colored red. See \code{
#' \link[DESeq2]{estimateDispersions}} for more details. For more comprehensive
#' statistical background, see the original DESeq paper (Anders & Huber, 2010), as
#' well as the DESeq2 paper (Love et al., 2014).
#'
#' \code{plot_dispersion} effectively combines \code{edgeR::\link[edgeR]{plotBCV}}
#' and \code{DESeq2::\link[DESeq2]{plotDispEsts}} into a single function that can
#' take either a \code{DGEList} or a \code{DESeqDataSet} as its argument and return
#' the matching figure. By default, dispersions are plotted under log10 transform.
#' They may also be displayed under square root transform, in which case the y-axis
#' may be interpreted as the biological coefficient of variation (McCarthy et al.,
#' 2012).
#'
#' @references
#' Anders, S. & Huber, W. (2010).
#' \href{https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106}{
#' Differential expression analysis for sequence count data}. \emph{Genome Biology},
#' 11:R106.
#'
#' McCarthy, D.J., Chen, Y., & Smyth, G.K. (2012).
#' \href{http://dx.doi.org/10.1093/nar/gks042}{Differential expression analysis of
#' multifactor RNA-Seq experiments with respect to biological variation}.
#' \emph{Nucleic Acids Res., 40}(10): 4288-4297.
#'
#' Love, M., Huber, W. & Anders, S. (2014).
#' \href{http://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8}{
#' Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2}.
#' \emph{Genome Biology, 15}(12): 550.
#'
#' Robinson, M.D. and Smyth, G.K. (2008).
#' \href{http://biostatistics.oxfordjournals.org/content/9/2/321}{Small-sample
#' estimation of negative binomial dispersion, with applications to SAGE data}.
#' \emph{Biostatistics, 9}(2):321-332.
#'
#' Robinson, M.D. & Smyth, G.K. (2007).
#' \href{http://bioinformatics.oxfordjournals.org/content/23/21/2881}{Moderated
#' statistical tests for assessing differences in tag abundance}.
#' \emph{Bioinformatics, 23}(21): 2881-2887.
#'
#' @examples
#' library(DESeq2)
#' dds <- makeExampleDESeqDataSet()
#' plot_dispersion(dds)
#'
#' @seealso
#' \code{\link[DESeq2]{plotDispEsts}, \link[edgeR]{plotBCV},
#' \link[DESeq2]{estimateDispersions}, \link[edgeR]{estimateDisp}}
#'
#' @export
#' @importFrom DESeq2 sizeFactors normalizationFactors estimateSizeFactors
#'   estimateDispersions counts dispersions
#' @importFrom edgeR calcNormFactors estimateDisp aveLogCPM
#' @import dplyr
#' @import ggplot2
#' @importFrom ggsci pal_d3
#' @importFrom plotly ggplotly
#'

plot_dispersion <- function(dat,
                            design = NULL,
                             trans = 'log',
                             title = NULL,
                            legend = 'outside',
                             hover = FALSE) {

  # Preliminaries
  if (is(dat, 'DESeqDataSet')) {
    if (is.null(sizeFactors(dat)) & is.null(normalizationFactors(dat))) {
      dat <- estimateSizeFactors(dat)
    }
    if (!is.null(design)) {
      dat <- estimateDispersions(dat, modelMatrix = design, quiet = TRUE)
    } else if (is.null(dispersions(dat))) {
      dat <- estimateDispersions(dat, quiet = TRUE)
    }
  } else if (is(dat, 'DGEList')) {
    if (all(dat$samples$norm.factors == 1L)) {
      dat <- calcNormFactors(dat)
    }
    if (any(c(is.null(dat$common.dispersion), is.null(dat$trended.dispersion),
              is.null(dat$tagwise.dispersion)))) {
      if (is.null(design) & !is.null(dat$design)) {
        design <- dat$design
      }
      dat <- estimateDisp(dat, design = design)
    }
  } else {
    stop('dat must be a DESeqDataSet or DGEList object.')
  }
  if (trans == 'log') {
    fn <- log10
    ylab <- expression(log[10]~'Dispersion')
  } else if (trans == 'sqrt') {
    fn <- sqrt
    ylab <- 'Biological Coefficient of Variation'
  } else {
    stop('trans must be either "log" or "sqrt".')
  }
  if (is.null(title)) {
    title <- 'Mean-Dispersion Plot'
  }
  if (!legend %in% c('outside', 'bottomleft', 'bottomright', 'topleft', 'topright')) {
    stop('legend must be one of "outside", "bottomleft", "bottomright", ',
         '"topleft", or "topright".')
  }

  # Tidy data
  if (is(dat, 'DESeqDataSet')) {
    keep <- mcols(dat)$baseMean > 0L
    dat <- dat[keep, , drop = FALSE]
    cnts <- DGEList(counts(dat))
    if (is.null(sizeFactors(dat))) {             # Need size factors for aveLogCPM calculation
      cnts <- calcNormFactors(cnts, method = 'RLE')
    } else {
      cnts$samples$norm.factors <- sizeFactors(dat)
    }
    df <- data_frame(Gene = rownames(dat),
                     Mean = aveLogCPM(cnts, prior.count = 1L,
                                      dispersion = dispersions(dat)),
                 Genewise = fn(mcols(dat)$dispGeneEst),
                      Fit = fn(mcols(dat)$dispFit),
                    Final = fn(dispersions(dat)),
                  Outlier = mcols(dat)$dispOutlier)
  } else {
    cmn <- fn(dat$common.dispersion)
    df <- data_frame(Gene = rownames(dat),
                     Mean = aveLogCPM(dat, prior.count = 1L),
                 Genewise = fn(dat$tagwise.dispersion),
                      Fit = fn(dat$trended.dispersion),
                  Tagwise = rep(TRUE, nrow(dat))) %>%
      arrange(Mean)
  }

  # Build plot
  size <- probe_ptsize(df)
  alpha <- probe_alpha(df)
  p <- ggplot(df) +
    labs(title = title, x = expression('Mean'~log[2]*'-CPM'), y = ylab) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (is(dat, 'DESeqDataSet')) {
    p <- p + suppressWarnings(
        geom_point(aes(Mean, Genewise, text = Gene, color = Outlier),
                   size = size, alpha = alpha)
      ) + geom_point(data = df %>% filter(Outlier == FALSE),
                     aes(Mean, Final, color = 'Final'), size = size, alpha = alpha) +
      geom_point(aes(Mean, Fit, color = 'Trend'), size = size) +
      scale_color_manual(name = 'Dispersion\n Estimate',
                       breaks = c(FALSE, 'Trend', 'Final', TRUE),
                       labels = c('Genewise', 'Trend', 'Final', 'Outlier'),
                       values = c('black', pal_d3()(4)[c(1:2, 4)]),
                        guide = guide_legend(override.aes = list(size = rep(1L, 4L))))
  } else {
    p <- p + suppressWarnings(
        geom_point(aes(Mean, Genewise, text = Gene, color = Tagwise),
                   size = size, alpha = alpha)
      ) + geom_hline(aes(color = 'Common', yintercept = cmn), size = 0.5) +
      geom_path(aes(Mean, Fit, color = 'Trend'), size = 0.5) +
      scale_color_manual(name = 'Dispersion\n Estimate',
                       breaks = c(TRUE, 'Common', 'Trend'),
                       labels = c('Genewise', 'Common', 'Trend'),
                       values = c(pal_d3()(2), 'black'),
                        guide = guide_legend(override.aes = list(
                          linetype = c('blank', rep('solid', 2L)),
                             shape = c(16L, NA, NA), size = rep(1L, 3L)
                        )))
  }
  p <- locate_legend(p, legend)

  # Output
  gg_out(p, hover, legend)

}

