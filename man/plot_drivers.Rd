% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_drivers.R
\name{plot_drivers}
\alias{plot_drivers}
\title{Plot Drivers of Omic Variation}
\usage{
plot_drivers(
  dat,
  clin,
  stat = "p",
  bivariate = TRUE,
  block = NULL,
  unblock = NULL,
  parametric = TRUE,
  kernel = NULL,
  kpar = NULL,
  top = NULL,
  n_pc = 5L,
  alpha = NULL,
  p_adj = NULL,
  r_adj = FALSE,
  label = FALSE,
  pal_tiles = "PiRdBr",
  lim = NULL,
  coord_equal = FALSE,
  title = "Variation By Feature",
  legend = "right",
  hover = FALSE
)
}
\arguments{
\item{dat}{Omic data matrix or matrix-like object with rows corresponding to
probes and columns to samples. It is strongly recommended that data be
filtered and normalized prior to plotting. Raw counts stored in 
\code{\link[edgeR]{DGEList}} or \code{\link[DESeq2]{DESeqDataSet}} objects 
are automatically extracted and transformed to the log2-CPM scale, with a
warning.}

\item{clin}{Data frame or matrix with rows corresponding to samples and
columns to technical and/or biological features to test for associations
with omic data.}

\item{stat}{Association statistic of choice. Currently supports \code{"p"} 
(-log \emph{p}-values) and \code{"r2"} (R-squared). Interpretations vary 
depending on whether covariates are included. See Details.}

\item{bivariate}{Test associations in isolation, or after adjusting for
all remaining covariates? If \code{FALSE}, then \code{clin} is treated as 
a design matrix against which each PC is sequentially regressed. See 
Details.}

\item{block}{String specifying the name of the column in which to find the
blocking variable, should one be accounted for. See Details.}

\item{unblock}{Column name(s) of one or more features for which the 
\code{block} covariate should not be applied, if one was supplied. See 
Details.}

\item{parametric}{Compute statistics using parametric association tests?
If \code{FALSE}, rank-based alternatives are used instead. Either a single
logical value, in which case it applies to all tests, or a logical vector
of length equal to \code{ncol(clin)}. See Details.}

\item{kernel}{The kernel generating function, if using KPCA. Options include
\code{"rbfdot"}, \code{"polydot"}, \code{"tanhdot"}, \code{"vanilladot"}, 
\code{"laplacedot"}, \code{"besseldot"}, \code{"anovadot"}, and 
\code{"splinedot"}. To run normal PCA, set to \code{NULL}.}

\item{kpar}{A named list of arguments setting parameters for the kernel
function. Only relevant if \code{kernel} is not \code{NULL}.}

\item{top}{Optional number (if > 1) or proportion (if < 1) of most variable
probes to be used for PCA.}

\item{n_pc}{Number of principal components to include in the figure.}

\item{alpha}{Optional significance threshold to impose on associations. 
Those with \emph{p}-values (optionally adjusted) less than or equal to 
\code{alpha} are outlined in black.}

\item{p_adj}{Optional \emph{p}-value adjustment for multiple testing. Options
include \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{
"bonferroni"}, \code{"BH"}, \code{"BY"}, and \code{"fdr"}. See \code{
\link[stats]{p.adjust}}.}

\item{r_adj}{Adjust partial R-squared? Only relevant if \code{stat = "r2"} 
and either \code{bivariate = FALSE} or \code{block} is non-\code{NULL}.}

\item{label}{Print association statistics over tiles?}

\item{pal_tiles}{String specifying the color palette to use for heatmap
tiles. Options include the complete collection of \code{\href{
https://bit.ly/2n7D6tF}{viridis}} palettes, as well as all sequential and
divergent color schemes available in \code{\href{
https://bit.ly/2ipuEjn}{RColorBrewer}}. Alternatively, a character vector 
of at least two colors.}

\item{lim}{Optional vector of length two defining lower and upper bounds for 
the scale range. Default is observed extrema for \code{stat = "p"} and the
unit interval for \code{stat = "r2"}.}

\item{coord_equal}{Plot tiles of equal width and height?}

\item{title}{Optional plot title.}

\item{legend}{Legend position. Must be one of \code{"bottom"}, \code{"left"},
\code{"top"}, \code{"right"}, \code{"bottomright"}, \code{"bottomleft"},
\code{"topleft"}, or \code{"topright"}.}

\item{hover}{Show association statistics by hovering mouse over tiles? If 
\code{TRUE}, the plot is rendered in HTML and will either open in your 
browser's graphic display or appear in the RStudio viewer.}
}
\description{
This function visualizes the strength of associations between the principal
components of an omic data matrix and a set of biological and/or technical
features.
}
\details{
Strength of association may be measured either by --log \emph{p}-values (if
\code{stat = "p"}) or R-squared (if \code{stat = "r2"}). The former may be
adjusted for multiple testing, while the latter can be adjusted for 
covariates.

If \code{bivariate = TRUE}, then association tests are performed between
each PC and each clinical covariate, optionally adjusting for a blocking 
variable (if \code{block} is non-\code{NULL}). If \code{bivariate = FALSE},
then all tests are partial association tests, in the sense that they control 
for all remaining covariates. 

When \code{bivariate = TRUE}, \code{block = NULL}, and \code{parametric = 
TRUE}, significance is computed from Pearson correlation tests (for 
continuous features) or ANOVA \emph{F}-tests (for categorical features). When 
\code{parametric = FALSE}, significance is computed from rank-based 
alternatives, i.e. Spearman correlation tests (for continuous features) or 
Kruskal-Wallis tests (for categorical features). 

When \code{bivariate = FALSE} or \code{block} is non-\code{NULL}, 
significance is computed from partial correlation tests for continuous data 
(Pearson if \code{parametric = TRUE}, Spearman if \code{parametric = FALSE}) 
or repeated measures ANOVA \emph{F}-tests (under rank-transformation if 
\code{parametric = FALSE}). In all cases, the alternative hypothesis assumes
a monotonic relationship between variables.

A blocking variable may be provided if samples violate the assumption of 
independence, e.g. for studies in which subjects are observed at multiple 
time points. If a blocking variable is identified, it will be regressed out 
prior to testing for all variables except those explicitly exempted by the 
\code{unblock} argument. When supplying a blocking variable, be careful to 
consider potential collinearities in the data. For instance, clinical 
features may be invariant with respect to subject, while subject may be 
nested within other variables like batch or treatment group. The \code{block} 
and \code{unblock} arguments are intended to help parse out these 
relationships.

Numeric and categorical features are tested differently. To protect against
potential mistakes (e.g., one-hot encoding a Boolean variable), \code{
plot_drivers} automatically prints a data frame listing the class of each 
feature.

If \code{kernel} is non-\code{NULL}, then KPCA is used instead of PCA. See
\code{\link{plot_kpca}} for more info. Details on kernel functions and their
input parameters can be found in \code{kernlab::\link[kernlab]{dots}}.
}
\examples{
library(SummarizedExperiment)
library(edgeR)
library(dplyr)
data(airway)
cnts <- assay(airway)
keep <- rowSums(cpm(cnts) > 1) >= 4
mat <- cpm(cnts[keep, ], log = TRUE)
clin <- colData(airway) \%>\%
  as_tibble(.) \%>\%
  select(cell, dex)
plot_drivers(mat, clin)

}
\seealso{
\code{\link{plot_pca}}, \code{\link{plot_kpca}}
}
