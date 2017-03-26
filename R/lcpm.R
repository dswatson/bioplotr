#' Transform Counts
#'
#' This convenient wrapper function converts raw counts to the log2-counts per million
#' scale, following normalization for library size and a minimal count shift.
#'
#' @param mat Probe by sample matrix of raw counts.
#' @param filter Optional vector of length 2 specifying the filter criterion. Each
#'   probe must have at least \code{filter[1]} log2-counts per million in at least
#'   \code{filter[2]} libraries to pass the expression threshold.
#' @param method Normalization method to be used. See Details.
#'
#' @details
#' \code{lcpm} applies the \code{\link[limma]{voom}} transformation to sequencing
#' data, converting counts to approximately normal distributions on the log2-CPM
#' scale. Data can now be modelled using traditional linear techniques (at least
#' once spot weights have been applied; see Law, et al. (2014)), or used for
#' unsupervised clustering analysis via PCA, MDS, or other methods.
#'
#' It is recommended that low count genes be filtered out prior to transformation.
#' There is no general algorithm for determining the most appropriate expression
#' filter for a given data set. As a rule of thumb, the \code{limma} authors advise
#' setting \code{filter[1]} to either 1, or 10 / (\emph{L} / 1,000,000), where
#' \emph{L} = the minimum library size for a given count matrix. The former
#' corresponds to a log2-CPM of 0, while the latter may be preferable for studies with
#' especially deep read coverage. For \code{filter[2]}, the authors recommend using
#' the number of replicates in the largest group, to guarantee that a gene is expresed
#' in at least one sample for any groupwise comparison. These are broad guidelines,
#' however, not strict rules. See Law, et al. (2016) for some worked through examples.
#'
#' \code{method = "TMM"} is the weighted trimmed mean of M-values (to the reference)
#' proposed by Robinson & Oshlack (2010), where the weights are from the delta method
#' on binomial data.
#'
#' \code{method = "RLE"} is the scaling factor method proposed by Anders & Huber
#' (2010). We call it "relative log expression", as median library is calculated from
#' the geometric mean of all columns and the median ratio of each sample to the median
#' library is taken as the scale factor.
#'
#' \code{method = "upperquartile"} is the upper-quartile normalization method of
#' Bullard, et al. (2010), in which the scale factors are calculated from the 75%
#' quantile of the counts for each library, after removing genes which are zero in all
#' libraries.
#'
#' If \code{method = "none"}, then the normalization factors are set to 1.
#'
#' For symmetry, normalization factors are adjusted to multiply to 1. The effective
#' library size is then the original library size multiplied by the scaling factor.
#'
#' Note that rows that have zero counts for all columns are trimmed before
#' normalization factors are computed. Therefore rows with all zero counts do not
#' affect the estimated factors.
#'
#' @return A numeric matrix of normalized counts on the log2-CPM scale.
#'
#' @references
#' Anders, S. & Huber, W. (2010). "Differential expression analysis for sequence
#' count data." \emph{Genome Biology}, 11:R106.
#' \url{https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106}
#'
#' Bullard, J.H., Purdom, E., Hansen, K.D. & Dudoit, S. (2010). "Evaluation of
#' statistical methods for normalization and differential expression in mRNA-Seq
#' experiments." \emph{BMC Bioinformatics}, 11:94.
#' \url{http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-94}
#'
#' Law, C.W., Alhamdoosh, M., Su, S., Smyth, G.K., & Ritchie, M.E. (2016). "RNA-seq
#' analysis is easy as 1-2-3 with limma, Glimma and edgeR." \emph{F1000 Research},
#' \emph{5}(1408).
#' \url{https://f1000research.com/articles/5-1408/v2}
#'
#' Law, C.W., Chen, Y., Shi, W., & Smyth, G.K. (2014). "voom: precision weights unlock
#' linear model analysis tools for RNA-seq read counts." \emph{Genome Biology},
#' \strong{15}:R29.
#' \url{https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29}
#'
#' Robinson, M.D. & Oshlack, A. (2010). "A scaling normalization method for
#' differential expression analysis of RNA-seq data." \emph{Genome Biology}, 11:R25.
#' \url{https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25}
#'
#' @examples
#' # Simulate count data
#' mat <- matrix(rnbinom(5000, mu = 4, size = 1), nrow = 1000, ncol = 5)
#'
#' # Plot raw counts
#' plot_density(mat)
#'
#' # Plot transformed counts
#' trans_mat <- lcpm(mat)
#' plot_density(trans_mat)
#'
#' @seealso
#' \code{\link[limma]{voom}, \link[edgeR]{cpm}}
#'
#' @export
#' @importFrom edgeR cpm DGEList calcNormFactors
#'

lcpm <- function(mat,
                 filter = NULL,
                 method = 'TMM') {

  # Preliminaries
  if (!is.null(filter)) {
    if (length(filter) != 2L) {
      stop('filter must be a vector of length 2.')
    }
  }

  # Filter
  if (is.null(filter)) y <- DGEList(mat)
  else {
    keep <- rowSums(cpm(mat) > filter[1]) >= filter[2]
    y <- DGEList(mat[keep, , drop = FALSE])
  }

  # Transform
  y <- calcNormFactors(y, method = method)
  lib.size <- with(y$samples, lib.size * norm.factors)
  y <- t(log2(t(y$counts + 0.5) / (lib.size + 1L) * 1e6L))
  colnames(y) <- colnames(mat)

  # Export
  return(y)

}
