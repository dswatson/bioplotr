#' Size Probewise Data Points
#'
#' This utility function assigns a reasonable point size parameter based on data
#' dimensionality.
#'
#' @param df Data frame to be passed to \code{ggplot}.
#'

probe_ptsize <- function(df) {
  if (nrow(df) <= 1e4L) {                        # miRNA, peptides
    out <- 0.75
  } else if (nrow(df) <= 1e5L) {                 # mRNA, genes, metabolites
    out <- 0.5
  } else {                                       # transcripts, exons, methylation loci, etc.
    out <- 0.25
  }
  return(out)
}

#' Set Transparency for Probewise Data Points
#'
#' This utility function assigns a reasonable alpha parameter based on data
#' dimensionality.
#'
#' @param df Data frame to be passed to \code{ggplot}.
#'

probe_alpha <- function(df) {
  if (nrow(df) <= 1e5L) {
    out <- 0.5
  } else {
    out <- 0.25
  }
  return(out)
}

#' Size Samplewise Data Points
#'
#' This utility function assigns a reasonable point size parameter based on data
#' dimensionality.
#'
#' @param df Data frame to be passed to \code{ggplot}.
#'

sample_ptsize <- function(df) {
  if (nrow(df) <= 10L) {
    out <- 3L
  } else if (nrow(df) <= 20L) {
    out <- 2L
  } else {
    out <- 1.5
  }
  return(out)
}

#' Set Transparency for Samplewise Data Points
#'
#' This utility function assigns a reasonable alpha parameter based on data
#' dimensionality.
#'
#' @param df Data frame to be passed to \code{ggplot}.
#'

sample_alpha <- function(df) {
  if (nrow(df) <= 20L) {
    out <- 1L
  } else {
    out <- 0.85
  }
  return(out)
}

#' Output Image
#'
#' This utility function locates the figure legend and prints a ggplot or
#' ggplotly figure.
#'
#' @param p A \code{ggplot2} object.
#' @param hover Add text tooltip using plotly?
#' @param loc String specifying legend location.
#'
#' @importFrom ggplot2 theme
#' @importFrom plotly ggplotly
#'

gg_out <- function(p,
                   hover,
                   loc = NULL) {

  # Locate legend
  if (loc == 'bottomleft') {
    p <- p + theme(legend.justification = c(0.01, 0.01),
                   legend.position = c(0.01, 0.01))
  } else if (loc == 'bottomright') {
    p <- p + theme(legend.justification = c(0.99, 0.01),
                   legend.position = c(0.99, 0.01))
  } else if (loc == 'topleft') {
    p <- p + theme(legend.justification = c(0.01, 0.99),
                   legend.position = c(0.01, 0.99))
  } else if (loc == 'topright') {
    p <- p + theme(legend.justification = c(0.99, 0.99),
                   legend.position = c(0.99, 0.99))
  }

  # Print figure
  if (!hover) {
    print(p)
  } else {
    if (loc == 'outside') {
      p <- ggplotly(p, tooltip = 'text', height = 600, width = 650)
    } else {
      p <- ggplotly(p, tooltip = 'text', height = 600, width = 600)
    }
    print(p)
  }

}

#' Standardize Matrix
#'
#' This utility function takes data objects from \code{limma}, \code{edgeR}, or
#' \code{DESeq2} pipelines and outputs a standard probe by sample matrix. Raw counts
#' are log2-CPM transformed with a warning.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples.
#'
#' @importFrom edgeR cpm calcNormFactors
#' @importFrom DESeq2 sizeFactors normalizationFactors estimateSizeFactors counts
#' @importFrom SummarizedExperiment assay
#' @importFrom limma getEAWP
#'

matrixize <- function(dat) {


  if (is(dat, 'DGEList')) {
    keep <- rowSums(dat$counts) > 1L             # Minimal count filter
    dat <- dat[keep, ]
    if (is.null(dat$samples$norm.factors) |      # Calculate size factors
        all(dat$samples$norm.factors == 1L)) {
      dat <- calcNormFactors(dat)
    }
    dat <- cpm(dat, log = TRUE, prior.count = 1L)
    warning('Transforming raw counts to log2-CPM scale.')
  } else if (is(dat, 'DESeqDataSet')) {
    if (is.null(sizeFactors(dat)) & is.null(normalizationFactors(dat))) {
      dat <- estimateSizeFactors(dat)            # Normalize counts
    }
    dat <- counts(dat, normalized = TRUE)
    keep <- rowMeans(dat) > 0L                   # Minimal count filter
    dat <- dat[keep, , drop = FALSE]
    dat <- cpm(dat, log = TRUE, prior.count = 1L)
    warning('Transforming raw counts to log2-CPM scale.')
  } else if (is(dat, 'DESeqTransform')) {
    dat <- assay(dat)
  } else {
    dat <- getEAWP(dat)$expr
    keep <- rowSums(is.finite(dat)) == ncol(dat)
    dat <- dat[keep, , drop = FALSE]
  }

  # Export
  return(dat)

}

#' Create Distance Matrix
#'
#' This utility function calculates distance based on a user defined measure
#' and optionally filters probes by leading log fold change.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples.
#' @param top Optional number (if > 1) or proportion (if < 1) of top probes to be used
#'   for distance calculations.
#' @param dist Distance measure to be used. Currently supports \code{"euclidean",
#'   "pearson", "MI",} or \code{"KLD"}.
#'
#' @importFrom wordspace dist.matrix
#' @importFrom bioDist MIdist KLdist.matrix
#' @importFrom KernSmooth dpih
#'

dist_mat <- function(dat,
                     top,
                     dist) {

  # Preliminaries
  if (!is.null(top)) {
    if (top > 1L) {
      if (top > nrow(dat)) {
        warning(paste('top exceeds nrow(dat), at least after removing probes with
                      missing values and/or applying a minimal expression filter.
                      Proceeding with the complete', nrow(dat), 'x', ncol(dat), 'matrix.'))
        top <- NULL
      }
    } else {
      top <- round(top * nrow(dat))
    }
    }

  # Create matrix
  if (is.null(top)) {
    if (dist == 'euclidean') {
      dm <- dist.matrix(t(dat), method = 'euclidean')
    } else if (dist == 'pearson') {
      dm <- 1 - cor(dat)
    } else if (dist == 'MI') {
      dm <- as.matrix(MIdist(t(dat)))
    } else if (dist == 'KLD') {
      dm <- as.matrix(KLdist.matrix(t(dat)))
    }
  } else {
    dm <- matrix(0L, nrow = ncol(dat), ncol = ncol(dat))
    for (i in 2L:ncol(dat)) {
      for (j in 1L:(i - 1L)) {
        if (dist == 'euclidean') {
          top_idx <- nrow(dat) - top + 1L
          dm[i, j] <- sqrt(sum(sort.int((dat[, i] - dat[, j])^2L,
                                        partial = top_idx)[top_idx:nrow(dat)]))
        } else {
          tops <- order((dat[, i] - dat[, j])^2, decreasing = TRUE)[seq_len(top)]
          if (dist == 'pearson') {
            dm[i, j] <- 1 - cor(dat[tops, i], dat[tops, j])
          } else if (dist == 'MI') {
            dm[i, j] <- max(as.matrix(MIdist(t(dat[tops, c(i, j)]))))
          } else if (dist == 'KLD') {
            dm[i, j] <- max(as.matrix(KLdist.matrix(t(dat[tops, c(i, j)]))))
          }
        }
      }
    }
    dm <- pmax(dm, t(dm))
  }

  # Export
  return(dm)

}


