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


