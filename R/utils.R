#' Utilities
#'
#' These utility functions assign reasonable size and alpha parameters based on
#' data dimensionality.
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

probe_alpha <- function(df) {
  if (nrow(df) <= 1e5L) {
    out <- 0.5
  } else {
    out <- 0.25
  }
  return(out)
}

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

sample_alpha <- function(df) {
  if (nrow(df) <= 20L) {
    out <- 1L
  } else {
    out <- 0.85
  }
  return(out)
}


