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

#' Locate Legend
#'
#' This utility function translates a user-supplied string into coordinates for
#' legend location.
#'
#' @param p A \code{ggplot2} object.
#' @param loc String specifying legend location.
#'
#' @importFrom ggplot2 theme
#'

locate_legend <- function(p, loc) {
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
  return(p)
}

#' Output Image
#'
#' This utility function prints a ggplot or ggplotly figure.
#'
#' @param p A \code{ggplot2} object.
#' @param hover Add text tooltip using plotly?
#' @param loc String specifying legend location.
#'
#' @importFrom plotly ggplotly
#'

gg_out <- function(p, hover, loc) {
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


