% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_taco.R
\name{plot_taco}
\alias{plot_taco}
\title{Taco Plot}
\usage{
plot_taco(dat, fdr = 0.05, title = "Taco Plot", legend = "right")
}
\arguments{
\item{dat}{Data frame representing the results of a test for differential
expression, such as the output of a call to \code{
limma::\link[limma]{topTable}}, \code{edgeR::\link[edgeR]{topTags}}, or
\code{DESeq2::\link[DESeq2]{results}}. Alternatively, any object coercable
to a data frame with columns for log fold changes, probewise means,
\emph{p}-values, and FDR. Missing values are silently removed.}

\item{fdr}{Significance threshold for declaring a probe differentially
expressed.}

\item{title}{Optional plot title.}

\item{legend}{Legend position. Must be one of \code{"bottom"}, \code{"left"},
\code{"top"}, \code{"right"}, \code{"bottomright"}, \code{"bottomleft"},
\code{"topleft"}, or \code{"topright"}.}
}
\description{
This function plots effect size against significance against probewise means
for a given test of differential expression.
}
\details{
A taco plot combines the elements of a volcano plot and an MD plot into a
single three-dimensional figure. Points are colored to distinguish between
those that do and do not meet a user-defined FDR threshold. These
figures help to evaluate the symmetry, magnitude, and significance of effects
in an omic experiment.
}
\examples{
library(limma)
DE_genes <- cbind(matrix(rnorm(50 * 5, mean = 5), nrow = 50, ncol = 5),
                  matrix(rnorm(50 * 5), nrow = 50, ncol = 5))
eset <- rbind(DE_genes, matrix(rnorm(4950 * 10), nrow = 4950, ncol = 10))
treat <- rep(c("A", "B"), each = 5)
des <- model.matrix(~ treat)
fit <- eBayes(lmFit(eset, des))
top <- topTable(fit, number = Inf)
plot_taco(top)

}
