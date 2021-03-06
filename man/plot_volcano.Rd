% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_volcano.R
\name{plot_volcano}
\alias{plot_volcano}
\title{Volcano Plot}
\usage{
plot_volcano(
  dat,
  y = "p",
  fdr = 0.05,
  lfc = NULL,
  probes = NULL,
  size = NULL,
  alpha = NULL,
  title = "Volcano Plot",
  legend = "right",
  hover = FALSE
)
}
\arguments{
\item{dat}{Data frame or similar object representing the results of a test
for differential expression, such as the output of a call to \code{
limma::\link[limma]{topTable}}, \code{edgeR::\link[edgeR]{topTags}}, or
\code{DESeq2::\link[DESeq2]{results}}. Alternatively, any object coercable
to a data frame with columns for \emph{p}-values, \emph{q}-values, and log
fold changes. Missing values are silently removed.}

\item{y}{Plot \emph{p}-values (\code{y = "p"}), \emph{q}-values 
(\code{y = "q"}), or log-odds (\code{y = "B"}) on the y-axis?}

\item{fdr}{Significance threshold for declaring a probe differentially
expressed.}

\item{lfc}{Optional effect size threshold for declaring a probe
differentially expressed.}

\item{probes}{Optional column number or name specifying where probe names are
stored, presuming they are not stored in \code{rownames(dat)}.}

\item{size}{Point size.}

\item{alpha}{Point transparency.}

\item{title}{Optional plot title.}

\item{legend}{Legend position. Must be one of \code{"bottom"}, \code{"left"},
\code{"top"}, \code{"right"}, \code{"bottomright"}, \code{"bottomleft"},
\code{"topleft"}, or \code{"topright"}.}

\item{hover}{Show probe name by hovering mouse over data point? If \code{
TRUE}, the plot is rendered in HTML and will either open in your browser's
graphic display or appear in the RStudio viewer. Probe names are extracted
from \code{dat}.}
}
\description{
This function plots effect size against significance for a given test of
differential expression.
}
\details{
Volcano plots visualize the relationship between each probe's effect size 
and significance for a given test of differential expression. Points are 
colored to distinguish between those that do and do not meet a user-defined 
FDR threshold. Up- and down-regulated genes may also be differentially 
colored if a minimum absolute fold change is supplied. These figures help to 
evaluate the symmetry, magnitude, and significance of effects in an omic 
experiment.
}
\examples{
# Simulated data
library(DESeq2)
dds <- makeExampleDESeqDataSet()
dds <- DESeq(dds)
res <- results(dds)
plot_volcano(res)

# Real data
data(airway)
dds <- DESeqDataSet(airway, design = ~ cell + dex)
dds <- DESeq(dds)
res <- results(dds)
plot_volcano(res)

}
