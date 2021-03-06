% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_tree.R
\name{plot_tree}
\alias{plot_tree}
\title{Tree Plots}
\usage{
plot_tree(dat, labels = TRUE, leaf_labels = TRUE, rotate = FALSE, title = NULL)
}
\arguments{
\item{dat}{An object representing a dendrogram (e.g. of class \code{
\link[stats]{hclust}}), regression tree (e.g. of class \code{\link[rpart]{
rpart}}), or other tree-like structure.}

\item{labels}{Plot labels? This refers to items for dendrograms, and break
points for regression or classification trees.}

\item{leaf_labels}{Plot leaf labels? Only relevant for regression and
classification trees.}

\item{rotate}{Rotate plot 90 degrees?}

\item{title}{Optional plot title.}
}
\description{
This function plots tree structures such as dendrograms and regression trees.
}
\details{
This function is essentially a wrapper for tools in the \code{\href{
https://github.com/andrie/ggdendro}{ggdendro}} package, with a few amendments
and extra options. See the package documentation for more info.
}
\examples{
hc <- hclust(dist(USArrests))
plot_tree(hc)

}
