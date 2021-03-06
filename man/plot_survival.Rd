% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_survival.R
\name{plot_survival}
\alias{plot_survival}
\title{Survival Curve(s)}
\usage{
plot_survival(
  fit,
  fun = NULL,
  CI = FALSE,
  censor = TRUE,
  pal_curves = "npg",
  title = ifelse(is.null(fit$strata), "Survival Curve", "Survival Curves"),
  leg.txt = NULL,
  legend = "right",
  ...
)
}
\arguments{
\item{fit}{An object of class \code{\link[survival]{survfit}}.}

\item{fun}{An arbitrary function defining a transformation of the survival
curve(s). Common transformations can be specified with a character
argument, e.g. \code{"event"} for cumulative events, \code{"cumhaz"} for
the cumulative hazard function, or \code{"pct"} for survival probability as
a percentage.}

\item{CI}{Plot confidence intervals?}

\item{censor}{Include tick-marks to indicate censored subjects?}

\item{pal_curves}{String specifying the color palette to use when plotting
multiple vectors. Options include \code{"ggplot"}, all qualitative color 
schemes available in \code{\href{https://bit.ly/2ipuEjn}{RColorBrewer}}, 
and the complete collection of \code{\href{http://bit.ly/2bxnuGB}{ggsci}} 
palettes. Alternatively, a character vector of colors with length equal to 
the number of strata in \code{fit}.}

\item{title}{Optional plot title.}

\item{leg.txt}{Optional legend title.}

\item{legend}{Legend position. Must be one of \code{"bottom"}, \code{"left"},
\code{"top"}, \code{"right"}, \code{"bottomright"}, \code{"bottomleft"},
\code{"topleft"}, or \code{"topright"}.}

\item{...}{Additional arguments to be passed to \code{\link[survminer]{
ggsurvplot}}.}
}
\description{
This function plots one or several survival curves for a given time-to-event
model.
}
\details{
Survival curves visualize the fit of a time-to-event model, e.g. a
Kaplan-Meier estimator or a Cox proportional hazards regression. They are an
essential tool in survival analysis, providing a simple and intuitive visual
summary of the relative risk associated with different conditions.

\code{plot_survival} is a lightweight wrapper for the \code{
\link[survminer]{ggsurvplot}} function from the \code{survminer} package,
which offers an impressive array of customization options for plotting
survival curves. \code{plot_survival} alters the default output of \code{
ggsurvplot} to align it with that of other figures generated by \code{
bioplotr}, but you can override those settings and/or take advantage of extra
\code{ggsurvplot} functionalities by passing additional arguments to \code{
plot_survival}.
}
\examples{
library(survival)
fit <- survfit(Surv(time, status) ~ sex, data = lung)
plot_survival(fit)

}
\references{
Andersen, P. & Gill, R. (1982).
\href{https://www.jstor.org/stable/2240714?seq=1#page_scan_tab_contents}{
Cox's regression model for counting processes, a large sample study}. \emph{
Annals of Statistics}, \strong{10}, 1100-1120.

Cox, D.R. (1972).
\href{https://www.jstor.org/stable/2985181?seq=1#page_scan_tab_contents}{
Regression Models and Life-Tables}. \emph{J. R. Stat. Soc., Series B, 34}(2):
187-220.

Kaplan, E.L. & Meier, P. (1958).
\href{https://www.jstor.org/stable/2281868?seq=1#page_scan_tab_contents}{
Nonparametric estimation from incomplete observations}. \emph{J. Amer. Stat.
Assn., 53}(282): 457-481.
}
