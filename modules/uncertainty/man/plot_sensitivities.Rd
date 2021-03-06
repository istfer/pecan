% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_sensitivities}
\alias{plot_sensitivities}
\title{Plot Sensitivities}
\usage{
plot_sensitivities(sensitivity.plot.inputs,
  prior.sensitivity.plot.inputs = NULL, ...)
}
\arguments{
\item{sensitivity.plot.inputs}{inputs}

\item{prior.sensitivity.plot.inputs}{priors}

\item{...}{arguments passed to \code{\link{plot_sensitivity}}}

\item{sensitivity.results}{list containing sa.samples and sa.splines}
}
\value{
list of plots, one per trait
}
\description{
Plot functions and quantiles used in sensitivity analysis
}
\details{
Generates a plot using \code{\link{plot_sensitivity}} for multiple traits.
}
