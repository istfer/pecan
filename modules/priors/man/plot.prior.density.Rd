% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot.prior.density}
\alias{plot.prior.density}
\title{Add Prior Density}
\usage{
\method{plot}{prior.density}(prior.density, base.plot = NULL,
  prior.color = "black")
}
\arguments{
\item{prior.density}{}

\item{base.plot}{a ggplot object (grob), created by \code{\link{create.base.plot}} if none provided}

\item{prior.color}{color of line to be plotted}
}
\value{
plot with prior density added
}
\description{
Plots a prior density from a parameterized probability distribution
}
\examples{
\dontrun{
plot.prior.density(pr.dens('norm', 0, 1))
}
}
\seealso{
\code{\link{pr.dens}}
}
\author{
David LeBauer
}
