\name{plot.rsp}
\alias{plot.rsp}
\title{
Plot posterior means and credible regions
}
\description{
This function plot posterior mean estimates per factor along with Highest Density Intervals, as well as simultaneous credible regions.
}
\usage{
\method{plot}{rsp}(x, prob, myCol, mfrow, subSet, simCR, HDI, ...)
}
\arguments{
  \item{x}{
An object of class \code{rsp}.
}
  \item{prob}{
Coverage probability of credible regions.
}
  \item{myCol}{
Vector of colours.
}
  \item{mfrow}{
Number of rows and columns in the resulting graphic.
}
  \item{subSet}{
Enable to plot a subset of factors.
}
  \item{simCR}{
Logical value for plotting simultaneous credible regions. Default: True.
}
  \item{HDI}{
Logical value for plotting Highest Density Intervals per factor loading. Default: True.
}
  \item{...}{
Ignored
}
}
\value{
A plot.
}
\author{
Panagiotis Papastasmoulis
}
\examples{
# load small mcmc sample of 100 iterations
#	with p=6 variables and q=2 factors.
data(small_posterior_2chains)
# post-process it
reorderedPosterior <- rsp_exact(
	lambda_mcmc = small_posterior_2chains[[1]])
# plot it
plot(reorderedPosterior, mfrow = c(1,2), prob=0.95)
}
