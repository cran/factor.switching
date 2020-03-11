\name{compareMultipleChains}
\alias{compareMultipleChains}
\title{
Compare multiple chains
}
\description{
Compares multiples chains after each one of them has been post-processed by the RSP algorithm, so that all of them are switched into a similar labelling.
}
\usage{
compareMultipleChains(rspObjectList, scheme, sa_loops, maxIter, threshold)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rspObjectList}{
A list consisting of \code{rsp} objects.
}
  \item{scheme}{
Character argument with possible values: "exact" (default), "partial" or "full".
}
  \item{sa_loops}{
Number of simulated annealing loops (only applicable when "exact" is disabled).
}
  \item{maxIter}{
Max number of iterations.
}
  \item{threshold}{
Threshold for convergence.
}
}
\value{
\code{reorderedChains}: an object of class \code{mcmc.list} containing all simultaneously processed chains.
}
\author{
Panagiotis Papastamoulis
}
\examples{
# load 2 chains each one consisting of a 
#	small mcmc sample of 100 iterations
#	with p=6 variables and q=2 factors.
data(small_posterior_2chains)
Nchains <- length(small_posterior_2chains)
reorderedPosterior <- vector('list',length=Nchains)
for(i in 1:Nchains){
	reorderedPosterior[[i]] <- rsp_exact( lambda_mcmc = small_posterior_2chains[[i]], 
			maxIter = 100, 
			threshold = 1e-6, 
			verbose=TRUE )
}
# make them comparable
makeThemSimilar <- compareMultipleChains(rspObjectList=reorderedPosterior)
}

