\name{factor.switching-package}
\alias{factor.switching-package}
\alias{factor.switching}
\docType{package}
\title{
\packageTitle{factor.switching}
}
\description{
\packageDescription{factor.switching}

There are three alternative schemes for minimizing the objective function. 
\enumerate{
   \item{Exact \code{\link{rsp_exact}}}
   \item{Partial Simulated Annealing \code{\link{rsp_partial_sa}}}
   \item{Full simulated annealing \code{\link{rsp_full_sa}}}
}

The exact algorithm solves \eqn{2^q} assignment problems per MCMC iteration, where \eqn{q} denotes the number of factors of the fitted model. For typical values of  the number of factors (e.g. \eqn{q<11}) the exact scheme should be preferred. Otherwise, the two approximate algorithms based on simulated annealing may be considered. The Partial simulated annealing is more efficient than the full simulated annealing scheme. 

In cases of parallel MCMC chains, applying the RSP algorithm for each chain separately will  identify the factor loadings within each chain. However, the results will not be comparable between chains. The comparison of multiple MCMC chains is doable via the \code{\link{compareMultipleChains}} function. 
}
\details{

The DESCRIPTION file:
\packageIndices{factor.switching}
}
\author{
Panagiotis Papastamoulis

Maintainer: Panagiotis Papastamoulis
}
\references{
  Papastamoulis, P. and Ntzoufras, I. (2020).
  On the identifiability of Bayesian Factor Analytic models.
  \emph{arXiv:2004.05105 [stat.ME]}.
}
\keyword{ package }
\seealso{
\code{\link{rsp_exact}}, \code{\link{plot.rsp}}, \code{\link{compareMultipleChains}}
}
\examples{
	# load 2 chains each one consisting of a 
	#	small mcmc sample of 100 iterations
	#	with p=6 variables and q=2 factors.
	data(small_posterior_2chains)
	Nchains <- length(small_posterior_2chains)
	reorderedPosterior <- vector('list',length=Nchains)
	# post-process the 2 chains
	for(i in 1:Nchains){
		reorderedPosterior[[i]] <- rsp_exact( lambda_mcmc = small_posterior_2chains[[i]], 
				maxIter = 100, 
				threshold = 1e-6, 
				verbose=TRUE )
	}
	#	plot posterior summary for chain 1:
	plot(reorderedPosterior[[1]])
	#	plot posterior summary for chain 2:
	plot(reorderedPosterior[[2]])
	# make them comparable
	makeThemSimilar <- compareMultipleChains(rspObjectList=reorderedPosterior)
	# plot the traces of both chains
	oldpar <- par(no.readonly =TRUE)
	par(mfcol=c(2,6),mar=c(4,4,2,1))
	plot(makeThemSimilar,auto.layout=FALSE,density=FALSE,
		ylim=c(-1.1,1.1),smooth=FALSE,col=c('red','blue'))
	legend('topright',c('post-processed chain 1', 
		'post-processed chain 2'),lty=1:2,col=c('red','blue'))
	par(oldpar)
	# you can also use the summary of mcmc.list
	summary(makeThemSimilar)
}
