\name{weighted_procrustes_switching}
\alias{weighted_procrustes_switching}
\title{
Weighted Orthogonal Procrustes rotations
}
\description{
Weighted Orthogonal Procrustes (WOP) post-processing (Assmann et al. 2016) augmented with a final varimax rotation as implemented in Papastamoulis and Ntzoufras (2020). The algorithm uses the \code{procrustes} function of the MCMCpack package. 
}
\usage{
weighted_procrustes_switching(lambda_mcmc, maxIter, threshold, verbose, 
    weight, printIter)
}
\arguments{
  \item{lambda_mcmc}{
Input matrix containing a MCMC sample of factor loadings. The column names should read as 'LambdaV1_1',..., 'LambdaV1_q', ..., 'LambdaVp_1',..., 'LambdaVp_q', where \eqn{p} and \eqn{q} correspond to the number of variables and factors, respectively. 
}
  \item{maxIter}{
Maximum number of iterations of the RSP algorithm. Default: 100.
}
  \item{threshold}{
Positive threshold for declaring convergence. The actual convergence criterion is \code{threshold} m p q with m denoting the number of MCMC iterations. Default: 1e-6.
}
  \item{verbose}{
Logical value indicating whether to print intermediate output or not. 
}
  \item{weight}{
This is argument is always set to TRUE. 
}
  \item{printIter}{
Print the progress of the algorithm when processing \code{printIter} MCMCdraws, per iteration. Default: 1000.
}
}
\value{
\item{lambda_reordered_mcmc }{Post-processed MCMC sample of factor loadings.}
\item{lambda_hat }{The resulting average of the post-processed MCMC sample of factor loadings.}
\item{objective_function }{A two-column matrix containing the time-to-reach and the value of the objective function for each iteration.}

}
\references{
  Assmann, C., Boysen-Hogrefem J. and Pape M. (2016). 
  Bayesian analysis of static and dynamic factor models: An ex-post approach towards the rotation problem.
  \emph{Journal of Econometrics}: 192 (1): Pages 190-206. 


   Martin AD, Quinn KM, Park JH (2011). 
   MCMCpack: Markov Chain Monte Carlo in R. 
   \emph{Journal of Statistical Software}: 42(9), 22.


  Papastamoulis, P. and Ntzoufras, I. (2020).
  On the identifiability of Bayesian Factor Analytic models.
  \emph{arXiv:2004.05105 [stat.ME]}.
  
}
\author{
Panagiotis Papastamoulis
}

\examples{
# load small mcmc sample of 100 iterations
#	with p=6 variables and q=2 factors.
data(small_posterior_2chains)
# post-process it
reorderedPosterior <- weighted_procrustes_switching(
	lambda_mcmc = small_posterior_2chains[[1]])
# summarize the post-processed MCMC sample with coda
summary(reorderedPosterior$lambda_reordered_mcmc)
}
