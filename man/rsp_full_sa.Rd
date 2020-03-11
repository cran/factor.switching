\name{rsp_full_sa}
\alias{rsp_full_sa}
\title{
Rotation-Sign-Permutation (RSP) algorithm (Full Simulated Annealing)
}
\description{
Rotation-Sign-Permutation (RSP) algorithm (Full Simulated Annealing).
}
\usage{
rsp_full_sa(lambda_mcmc, maxIter = 1000, threshold = 1e-06, verbose = TRUE, 
	sa_loops, rotate = TRUE, increaseIter = FALSE, 
	temperatureSchedule = NULL, printIter = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda_mcmc}{
Input matrix containing a MCMC sample of factor loadings. The column names should read as 'LambdaV1_1',..., 'LambdaV1_q', ..., 'LambdaVp_1',..., 'LambdaVp_q', where \eqn{p} and \eqn{q} correspond to the number of variables and factors, respectively. 
}
  \item{maxIter}{
Maximum number of iterations of the RSP algorithm. Default: 1000.
}
  \item{threshold}{
Positive threshold for declaring convergence. The actual convergence criterion is \code{threshold} m p q with m denoting the number of MCMC iterations. Default: 1e-6.
}
  \item{verbose}{
Logical value indicating whether to print intermediate output or not.
}
  \item{sa_loops}{
Number of simulated annealing loops per MCMC draw.
}
  \item{rotate}{
Logical. Default: TRUE.
}
  \item{increaseIter}{
Logical.
}
  \item{temperatureSchedule}{
Single valued function describing the temperature cooling schedule for the simulated annealing loops.
}
  \item{printIter}{
Print the progress of the algorithm when processing \code{printIter} MCMCdraws, per iteration. Default: 1000.
}
}
\details{
If necessary, more details than the description above.
}
\value{
\item{lambda_reordered_mcmc }{Post-processed MCMC sample of factor loadings.}
\item{sign_vectors }{The final sign-vectors.}
\item{permute_vectors }{The final permutations.}
\item{lambda_hat }{The resulting average of the post-processed MCMC sample of factor loadings.}
\item{objective_function }{A two-column matrix containing the time-to-reach and the value of the objective function for each iteration.}
}
\references{
put references to the literature/web site here.
}
\author{
Panagiotis Papastamoulis
}
\note{
further notes.
}

\examples{
# load small mcmc sample of 100 iterations
#	with p=6 variables and q=2 factors.
data(small_posterior_2chains)
# post-process it
reorderedPosterior <- rsp_partial_sa(
	lambda_mcmc = small_posterior_2chains[[1]], sa_loops=5)
#	sa_loops should be larger in general
# summarize the post-processed MCMC sample with coda
summary(reorderedPosterior$lambda_reordered_mcmc)

}
