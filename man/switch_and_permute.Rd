\name{switch_and_permute}
\alias{switch_and_permute}
\title{
Apply sign switchings and column permutations
}
\description{
Help function, not really meant to be used by the average user.
}
\usage{
switch_and_permute(lambda_mcmc, switch_vectors, permute_vectors)
}
\arguments{
  \item{lambda_mcmc}{
MCMC input.
}
  \item{switch_vectors}{
Sign vectors.
}
  \item{permute_vectors}{
Permutation vectors.
}
}
\value{
reordered \code{lambda_mcmc} according to sign and permutations provided.
}
\author{
Panagiotis Papastamoulis
}

