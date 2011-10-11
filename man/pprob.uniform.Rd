\name{pprob.uniform}
\alias{pprob.uniform}
\title{Bayesian diagnostic test for multiple testing p-values.}
\description{
This function accepts a vector of simulated null p-values from a single simulated study. 
The null p-values should representa subset of all the simulated p-values corresponding 
to the tests with no signal. 
}
\usage{
  pprob.uniform(p,alpha=c(0.1,10),beta=c(0.1,10),eps=1e-10)
}
\arguments{
  \item{p}{An vector of null p-values from a single simulated study.}
  \item{alpha}{The range of the first parameter for the prior on the beta distribution.}
  \item{beta}{The range of the second parameter for the prior on the beta distribution.}
  \item{eps}{Maximum integration error when computing the posterior distribution.}	
}
\details{
  The pprob.uniform function calculates the posterior probability that a set of null p-values
  come from the uniform distribution as described in Leek and Storey (2009). The p-values
  should be simulated from a realistic distribution and only the null p-values should
  be passed to the pprob.uniform function. 
}

\value{
  \item{pp}{The posterior probability that p is a sample from the uniform distribution.}
}
\references{
  J.T. Leek and J.D. Storey, "The Joint Null Distribution of Multiple Hypothesis Tests."
}

\author{Jeffrey T. Leek \email{jleek@jhsph.edu}}
\seealso{\code{\link{dks}}, \code{\link{dks.pvalue}}, \code{\link{pprob.dist}},\code{\link{cred.set}}}
\examples{
  ## Load data
  data(dksdata) 
  pp <- pprob.uniform(P[,1])
  hist(pp)
}

\keyword{misc}
