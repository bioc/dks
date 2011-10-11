\name{dks}
\alias{dks}
\title{Frequentist and Bayesian diagnostic tests for multiple testing p-values.}
\description{
This function accepts a matrix of simulated null p-values where each column corresponds
to the p-values from a single simulated study. The null p-values should represent
a subset of all the simulated p-values corresponding to the tests with no signal. 
}
\usage{
  dks(P,alpha=c(0.1,10),beta=c(0.1,10),plot=TRUE,eps=1e-10)
}
\arguments{
  \item{P}{An m0 x B matrix of null p-values, each column corresponds to the p-values from a single simulated study.}
  \item{alpha}{The range of the first parameter for the prior on the beta distribution.}
  \item{beta}{The range of the second parameter for the prior on the beta distribution.}
  \item{plot}{Should diagnostic plots be displayed.}
  \item{eps}{Maximum integration error when computing the posterior distribution.}	
}
\details{
  The dks function performs the Bayesian and Frequentist diagnostic tests outlined
  in Leek and Storey (2009). The result of the function is a double Kolmogorov-Smirnov p-value as well
  as posterior probability of uniformity estimates for each of the studies. The p-values
  should be simulated from a realistic distribution and only the null p-values should
  be passed to the dks function. 
}

\value{
  \item{dkspvalue}{The double Kolmogorov-Smirnov p-value.}
  \item{postprob}{A B-vector of the posterior probability that each study's null p-values are uniform.} 
}
\references{
  J.T. Leek and J.D. Storey, "The Joint Null Distribution of Multiple Hypothesis Tests."
}

\author{Jeffrey T. Leek \email{jleek@jhsph.edu}}
\seealso{\code{\link{pprob.uniform}}, \code{\link{dks.pvalue}}, \code{\link{pprob.dist}},\code{\link{cred.set}}}
\examples{

  ## Load data
  data(dksdata) 
  
  ## Perform the diagnostic tests with plots
  dks1 <- dks(P)
  dks1$dkspvalue
}

\keyword{misc}
