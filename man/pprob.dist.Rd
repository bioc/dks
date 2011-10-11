\name{pprob.dist}
\alias{pprob.dist}
\title{The posterior distribution for the hyper-parameters of the Beta distribution.}
\description{
This function accepts a vector of simulated null p-values from a single simulated study. 
The null p-values should represent a subset of all the simulated p-values corresponding 
to the tests with no signal. The result is an estimated posterior distribution for the
parameters of the Beta distribution. A posterior centered at (1,1) suggests a uniform
distribution.
}
\usage{
  pprob.dist(p,alpha=c(0.1,10),beta=c(0.1,10),delta=0.10,eps=1e-10)
}
\arguments{
  \item{p}{An vector of null p-values from a single simulated study.}
  \item{alpha}{The range of the first parameter for the prior on the beta distribution.}
  \item{beta}{The range of the second parameter for the prior on the beta distribution.}
  \item{delta}{The grid size, the posterior is calculated over the range of the parameters at grid points separated by delta.}
  \item{eps}{Maximum integration error when computing the posterior distribution.}	
}
\details{
  The pprob.dist function calculates the posterior probability for the parameters of 
  the beta distribution given the sample p. The prior is assumed to be uniform on 
  the range specified by the user. A posterior distribution is returned in the form
  of a matrix, where element (i,j) is the posterior at (alpha[1] + i*delta, beta[1] + j*delta).
  The null p-values should be simulated from a realistic distribution and only the null p-values should
  be passed to the pprob.dist function. 
}

\value{
  \item{dist}{The posterior distribution in the form of a matrix.}
}
\references{
  J.T. Leek and J.D. Storey, "The Joint Null Distribution of Multiple Hypothesis Tests."
}

\author{Jeffrey T. Leek \email{jleek@jhsph.edu}}
\seealso{\code{\link{dks}}, \code{\link{dks.pvalue}}, \code{\link{pprob.uniform}},\code{\link{cred.set}}}
\examples{

  ## Load data
  data(dksdata) 

  ## Calculate the posterior distribution
  dist1 <- pprob.dist(P[,1])

  delta <- 0.1

  ## Plot the posterior distribution
  alpha <- seq(0.1,10,by=delta)
  beta <- seq(0.1,10,by=delta)
  image(log10(alpha),log10(beta),dist1,xaxt="n",yaxt="n",xlab="Alpha",ylab="Beta")
  axis(1,at=c(-2,-1,0,1,2),labels=c("10^-2","10^-1","10^0","10^1","10^2"))
  axis(2,at=c(-2,-1,0,1,2),labels=c("10^-2","10^-1","10^0","10^1","10^2"))
  points(0,0,col="blue",cex=1,pch=19)	

}

\keyword{misc}
