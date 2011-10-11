\name{cred.set}
\alias{cred.set}
\title{Calculate a credible set for the posterior distribution on the Beta hyperparameters.}
\description{
This function accepts a distribution calculated with pprob.dist and calculates
a credible set of the specified level for the hyperparameters. If the credible set
includes the value (1,1) the sample is likely to be uniform. 
}
\usage{
  cred.set(dist,delta=NULL,level=0.95)
}
\arguments{
  \item{dist}{The posterior distribution for the hyperparameters computed with pprob.dist.}
  \item{delta}{The grid size, must match the grid size from pprob.dist.}
  \item{level}{The level of the credible set.}	
}
\details{
  The cred.set function calculates a credible set of the specified level based on the
  distribution calculated with pprob.dist. The grid size, delta, should match the grid
  size from the call to pprob.dist. The result is a matrix of the same size as dist
  which indicates whether each point is in the credible set. 
}

\value{
  \item{cred}{The credible set for the hyper-parameters of the beta distribution.}
  \item{level}{The user specified level of the set.} 
  \item{elevel}{The empirical level of the set, the smaller delta is, the closer elevel will be to level.}
}
\references{
  J.T. Leek and J.D. Storey, "The Joint Null Distribution of Multiple Hypothesis Tests."
}

\author{Jeffrey T. Leek \email{jleek@jhsph.edu}}
\seealso{\code{\link{dks}}, \code{\link{dks.pvalue}}, \code{\link{pprob.dist}},\code{\link{cred.set}}}
\examples{
  ## Load data
  data(dksdata) 

  ## Calculate the posterior distribution
  dist1 <- pprob.dist(P[,1])

  delta = 0.1
  ## Calculate a 95% credible set
  cred1 <- cred.set(dist1,delta=0.1)

  ## Plot the posterior and the credible set
  
  alpha <- seq(0.1,10,by=delta)
  beta <- seq(0.1,10,by=delta)

  par(mfrow=c(1,2))
  image(log10(alpha),log10(beta),dist1,xaxt="n",yaxt="n",xlab="Alpha",ylab="Beta")
  axis(1,at=c(-2,-1,0,1,2),labels=c("10^-2","10^-1","10^0","10^1","10^2"))
  axis(2,at=c(-2,-1,0,1,2),labels=c("10^-2","10^-1","10^0","10^1","10^2"))
  points(0,0,col="blue",cex=1,pch=19)	

  image(log10(alpha),log10(beta),cred1$cred,xaxt="n",yaxt="n",xlab="Alpha",ylab="Beta")
  axis(1,at=c(-2,-1,0,1,2),labels=c("10^-2","10^-1","10^0","10^1","10^2"))
  axis(2,at=c(-2,-1,0,1,2),labels=c("10^-2","10^-1","10^0","10^1","10^2"))
  points(0,0,col="blue",cex=1,pch=19)	
}

\keyword{misc}
