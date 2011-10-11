\name{dks.pvalue}
\alias{dks.pvalue}
\title{Frequentist diagnostic test for multiple testing p-values.}
\description{
This function accepts a matrix of simulated null p-values where each column corresponds
to the p-values from a single simulated study. The null p-values should represent
a subset of all the simulated p-values corresponding to the tests with no signal.
}
\usage{
  dks.pvalue(P)
}
\arguments{
  \item{P}{An m0 x B matrix of null p-values, each column corresponds to the p-values from a single simulated study.}	
}
\details{
  The dks.pvalue function performs the double Kolmogorov-Smirnov test outlined in Leek and Storey (2009). 
  The p-values should be simulated from a realistic distribution and only the null p-values should
  be passed to the dks.pvalue function. 
}

\value{
  \item{dkspvalue}{The double Kolmogorov-Smirnov p-value.}
  \item{kspvalue}{A B-vector of the Kolmogorov-Smirnov p-values one for each test.} 
}
\references{
  J.T. Leek and J.D. Storey, "The Joint Null Distribution of Multiple Hypothesis Tests."
}

\author{Jeffrey T. Leek \email{jleek@jhsph.edu}}
\seealso{\code{\link{pprob.uniform}}, \code{\link{dks}}, \code{\link{pprob.dist}},\code{\link{cred.set}}}
\examples{
  ## Load data
  data(dksdata) 

  ## Calculate the double KS p-value
  dksp <- dks.pvalue(P)
  dksp$dkspvalue

  ## Histogram of the distribution of KS test p-values
  hist(dksp$kspvalue)
}

\keyword{misc}
