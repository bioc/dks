
##
## Calculate the posterior probability a sample comes from the uniform distribution
##

pprob.uniform <- function(p,alpha=c(0.1,10),beta=c(0.1,10),eps=1e-10){
  p <- pmax(pmin(p,1-1e-11),1e-11)
  m <- length(p)
  s1 <- sum(log(p))
  s2 <- sum(log(1-p))
  f <- function(x){
    return(exp(m*log(gamma(x[1] + x[2])/(gamma(x[1])*gamma(x[2]))) + s1 * (x[1]-1) + s2*(x[2]-1)))
  }
  xx <- f(beta.mom(p))
  if(is.infinite(xx)){
    pp <- 0
    return(pp)
  }else{
    pp <- 1/2/(1/2 + 1/2*adaptIntegrate(f=f,lowerLimit=c(alpha[1],beta[1]),upperLimit=c(alpha[2],beta[2]),tol=1e-7,fDim=1)$integral)
  }
  return(pp)
}

##
## Calculate the method of moments estimates for the beta distribution
##

beta.mom <- function(p){
  m <- mean(p)
  v <- 1/length(p)*sum((p-m)^2)
  alpha <- m*((m*(1-m))/v - 1)
  beta <- (1-m)*((m*(1-m))/v - 1)
  return(c(alpha,beta))
}


##
## Calculate the double KS p-value
##

dks.pvalue <- function(P){
  m <- dim(P)[1]
  B <- dim(P)[2]
  ksp <- apply(P,2,function(x){ks.test(x,"punif")$p.value})
  dksp <- ks.test(ksp,"punif")$p.value
  return(list(dkspvalue=dksp,kspvalue=ksp))
}


##
## Calculate the posterior distribution
##
##

pprob.dist <- function(p,alpha=c(0.1,10),beta=c(0.1,10),delta=0.10,eps=1e-10){
  p <- pmax(pmin(p,1-1e-11),1e-11)
  m <- length(p)
  s1 <- sum(log(p))
  s2 <- sum(log(1-p))
  f <- function(x){
    return(exp(m*log(gamma(x[1] + x[2])/(gamma(x[1])*gamma(x[2]))) + s1 * (x[1]-1) + s2*(x[2]-1)))
  }
  nconst <- adaptIntegrate(f=f,lowerLimit=c(alpha[1],beta[1]),upperLimit=c(alpha[2],beta[2]),tol=1e-7,fDim=1)$integral
  alpha <- seq(alpha[1],alpha[2],by=delta)
  beta <- seq(beta[1],beta[2],by=delta)
  est <- matrix(0,nrow=length(alpha),ncol=length(beta))
  for(i in 1:length(alpha)){
    for(j in 1:length(beta)){
      est[i,j] <-  f(c(alpha[i],beta[j]))
    }
  }
  dist <- est/nconst
  return(dist=est/nconst)
}



##
## Calculate a credible set
##
##

cred.set <- function(dist,delta=NULL,level=0.95){
  if(is.null(delta)){stop("Grid width must be specified")}
  m <- dim(dist)[1]
  xx <- which.max((cumsum(rev(sort(dist*delta^2))) > level))
  cred <- (dist*delta^2 >= rev(sort(dist*delta^2))[xx])
  elevel <-sum(dist[cred]*delta^2)
  return(list(cred=cred,elevel=elevel,level=level))
}


##
## Global Function
##
##

dks <- function(P,alpha=c(0.1,10),beta=c(0.1,10),plot=TRUE,eps=1e-10){
  if(plot){cat("\nComputing DKS P-value...")}
  dksp <- dks.pvalue(P)
  if(plot){cat("Computing Posterior Probabilities...")}
  B <- ncol(P)
  m <- nrow(P)
  pp <- rep(0,B)
  for(i in 1:B){pp[i] <- pprob.uniform(P[,i],alpha=alpha,beta=beta,eps=eps)}
  if(plot){
    par(mfrow=c(1,3),mar=c(5,5,5,5))
    plot(1:100/101,1:100/101,xlab="Uniform Quantiles",ylab="Empirical Quantiles",type="n",cex.lab=1.5,xlim=c(0,1),ylim=c(0,1))
    for(i in 1:B){
      lines((1:m)/(m+1),sort(P[,i]),col="lightblue")
    }
    abline(c(0,1),col="black",lty=2)
    plot(sort(dksp$kspvalue),1:B/B,xlab="KS P-value",ylab="Empirical CDF",lwd=3,col="grey",type="n",cex.lab=1.5,xlim=c(0,1),ylim=c(0,1),main=paste("DKS P-value: ",round(dksp$dkspvalue,3)))
    cdfplot(sort(dksp$kspvalue),1:B/B,col="lightblue",lwd=3)
    abline(c(0,1),col="black",lty=2)
    crit.val <- qksone(0.95)$root
    ub <- pmin(seq(1,B)/B + crit.val/sqrt(B), 1)
    lb <- pmax(seq(1,B)/B - crit.val/sqrt(B), 0)
    cdfplot(sort(dksp$kspvalue),ub,col="grey",lwd=3)
    cdfplot(sort(dksp$kspvalue),lb,col="grey",lwd=3)
    hist(pp,col="lightblue",border="grey",xlim=c(0,1),xlab="Posterior Probability Uniform",ylab="Frequency",main="")
  }
  return(list(dkspvalue=dksp$dkspvalue,postprob=pp))
}


cdfplot<- function(x,y,lty=2,...){
  m <- length(x)
  for(i in 1:(m-1)){
    segments(x[i],y[i],x[(i+1)],y[i],lty=lty,...)
    segments(x[(i+1)],y[i],x[(i+1)],y[(i+1)],lty=lty,...)
  }
}


pksone <- function(x) {
  k <- seq(1:20)
  1 + 2 * sum( (-1)^k * exp(- 2 * k^2 * x^2))
}
qksone <- function(p) {
  foo <- function(x) pksone(x) - p
  uniroot(foo, c(0.5, 10))
}
