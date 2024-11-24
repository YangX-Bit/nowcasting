######################################################################
# File containing helper functionality to show the predictive PMF
# for N(t,Inf) for a given t based on information available
# at time T. C.f. Fig. 4 in the manuscript.
######################################################################

#Necessary libraries
require("RColorBrewer")

"%without%" <- function(x, y) x[!x %in% y] #--  x without y


######################################################################
# Plot the PMF of all forecast distributions contained in the
# nowcast object
#
# Parameters:
#  nc - nowcast object (class stsBP)
#  s  - "now"
#  t  - which time point
#  conf.level - confidence level for the confidence intervals to plot
#
# Note: The function outputs some values on the screen.
######################################################################

plotPMF <- function(nc,s,t,xlim=NULL,conf.level=0.95,without="",xlab="Hospitalizations",main=paste("now=",s,", when=",t,sep="")) {
  #Check that input is valid.
  if ((length(s) != 1) | (length(t) != 1)) {
    stop("Plotting only works for a single date s.")
  }
    
  i <- which(epoch(nc) == t)
  yts <- as.numeric(observed(nc)[i,])
  sts.truth <- nc@truth
  ytinf <- observed(sts.truth)[i,]
  
  yt.support <- nc@control$yt.support
  y.prior.max <- nc@control$y.prior.max
  
  cat("Observed value at y_{t,s}=",yts,"\n")
  cat("True value y_t =", ytinf,"\n")

  alpha <- 1-conf.level
  
  #summaries of the posterior
  summary.post <- function(post.pmf) {
    cat(conf.level*100, " equal tailed CI:\n")
    print(yt.support[c(which.max(cumsum(post.pmf) > alpha/2),which.max(cumsum(post.pmf) > 1-alpha/2))])
    cat("Mean, median and mode of the posterior: ",
        c(mean.post=sum(0:y.prior.max*post.pmf),
          median.post=which.max(cumsum(post.pmf)>0.5),
          mode.post=which.max(post.pmf)),"\n")
    invisible()
  }
  
#  cat(paste("--------- s=",s," t=",t," Delta=",timeDelay(t,s),"-----------\n",sep=""))
  PMF <- nc@predPMF[[as.character(t)]]
  which <- names(PMF) %without% without
  
  for (n in names(PMF[which])) {
    cat("-------------\n",n,":\n-------------\n")
    summary.post(PMF[[n]])
  }

  #Last observation with a mass bigger than 1e-5
  y.boring <- min(unlist(lapply(PMF[which], function(x) which.max( rev(x)>1e-5))))
  xlim.max <- y.prior.max - y.boring + 1
  if (is.null(xlim)) {
    xlim <- c(0,xlim.max)
  }

  ##   #Show prior of bayesian procedures
  ## if (all(c("mean.lambda","var.lambda") %in% names(attributes(nc@control$N.tInf.prior)))) {
  ##   mean <- attr(  nc@control$N.tInf.prior,"mean.lambda")
  ##   size <- mean^2/(attr(  nc@control$N.tInf.prior,"var.lambda") - mean)
  ##   nbprior <- dnbinom(nc@control$N.tInf.support,mu=mean,size=size)
  ## } else {
  ##   nbprior <- NULL
  ## }
  nbprior <- NULL
  #choose palette
  nLines <- length(which) + (!is.null(nbprior))
  if (getOption("color")) {
    pal <- brewer.pal(n=nLines,"Set1")
    lty <- rep(1,nLines)
  } else {
    pal <- grey(seq(0.1,0.6,length=nLines))
    lty <- seq_len(nLines)
  }
  

  #Create a plot
  plot(yt.support,PMF[[which[1]]],type="l",xlab=xlab,ylab="PMF",lwd=2,ylim=c(0,max(unlist(PMF[which]),na.rm=TRUE)),xlim=xlim,main=main,col=pal[1],lty=lty[1])
  #Add others
  for (i in 1+seq_len(length(which[-1]))) {
    lines(yt.support,PMF[[which[i]]],col=pal[i],lwd=2,lty=lty[i])
  }
  #Add prior if its there
  if (!is.null(nbprior)) {
    lines(nc@control$N.tInf.support,nbprior,col=pal[length(which)+1],lwd=2,lty=lty[length(which)+1])
  }
  
  #Show current available data
  #axis(1,at=yts,label=expression(y["t,s"]),col="black")
  axis(1,at=yts,label="",col="black",tcl=-0.8,lwd=3)
  axis(1,at=yts,label=expression(N(t,T)),col="black",line=1.2,tick=FALSE)
  legend(x="topright",c(which,if (!is.null(nbprior)) "prior" else NULL),col=pal,lwd=2,lty=lty)

  #Show truth
  axis(1,at=ytinf,label="",col="black",tcl=-0.8,lwd=3)
  axis(1,at=ytinf,label=expression(N(t,infinity)),col="black",line=1.2,tick=FALSE)
  lines(rep(ytinf,2),c(0,1e99),lty=2)
  invisible()
}
