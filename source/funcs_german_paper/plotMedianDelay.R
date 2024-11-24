######################################################################
# Show empirical delay distribution as a function of occurence time t.
# In case a discrete time hazard model was fit the posterior median
# of this fit can be illustrated as well.
#
# Author: Michael Höhle <http://www.math.su.se/~hoehle>
# Date:   Dec 2013
######################################################################


######################################################################
#Convert discrete time hazards to PMF
#Parameters:
#  haz - vector with entries for (0,...,Dmax)
######################################################################

haz2pmf <- function(haz) {
    PMF <- 0*haz
    for (i in 0:(length(haz)-1)) {
        PMF[i+1] <- haz[i+1] * (1-sum(PMF[seq(i)]))
    }
    return(PMF)
}

######################################################################
# Find a quantile of a discrete random variable with support on
# 0,...,D and which has a PMF given by the vector prob. We
# define the q quantile as \min_{x} F(x) \geq q.
#
# Parameters:
#   prob - vector on 0,..,D containing the PMF
#   q    - quantile to compute
######################################################################

pmfQuantile <- function(prob,q=0.5) {
    which.max(cumsum(prob) >= q)-1
}


######################################################################
# Show empirical and model based median of delay distribution as a
# function of occurence time t.
#
# Parameters:
#  w  - half-width of moving window
#  rT.truth - reporting triangle as it would be at the end
#  date - vector of dates where to show the result
#  modelQuantiles - which model quantiles to show
######################################################################

plotMedianDelay <- function(nc, w=1, rT.truth, dates, modelQuantiles=0.5) {
    #Function requires bayes.trunc.ddcp output in order to show the model output
    #Possibly change code so this is not necessary anymore.
    if (! ("bayes.trunc.ddcp" %in% names(delayCDF(nc)))) {
      stop("No bayes.trunc.ddcp output present in nc object.")
    }

    #Determine max delay from reporting triangle.
    D <- ncol(rT.truth) - 1
    res <- matrix(NA, nrow=length(dates), ncol=D+1)

    #which data variables are actually in rT.truth
    isThere <- !is.na(sapply(dates, function(date) pmatch(as.character(date),rownames(rT.truth))))
    idx <- which(isThere)
    
    #Loop over all time points.
    for (i in (w+min(idx)):(max(idx)-w)) {
        now <- dates[i]
        idx <- pmatch(as.character(now),rownames(rT.truth))
        subset <- rT.truth[idx + c(-w:w),,drop=FALSE] 
        res[i,] <- colSums(subset) / sum(subset)
    }

    #A slightly modified function to determine quantiles, which can
    #handle NAs (if no case at all)
    quantile <- function(q) {
        apply(res, 1, function(x) {
            if (all(is.na(x))) return(NA) else return(which.max(cumsum(x) >= q) - 1)
        })
    }

    #Find 10%, 50% and 90% quantiles
    quants <- sapply(c(0.1,0.5,0.9), quantile)

    
    ## #Show model based estimates
#    ddcp.model <- attr(delayCDF(nc)[["bayes.trunc.ddcp"]],"model")
#    post.median <- ddcp.model$post.median
#    W <- ddcp.model$W
    
    ## #Extract parameters
    ## gamma.red <- post.median[grep("gamma",names(post.median))]
    ## eta <- matrix(post.median[grep("^eta",names(post.median))])
    ## #Map from reduced set to full set
    ## gamma <- gamma.red[round( (0:(D-1))/2 - 0.4) + 1]

    #Prepare result of model estimated median
#    quants.model <- matrix(NA, nrow=length(dates),ncol=length(modelQuantiles),dimnames=list(as.character(dates),modelQuantiles))
#    pmf         <- matrix(NA, nrow=length(dates),ncol=D+1,dimnames=list(as.character(dates),paste("delay",0:D,sep="")))

    ## #Determine median if part of model
    ## for (t in 1:length(dates)) {
    ##     if (as.character(dates[t]) %in% rownames(W)) {
    ##         lin.pred <- ( gamma + t(eta) %*% W[t,,0:D])
    ##         pmf[t,] <- haz2pmf(c(plogis(lin.pred),1))
    ##         quants.model[t,] <- sapply(modelQuantiles, function(q) pmfQuantile( pmf[t,],q=q))
    ##     }
    ## }

    #Alternative, because the CDF is already there
    
    cdf <- cbind(0,delayCDF(nc)[["bayes.trunc.ddcp"]])
    pmf <- t(apply(cdf,1,diff))

    #Determine quantiles
    quants.model <- matrix(NA, nrow=length(dates),ncol=length(modelQuantiles),dimnames=list(as.character(dates),modelQuantiles))
    for (t in 1:length(dates)) {
      quants.model[t,] <- sapply(modelQuantiles, function(q) pmfQuantile( pmf[t,],q=q))
    }
      
    #Make sure the NAs in the beginning agree
    i <- 1
    while (all(is.na(quants[i,]))) {quants.model[i,] <- NA ; i <- i + 1}
    
    #Make a plot (use plot.Dates instead of matplot)
    plot(dates, quants[,2],xlab="Time of occurence",ylab="delay (days)",ylim=c(0,15),type="l",col=1,lty=c(1),lwd=4)
    matlines(dates, quants[,c(1,3)],type="l",col=1,lty=c(2,3),lwd=c(1,1))
    matlines(dates, quants.model, col="gray",lwd=ifelse(modelQuantiles==0.5,2,1),lty=ifelse(modelQuantiles==0.5,1,2))
    
    #Show lines for breakpoints (if available)
    ddcp.model <- attr(delayCDF(nc)[["bayes.trunc.ddcp"]],"model")
    changePoints <- as.Date(colnames(ddcp.model$W))    
    for (i in 1:length(changePoints)) {
        axis(1,at=changePoints[i], changePoints[i], las=1, cex.axis=0.7,line=-2.5)
        lines( rep(changePoints[i],2),c(0,1e99),lty=2)
    }

    #Make a legend
    legend(x="bottomleft",c(expression(q[0.1](T)),expression(q[0.5](T)),expression(q[0.9](T)),expression(q[0.5]^"ddcp"(T))),lty=c(2,1,3,1),col=c(1,1,1,"gray"),lwd=c(1,4,1,2))

    #Add title
    title(nc@control$now)
    
    #Done
    invisible(pmf)
}
