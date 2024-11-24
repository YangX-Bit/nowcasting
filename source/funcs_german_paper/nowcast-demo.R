######################################################################
# This file contains code to reproduce some of the analyses and
# graphics from the 2014 Biometrics Manuscript 
# Bayesian Nowcasting during the STEC O104:H4 Outbreak in Germany, 2011
# by Höhle, M and an der Heiden, M.
#
# Author: Michael Höhle <http://www.math.su.se/~hoehle>
# Date:   12 Apr 2014
#
# Disclaimer: The code is available under the GNU GPL 2.0 license.
# Future developments of the code and plotting routines will
# occur within the surveillance R package. See the help of the
# function 'nowcast' in the package for more details.
######################################################################

#Load surveillance package. The nowcast functionality is available
#from version 1.8-0. Currently, this version is not on CRAN but only
#on R-forge -- see http://surveillance.r-forge.r-project.org for
#details on how to install this version.
library("surveillance")
library("rjags")


#Set working directoy to the directory where this file resides. Often automatic.
setwd("/Users/xiaoy0a/Desktop/Task/SecondProject/biom12194-sm-0001-suppdatacode/R")

#Source in extra code. This will go into the surveillance package
#in the future. Currently the functions are somewhat specific
#to the STEC O104:H4 application.
source("../R/animate-nowcast.R")
source("../R/plotPMF.R")
source("../R/plotMedianDelay.R")
source("../R/plotReportingTriangle.R")

######################################################################
#Load the data, which are courtesy of the Robert Koch Institute,
#Berlin, Germany.
#
#Note: The data are such that the resulting reporting triangle
#corresponds to Fig.1 of the Web appendix. This means that the reports
#which arrived with a delay longer than 15 days are set to have have
#arrived after 15 days. Altogether, this gives small discrepancies
#with the results in the manuscript. However, as mentioned in the paper,
#longer delays were not very relevant for the nowcasting.
######################################################################
data("husO104Hosp")

######################################################################
#Set some variables for the processing
######################################################################
options(color=TRUE) #Say output of figures can be in color (b/w for paper)


#======================================================================
# Begin nowcasting demo
#======================================================================

#Fix maximum delay
D <- 15 

#Extract the reporting triangle at a specific day
t.repTriangle <- as.Date("2011-07-04")
nc <- nowcast(now=t.repTriangle,when=t.repTriangle,dEventCol="dHosp",dReportCol="dReport",data=husO104Hosp,D=D,method="lawless")

?nowcast
nc_2 <- nowcast(now=t.repTriangle,when=as.Date("2011-06-15"),dEventCol="dHosp",dReportCol="dReport",data=husO104Hosp,D=D,method="lawless")
plotReportingTriangle(nc_2)
#Show the reporting triangle (this is Fig. 1 from the Web Appendix)
plotReportingTriangle(nc)

#Extract it manually and show last 20 days of the outbreak
tail(reportingTriangle(nc),n=20)

######################################################################
#Do the calculations of the delay distribution using the equations
#of Lawless (1994). This corresponds to Table 1 from the Web Appendix
#of Hoehle and an der Heiden (2014)
######################################################################

#Extract aggregates from reporting triangle
n.x <- matrix(attr(reportingTriangle(nc),"n.x"),ncol=1)
N.x <- matrix(attr(reportingTriangle(nc),"N.x"),ncol=1)
#Compute MLEs of reverse time hazards
g.hat <- ifelse( !is.na(n.x/N.x), n.x/N.x, 0)
#Computing resulting probabilities
pd <- g.hat[(0:D)+1] * c(sapply(0:(D-1), function(d) prod(1-g.hat[(d+1):D+1])),1)
#Make a table with the numbers (corresponds to Tab. 1 in the Web Appendix)
data.frame(n.x=n.x,N.x,g.hat=g.hat,pd=pd,row.names=as.character(0:D))

######################################################################
# Show predictive PMF for a specific date  using several different
# nowcast procedures
######################################################################

#Setup the data
now <- as.Date("2011-06-02")
when <- as.Date("2011-05-28")
#Dates where the corresponding sts object should be available
dateRange <- seq(as.Date("2011-05-01"),as.Date("2011-07-06"),by="1 day")

#Setup up the control object for the different nowcast procedures
nc.control <- list(predPMF  =TRUE,  #compute and store the predictive PMF for each time point and method
                   N.tInf.max=300,  #Support of each N.tInf is 0,\ldots,N.tInf.max
                   score=F,      #Compute scoring rules for each predictive PMF (FALSE is faster)
                   dRange=dateRange,#Dates which constitutes time points 1,...,length(dateRange)
                   #Specification of the hierarchical Bayes model
                   ddcp=list(ddChangepoint=as.Date(c("2011-05-23")),
                       logLambda="tps",
                       tau.gamma=1,
                       mcmc=c(burnin=1000,sample=10000,thin=1,
                              adapt=1000,store.samples = F)),
                   #Specification for bayes.notrunc and bayes.trunc (especially prior)
                   N.tInf.prior=structure("poisgamma",mean.lambda=50,var.lambda=3000))
                  

#Nowcast using a set of methods (no moving window)
nc2 <- nowcast(now=now,when=when,dEventCol="dHosp",dReportCol="dReport",data=husO104Hosp,
               D=D,m=NULL,
               #method=c("lawless","bayes.notrunc","bayes.trunc","bayes.trunc.ddcp"),
               method=c("bayes.trunc.ddcp"),
               control=nc.control)

#Compare the predictive distributions for N(t=when,Inf) between procedures
plotPMF(nc2,s=now,t=when,main=when,xlab="")

#Show CDF of delay distribution for bayes.trunc
delayCDF(nc2)[["bayes.trunc"]]
#For bayes.trunc.ddcp the CDF changes over time. The delayCDF is the posterior median, but
#all specific of the model are included in the attribute "model"
class(delayCDF(nc2)[["bayes.trunc.ddcp"]])
attr(delayCDF(nc2)[["bayes.trunc.ddcp"]],"model")

#Show empirical delay distribution (median and quantiles) together with the fit
#of the bayes.trunc.ddcp procedure based on the data available up to "now" of nc2.
plotMedianDelay(nc2,w=2,rT.truth=reportingTriangle(nc), dates=attr(reportingTriangle(nc2),"t02s"))


#Scoring info: Only available if nc.control$score == TRUE
dimnames(score(nc2))
#Show scores for the single time point
score(nc2)[as.character(when),,]


######################################################################
# Nowcasting an entire sequence of dates (we choose only
# one fast method here for illustration purposes) and animate
# result.
######################################################################

#=============================
#Set up some variables for use
#=============================

#Which method to use, could e.g. also be "bayes.trunc.ddcp" (which is much slower!)
method <- "bayes.trunc.ddcp" 
#Number of nowcasts to do back in time
k <- 10                 
#First safe time point back in time to do nowcasting for
safePredictLag <- 3     
#Range of values to do the scoring for.
scoreRange <- seq(as.Date("2011-06-02"),as.Date("2011-07-04"),by="1 day")

#=============================================
#Do nowcasting for an entire set of timepoints
#=============================================

nowcastList <- list()
#Nowcast all time point within the specified range. This might take a while (!)
for (i in 1:length(scoreRange)) {
  #What's "today"
  now <- scoreRange[i]
  #Show some status information
  cat(paste("====================\nnow=",now," (",i,"/",length(scoreRange),")\n====================\n",sep=""))
  #Which time points to do the nowcast for
  when <- seq(now-k-safePredictLag+1, now-safePredictLag, by="1 day")
  #Nowcast
  # nc.one <- nowcast(now=now,when=when,dEventCol="dHosp",dReportCol="dReport",
  #                   data=husO104Hosp,D=D,m=NULL,method=method,control=nc.control)
  # nowcastList[as.character(now)] <- nc.one
  print(when)
}


#Make an animation similar to Web Animation 3 of the manuscript.
#Other output options (html or shockwave movie are possible using the 'animation' package)
getwd()  #Show current working directory (has to be a writeable directory)
fileName <- "NowcastAnimation.pdf" #File name, which will be stored in the current directory
#Note: On windows it is not possible to modify a file while it is open in some
#other program.

pdf(fileName,width=8,height=5,onefile=TRUE)
animate.nowcast(linelist=husO104Hosp,
                dEventCol="dHosp",dReportCol="dReport",
                aggregate.by="1 day",
                nowcasts=nowcastList,
                method=method,
                control=list(dRange=range(dateRange),
                  anim.dRange=range(scoreRange),
                  plot.dRange=as.Date(c("2011-05-01","2011-07-04")),
                  consistent=FALSE,hookFun=function(range, ymax) {},
                  sys.sleep=0,
                  safeguard=safePredictLag,cex.names=0.7))
dev.off()

#-> You can open resulting PDF file in Acrobat Reader or similar PDF file viewer
#and then scroll through the pages (one page equals one day).


######################################################################
#Do just one time point with bayes.trunc.ddcp (it's rather slow).
######################################################################

now <- as.Date("2011-06-02")
method <- "bayes.trunc.ddcp"
when <- seq(now-k-safePredictLag+1, now-safePredictLag, by="1 day")
nowcastList[as.character(now)] <- nc.one <- nowcast(now=now,when=when,dEventCol="dHosp",dReportCol="dReport",data=husO104Hosp,D=D,m=NULL,method=method,control=nc.control)

#Show the result for the one time point
animate.nowcast(linelist=husO104Hosp,
                dEventCol="dHosp",dReportCol="dReport",
                aggregate.by="1 day",
                nowcasts=nowcastList,
                method=method,
                control=list(dRange=range(dateRange),
                  anim.dRange=c(now,now),  #just one time point
                  plot.dRange=as.Date(c("2011-05-01","2011-07-04")),
                  consistent=FALSE,hookFun=function(range, ymax) {},
                  sys.sleep=0,
                  safeguard=safePredictLag,cex.names=0.7))


#======================================================================
# End nowcasting demo
#======================================================================

