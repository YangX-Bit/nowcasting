#Shows the rows and columns of the reporting trianle
plotReportingTriangle <- function(nc) {
    #Labelling of the axis - depends on day of the month
    range.dates <- attr(reportingTriangle(nc),"t02s")
    names.arg <- format(range.dates,"%a %d %b")
    #Maximum delay used in the calculation
    D <- attr(reportingTriangle(nc),"D")
    
    #Size of text
    theCex <- 0.9
    theFont <- 2
    


    #Make a plot of the triangle
    plot(NA,xlim=c(0,D+2),ylim=range(1:length(range.dates)),axes=FALSE,ylab="",xlab="Delay")
    axis(1, at = 0:D, cex.axis = theCex, font = theFont)  # x axis
    axis(2, at = 1:length(range.dates), label = names.arg, las = 1, cex.axis = theCex, font = theFont)  # y
    axis(1,at=D+2, lab=expression(Y(t,D)),cex.axis=theCex)
    lines( rep(D+1,2), c(-1e99,1e99),lty=2)
    
    #Add counts as text
    for (i in 1:length(range.dates)) {
        for (j in 0:D) {
            text(j,i,reportingTriangle(nc)[i,j+1],cex=theCex)
        }
        #Show the row sum
        text(D + 2, i, attr(reportingTriangle(nc),"N.tT")[i],cex=theCex)
    }

    invisible()
}

