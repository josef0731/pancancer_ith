qvalue <- function(p, alpha=NULL, lam=NULL, robust=F)
{
#This is a function for estimating the q-values for a given set of p-values. The
#methodology comes from a series of recent papers on false discovery rates by John
#D. Storey et al. See http://www.stat.berkeley.edu/~storey/ for references to these
#papers. This function was written by John D. Storey. Copyright 2002 by John D. Storey.
#All rights are reserved and no responsibility is assumed for mistakes in or caused by
#the program.
#
#Input
#=============================================================================
#p: a vector of p-values (only necessary input)
#alpha: a level at which to control the FDR (optional)
#lam: the value of the tuning parameter to estimate pi0 (optional)
#robust: an indicator of whether it is desired to make the estimate more robust 
#        for small p-values (optional)
#
#Output
#=============================================================================
#remarks: tells the user what options were used, and gives any relevant warnings
#pi0: an estimate of the proportion of null p-values
#qvalues: a vector of the estimated q-values (the main quantity of interest)
#pvalues: a vector of the original p-values
#significant: if alpha is specified, and indicator of whether the q-value fell below alpha 
#    (taking all such q-values to be significant controls FDR at level alpha)

#This is just some pre-processing
    m <- length(p)
#These next few functions are the various ways to automatically choose lam
#and estimate pi0
    if(!is.null(lam)) {
        pi0 <- mean(p>lam)/(1-lam)
        remark <- "The user prespecified lam in the calculation of pi0."
    }
    else{
        remark <- "A smoothing method was used in the calculation of pi0."
        lam <- seq(0,0.95,0.01)
        pi0 <- rep(0,length(lam))
        for(i in 1:length(lam)) {
        pi0[i] <- mean(p>lam[i])/(1-lam[i])
        }
        spi0 <- smooth.spline(lam,pi0,df=3,w=(1-lam))
        pi0 <- predict(spi0,x=0.95)$y
    }
#The q-values are actually calculated here
    u <- order(p)
    v <- rank(p)
    qvalue <- pi0*m*p/v
    if(robust) {
        qvalue <- pi0*m*p/(v*(1-(1-p)^m))
        remark <- c(remark, "The robust version of the q-value was calculated. See Storey JD (2002) JRSS-B 64: 479-498.")
    }
    qvalue[u[m]] <- min(qvalue[u[m]],1)
    for(i in (m-1):1) {
    qvalue[u[i]] <- min(qvalue[u[i]],qvalue[u[i+1]],1)
    }
#Here the results are returned
    if(!is.null(alpha)) {
        output = list(remarks=remark, pi0=pi0, qvalue=qvalue, significant=(qvalue <= alpha), pvalue=p)
    }
    else {
        output = list(remarks=remark, pi0=pi0, qvalue=qvalue, pvalue=p)
    }
    return(output)
}