##' Estimates the expected Z Score for a single SNP
##'
##' Assumes the input CVs, and the relationship, gamma, between them and the trait of interest
##'
##' Assumes we have already generated GenoProbXW for all X
##' @title estimate Z score at a single SNP
##' @export
##' @param N0 The number of Y=0
##' @param N1	The number of Y=1
##' @param Ufactor	The constant factor used to compute the expectation of U
##' @param powerfactor	The constant factor used to compute the expectation of the genotype of X to some power
##' @param freq Frequencies of SNP appearances (computed using snphap)
##' @param GenoProbXW An object giving the probability of seeing each {X,W} genotype vector
##' @return The expected Z Score for SNP X, assuming the causal SNPs are W
##' @author Mary Fortune and Chris Wallace
est_zscore<-function(N0,N1,Ufactor,powerfactor,freq,GenoProbXW){
    ## uses Rcpp file ../src/est_zscore.cpp
  N<-N0+N1
  #Contribution of each W=w to the Sum
  #P(X=1 AND W=w)
  PX1W<-GenoProbXW[[2]]
  #P(X=2 AND W=w)
  PX2W<-GenoProbXW[[3]]
    zscore(N0,N1,Ufactor,powerfactor,PX1W,PX2W)
}


#returns P(X=x, W=w)
find_PXaW_MK<-function(x,w,GenoProbXW){
  XW<-paste0(x,paste0(w,collapse=""),collapse="")
  if (XW %in% names(GenoProbXW)){
    return(GenoProbXW[XW])
  }else{
    return(0)
  }	
  
}


##' Wrapper function to run est_zscore for all snps in snps
##'
##' Assumes we have a list, GenoProbList, giving the GenoProb values for each X. 
##' @title estimate Z score at a single SNP
##' @export
##' @param N0 The number of Y=0
##' @param N1	The number of Y=1
##' @param snps The snps at which we wish to compute the expected Z Score
##' @param W	The true causal SNPs (these need not be in "snps")
##' @param gamma1	The log odds ratios of effect of the true causal SNPs (not including gamma0, the intercept term)
##' @param freq Frequencies of SNP appearances (computed using snphap)
##' @param GenoProbList An list of objects giving the probability of seeing each {X,W} genotype vector
##' @return The expected Z Score for SNP X, assuming the causal SNPs are W
##' @author Mary Fortune and Chris Wallace
est_statistic<-function(N0,N1,snps,W,gamma1,freq,GenoProbList){
                                        #check that we have SNPs X and W in the reference dataset
    if (!all(c(snps,W) %in% colnames(freq)))
        stop("SNPs of interest not present in reference dataset.")
    if(length(gamma1)!=length(W))
        stop("length mismatch: gamma1 and W")
    if(length(GenoProbList)!=length(snps))
        stop("GenoProbList should have same length and order as snps")
    g0 <- compute_gamma0(N0=N0,N1=N1,W=W,gamma.CV=gamma1,freq=freq)
    ## compute P(Y=1 | W=w)
    N<-N0+N1
    expeta<-exp(g0+rowSums(sweep((hcube(rep(3,length(W)))-1),MARGIN=2,gamma1,`*`)))
                                        #compute the constant factors we will multiply by
    Ufactor<-N0*(N-1)*(N0*expeta-N1)/(N^2)
    powerfactor<-N0*(expeta+1)/N
    sapply(seq_along(snps), function(ii) {
        est_zscore(N0,N1,Ufactor,powerfactor,freq,GenoProbList[[ii]])
    })
}


