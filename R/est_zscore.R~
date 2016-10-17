##' Estimates the expected Z Score for a single SNP, assuming the input CVs, and the relationship, gamma, between them and the trait of interest
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
##' @author Mary Fortune
est_zscore<-function(N0,N1,Ufactor,powerfactor,freq,GenoProbXW){
  #Compute U
  #####
  N<-N0+N1
  #Contribution of each W=w to the Sum
  #P(X=1 AND W=w)
  PX1W<-GenoProbXW[[2]]
  #P(X=2 AND W=w)
  PX2W<-GenoProbXW[[3]]
  #Reduce compultational time marginally 
	UfactorPX2W<-sum(Ufactor*PX2W)
	UfactorPX1W<-sum(Ufactor*PX1W)
	powerfactorPX2W<-sum(powerfactor*PX2W)
	powerfactorPX1W<-sum(powerfactor*PX1W)
  # compute E(X=x|Y=y) 
  #find U
  U<-2*UfactorPX2W+UfactorPX1W
  # compute E(X^a) for various powers a
  EX<-2*powerfactorPX2W+powerfactorPX1W
  EX2<-4*powerfactorPX2W+powerfactorPX1W
  EX3<-8*powerfactorPX2W+powerfactorPX1W
  EX4<-16*powerfactorPX2W+powerfactorPX1W
  #find expected value of VX
  EVX<-EX2-(EX^2)
  #find expected value of VX^2
  Term1<-N*EX4+N*(N-1)*(EX2^2)
  Term2<-EX4+2*(N-1)*EX3*EX+(N-1)*(EX2^2)+(N-1)*(N-2)*EX2*(EX^2)
  Term3<-(EX4+4*(N-1)*EX3*EX+6*(N-1)*(EX2^2)+6*(N-1)*(N-2)*EX2*(EX^2)+(N-1)*(N-2)*(N-3)*(EX^4))/N
  EVX2<-(1/(N-1))^2 * ( Term1 - (2*Term2) + (Term3))
  #distribution of V
  a<-(2*EVX2-EVX*EVX)/(EVX2-EVX*EVX)
  b<-(EVX2*EVX)/(EVX2-EVX*EVX)
  EinvsqrtVX<-(b^-0.5)*exp(lgamma((2*a+1)/2)-lgamma(a))
  #Value of Z
  Z<-U*EinvsqrtVX*(N/(N0*N1))^0.5
  return(Z)
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


#wrapper function to run est_zscore for all snps in snps
#assumes we have a list, GenoProbList, giving the GenoProb values for each X. 
##' @title estimate Z score at a single SNP
##' @export
##' @param N0 The number of Y=0
##' @param N1	The number of Y=1
##' @param snps The snps at which we wish to compute the expected Z Score
##' @param W	The true causal SNPs (these need not be in "snps")
##' @param gamma	The odds ratios of effect of the true causal SNPs (including gamma0, the intercept term)
##' @param freq Frequencies of SNP appearances (computed using snphap)
##' @param GenoProbList An list of objects giving the probability of seeing each {X,W} genotype vector
##' @return The expected Z Score for SNP X, assuming the causal SNPs are W
##' @author Mary Fortune
est_statistic<-function(N0,N1,snps,W,gamma,freq,GenoProbList){
  #check that we have SNPs X and W in the reference dataset
  if (length(which(!(c(snps,W) %in% colnames(freq))))>0 ){stop("SNPs of interest not present in reference dataset.")}
  Est_stat<-rep(0,length(snps))
  # compute P(Y=1 | W=w)
  N<-N0+N1
  expeta<-exp(gamma[1]+rowSums(sweep((hcube(rep(3,length(W)))-1),MARGIN=2,gamma[-1],`*`)))
  #compute the constant factors we will multiply by
  Ufactor<-N0*(N-1)*(N0*expeta-N1)/(N^2)
  powerfactor<-N0*(expeta+1)/N
  for (ii in 1:length(snps)){
    Est_stat[ii]<-est_zscore(N0,N1,Ufactor,powerfactor,freq,GenoProbList[[ii]])
  }
  return(Est_stat)
}
