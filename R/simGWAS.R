##' @title Compute expected Z Score
##' @export
##' @param N0 The number of Y=0
##' @param N1	The number of Y=1
##' @param snps The snps at which we wish to compute the expected Z Score
##' @param W	The true causal SNPs (these need not be in "snps")
##' @param gamma.CV	The odds ratios of effect of the true causal SNPs
##' @param freq Frequencies of SNP appearances (computed using snphap)
##' @return The expected Z Score for all snps in snps, assuming the causal SNPs are W
##' @author Mary Fortune
expected_z_score<-function(N0,N1,snps,W,gamma.CV,freq){
	GenoProbList<-make_GenoProbList(snps,W,freq)
	gamma.sim<-c(compute_gamma0(N0,N1,W,gamma.CV,freq),gamma.CV)
	exp_z_score<-est_statistic(N0,N1,snps,W,gamma.sim,freq,GenoProbList)
	return(exp_z_score)
}

##' @title Compute a simulated Z Score
##' @export
##' @param N0 The number of Y=0
##' @param N1	The number of Y=1
##' @param snps The snps at which we wish to compute the expected Z Score
##' @param W	The true causal SNPs (these need not be in "snps")
##' @param gamma.CV	The odds ratios of effect of the true causal SNPs
##' @param freq Frequencies of SNP appearances (computed using snphap)
##' @param df_control A reference set of control samples
##' @author Mary Fortune
simulated_z_score<-function(N0,N1,snps,W,gamma.CV,freq,df_control){
	GenoProbList<-make_GenoProbList(snps,W,freq)
	gamma.sim<-c(compute_gamma0(N0,N1,W,gamma.CV,freq),gamma.CV)
	exp_z_score<-est_statistic(N0,N1,snps,W,gamma.sim,freq,GenoProbList)
	XX<-new("SnpMatrix", as.matrix(df_control))
	LD <- snpStats::ld(XX,XX,stat="R",symmetric=TRUE)
	LD<-as.matrix(make.positive.definite(LD))
	sim_z_score<-c(rmvnorm(n=1,mean=exp_z_score,sigma=LD))
	return(sim_z_score)
}
