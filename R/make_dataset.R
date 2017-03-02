##' @title Make simulated data which has N0 cases, N1 controls and is consitent with the causal model {W,gamma}
##' @param df dataset of genotypes to sample from (the control samples)
##' @param N0 number of samples with Y=0
##' @param N1 number of samples with Y=1
##' @param CV which SNPs are causal?
##' @param PWgY0 the distribution of causal snps in controls
##' @param PWgY1 the distribution of causal snps in cases
##' @return simulated genotype data for N0 controls and N1 cases
##' @author Mary Fortune
make_dataset<-function(df,N0,N1,CV,PWgY0,PWgY1){
    m <- length(CV)
	#make sure we have genotypes in {0,1,2}
	if (isTRUE(all.equal(sort(unique(df[,2])),1:3))){
		df[,-1]<-df[,-1]-1
	}
	#compute the weights for sampling the cases
	ref_data_W<-df[,CV+1]
	Y1weights<-rep(0,nrow(df))
	#which order do the Ws come in the vectors of probabilities above?
	orderW<-apply(hcube(rep(3,m))-1, 1,function(x){paste(x,collapse="")})
	for (ii in 1:nrow(df)){
		if (length(CV)==1){
			Wvalue<-paste(ref_data_W[ii])
		}else{
			Wvalue<-paste(ref_data_W[ii,],collapse="")
		}
		whichw<-which(orderW==Wvalue)
		Y1weights[ii]<-PWgY1[whichw]/PWgY0[whichw]
	}	
	rows.index0<-sample(nrow(df), size=N0, replace = TRUE)
	rows.index1<-sample(nrow(df), size=N1, replace = TRUE, prob=Y1weights)
	dfnew_Y0<-cbind(rep(0,N0),df[rows.index0,-1])
	dfnew_Y1<-cbind(rep(1,N1),df[rows.index1,-1])
	colnames(dfnew_Y0)[1]<-"Y"
	colnames(dfnew_Y1)[1]<-"Y"
	dfnew<-rbind(dfnew_Y0,dfnew_Y1)
	return(list(dfnew,c(rows.index0,rows.index1)))
}
