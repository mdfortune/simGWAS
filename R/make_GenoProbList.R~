Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

Rcpp::sourceCpp('/home/mdf34/simgwas_package/test_prop.cpp')
Rcpp::sourceCpp('/home/mdf34/simgwas_package/algebra.cpp')
vals_1SNP<-expand.grid(1:2)
vals_2SNP<-expand.grid(1:2,1:2)
vals_3SNP<-expand.grid(1:2,1:2,1:2)


#To know which outputs to give for X=0, X=1,X=2
geno_2SNP<-hcube(rep(3,2))-1
geno_3SNP<-hcube(rep(3,3))-1

which_X0_2SNP<-which(geno_2SNP[,1]==0)
which_X1_2SNP<-which(geno_2SNP[,1]==1)
which_X2_2SNP<-which(geno_2SNP[,1]==2)
which_X0_3SNP<-which(geno_3SNP[,1]==0)
which_X1_3SNP<-which(geno_3SNP[,1]==1)
which_X2_3SNP<-which(geno_3SNP[,1]==2)


make_GenoProbList<-function(snps,W,freq){
	nsnps<-length(snps)
	m<-length(W)
	GenoProbList<-vector("list", nsnps) 
	for (ii in 1:nsnps){
		X<-snps[ii]
		GenoProbList[[ii]]<-extractsnps(X,W,freq)
	}	
	return(GenoProbList)
}


make_GenoProbList_Wsize1<-function(snps,W,freq){
	# A specific function to make the GenoProbList when W has size 1
	which.snps<-which(colnames(freq) %in% snps)
	which.W<-which(colnames(freq)==W)
	nsnps<-length(snps)
	GenoProbList<-vector("list", nsnps) 
	#compute the haplotypes
	x<-rbind(which.snps,rep(which.W[1],nsnps))
	vals_2SNP<-expand.grid(1:2,1:2)
	hapProb<-combinationRefs(as.matrix(freq[,-ncol(freq)]),x,as.matrix(vals_2SNP),c(freq[,ncol(freq)]))
	#compute the genotypes
	genoProb<-t(apply(hapProb,1,function(hap_Probs){MatrixVector(geno_matrix_2snps,c(hap_Probs%*%t(hap_Probs)),verbose=F)}))
	for (ii in 1:nsnps){
		if (snps[ii]==W){	
			genotype<-genoProb[ii,c(1,5,9)]
			GenoProbList[[ii]]<-list(c(genotype[1],0,0),c(0,genotype[2],0),c(0,0,genotype[3]))
		}
		else{
			genotype<-genoProb[ii,]
			GenoProbList[[ii]]<-list(genotype[which_X0_2SNP],genotype[which_X1_2SNP],genotype[which_X2_2SNP])
		}
	}	
	return(GenoProbList)
}

make_GenoProbList_Wsize2<-function(snps,W,freq){
	# A specific function to make the GenoProbList when W has size 2
	which.snps<-which(colnames(freq) %in% snps)
	nsnps<-length(snps)
	GenoProbList<-vector("list", nsnps) 
	#compute the haplotypes
	x<-rbind(which.snps,rep(which(colnames(freq)==W[1]),nsnps),rep(which(colnames(freq)==W[2]),nsnps))
	vals_3SNP<-expand.grid(1:2,1:2,1:2)
	hapProb<-combinationRefs(as.matrix(freq[,-ncol(freq)]),x,as.matrix(vals_3SNP),c(freq[,ncol(freq)]))
	#compute the genotypes
	genoProb<-t(apply(hapProb,1,function(hap_Probs){MatrixVector(geno_matrix_3snps,c(hap_Probs%*%t(hap_Probs)),verbose=F)}))
	for (ii in 1:nsnps){
		if (snps[ii]==W[1]){	
			whichX=1
			genotype<-genoProb[ii,c(1,5,9,10,14,18,19,23,27)]
			GenoProbList[[ii]]<-list(genotype*as.integer(geno_2SNP[,whichX]==0),genotype*as.integer(geno_2SNP[,whichX]==1),genotype*as.integer(geno_2SNP[,whichX]==2))
		}else if (snps[ii]==W[2]){
			whichX=2
			genotype<-genoProb[ii,c(1,4,7,11,14,17,21,24,27)]
			GenoProbList[[ii]]<-list(genotype*as.integer(geno_2SNP[,whichX]==0),genotype*as.integer(geno_2SNP[,whichX]==1),genotype*as.integer(geno_2SNP[,whichX]==2))
		}
		else{
			genotype<-genoProb[ii,]
			GenoProbList[[ii]]<-list(genotype[which_X0_3SNP],genotype[which_X1_3SNP],genotype[which_X2_3SNP])
		}
	}	
	return(GenoProbList)
}

