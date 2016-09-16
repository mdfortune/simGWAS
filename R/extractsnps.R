library(dplyr)


#computes a matrix telling us which haplotype pairs correspond to which genotypes
which_genotypes<-function(nsnps){
	if (nsnps<2){return(0)}
	npheno<-2^nsnps
	#tells you which haplotype pairs correspond to which SNPs, for number of SNPS = n
	haptype<-(hcube(rep(2,nsnps))-1)
	Genotype<-apply(expand.grid(1:npheno,1:npheno),1,function(i) (colSums(haptype[i,])))
	Genotype<-apply(Genotype,2,function(x){paste(x,collapse="")})
	#We order the genotypes by the order we analyse them in the est_zscore code
	unique_Genotype<-apply(hcube(rep(3,nsnps))-1, 1,function(x){paste(x,collapse="")})
	geno_matrix<-matrix(0,length(unique_Genotype),npheno^2)
	for (ii in 1:length(unique_Genotype)){
		geno_matrix[ii,which(Genotype==unique_Genotype[ii])]<-1
	}
	rownames(geno_matrix)<-unique_Genotype
	return(geno_matrix)
}


geno_matrix_2snps<-which_genotypes(2)
geno_matrix_3snps<-which_genotypes(3)
geno_matrix_4snps<-which_genotypes(4)
geno_matrix_5snps<-which_genotypes(5)
geno_matrix_6snps<-which_genotypes(6)

#To know which outputs to give for X=0, X=1,X=2
geno_2SNP<-hcube(rep(3,2))-1
geno_3SNP<-hcube(rep(3,3))-1
geno_4SNP<-hcube(rep(3,4))-1
geno_5SNP<-hcube(rep(3,5))-1
geno_6SNP<-hcube(rep(3,6))-1

which_X0_2SNP<-which(geno_2SNP[,1]==0)
which_X1_2SNP<-which(geno_2SNP[,1]==1)
which_X2_2SNP<-which(geno_2SNP[,1]==2)

which_X0_3SNP<-which(geno_3SNP[,1]==0)
which_X1_3SNP<-which(geno_3SNP[,1]==1)
which_X2_3SNP<-which(geno_3SNP[,1]==2)

which_X0_4SNP<-which(geno_4SNP[,1]==0)
which_X1_4SNP<-which(geno_4SNP[,1]==1)
which_X2_4SNP<-which(geno_4SNP[,1]==2)

which_X0_5SNP<-which(geno_5SNP[,1]==0)
which_X1_5SNP<-which(geno_5SNP[,1]==1)
which_X2_5SNP<-which(geno_5SNP[,1]==2)

which_X0_6SNP<-which(geno_6SNP[,1]==0)
which_X1_6SNP<-which(geno_6SNP[,1]==1)
which_X2_6SNP<-which(geno_6SNP[,1]==2)




extractsnps<-function(X,W,freq){
	if (X %in% W){
		snp.int<-W
	}else{
		snp.int<-c(X,W)
	}	
	nsnps<-length(snp.int)
	#the haplotypes for our snps of interest
	haptype<-summarise(group_by_(freq,.dots=snp.int),totalProb=sum(Probability))
	#get the genotype matrix
	if (nsnps==1){
		return(list(c(haptype$totalProb[1]^2,0,0),c(0,2*haptype$totalProb[1]*haptype$totalProb[2],0),c(0,0,haptype$totalProb[2]^2)))
	}
	else if (nsnps==2){
		genotype<-extractsnps_2snps(haptype)
		if (X %in% W){
			whichX<-which(W==X)
			return(list(genotype*as.integer(geno_2SNP[,whichX]==0),genotype*as.integer(geno_2SNP[,whichX]==1),genotype*as.integer(geno_2SNP[,whichX]==2)))
		}else{
			return(list(genotype[which_X0_2SNP],genotype[which_X1_2SNP],genotype[which_X2_2SNP]))
		}
	}else if (nsnps==3){	
		genotype<-extractsnps_3snps(haptype)
		if (X %in% W){
			whichX<-which(W==X)
			return(list(genotype*as.integer(geno_3SNP[,whichX]==0),genotype*as.integer(geno_3SNP[,whichX]==1),genotype*as.integer(geno_3SNP[,whichX]==2)))
		}else{
			return(list(genotype[which_X0_3SNP],genotype[which_X1_3SNP],genotype[which_X2_3SNP]))
		}
	}else if (nsnps==4){	
		genotype<-extractsnps_4snps(haptype)
		if (X %in% W){
			whichX<-which(W==X)
			return(list(genotype*as.integer(geno_4SNP[,whichX]==0),genotype*as.integer(geno_4SNP[,whichX]==1),genotype*as.integer(geno_4SNP[,whichX]==2)))
		}else{
			return(list(genotype[which_X0_4SNP],genotype[which_X1_4SNP],genotype[which_X2_4SNP]))
		}
	}else if (nsnps==5){	
		genotype<-extractsnps_5snps(haptype)
		if (X %in% W){
			whichX<-which(W==X)
			return(list(genotype*as.integer(geno_5SNP[,whichX]==0),genotype*as.integer(geno_5SNP[,whichX]==1),genotype*as.integer(geno_5SNP[,whichX]==2)))
		}else{
			return(list(genotype[which_X0_5SNP],genotype[which_X1_5SNP],genotype[which_X2_5SNP]))
		}
	}else if (nsnps==6){	
		genotype<-extractsnps_6snps(haptype)
		if (X %in% W){
			whichX<-which(W==X)
			return(list(genotype*as.integer(geno_6SNP[,whichX]==0),genotype*as.integer(geno_6SNP[,whichX]==1),genotype*as.integer(geno_6SNP[,whichX]==2)))
		}else{
			return(list(genotype[which_X0_6SNP],genotype[which_X1_6SNP],genotype[which_X2_6SNP]))
		}
	}else{
		stop("Causal model has too many SNPs")
	}
	#output the relevant portions of the probability column
	
}

extractsnps_2snps<-function(haptype){
	hap_Probs<-haptype$totalProb
	names(hap_Probs)<-apply((haptype[,-3]-1),1,paste,collapse="")
	#hap_Probs<-hap_Probs[apply(hcube(rep(2,2))-1, 1,function(x){paste(x,collapse="")})]
	hap_Probs<-hap_Probs[c( "00", "10", "01", "11")]
	hap_Probs[is.na(hap_Probs)]<-0
	#Probs<-c(hap_Probs%*%t(hap_Probs))
	#genotype<-c(geno_matrix_2snps%*%c(hap_Probs%*%t(hap_Probs)))
	genotype<-MatrixVector(geno_matrix_2snps,c(hap_Probs%*%t(hap_Probs)),verbose=F)
	#return(cbind(hcube(rep(3,2))-1,genotype))
	return(genotype)
}


extractsnps_3snps<-function(haptype){
	hap_Probs<-haptype$totalProb
	names(hap_Probs)<-apply((haptype[,-4]-1),1,paste,collapse="")
	#add extra columns with probability 0 for the haplotypes which do not appear in our dataset
	#hap_Probs<-hap_Probs[apply(hcube(rep(2,3))-1, 1,function(x){paste(x,collapse="")})]
	hap_Probs<-hap_Probs[c( "000", "100", "010", "110", "001", "101", "011", "111")]
	hap_Probs[is.na(hap_Probs)]<-0
	#compute genotype probabilities (may be some duplicates)
	#Probs<-c(hap_Probs%*%t(hap_Probs))
	#sum together hap pairs corresponding to the same genotype
	#genotype<-c(geno_matrix_3snps%*%c(hap_Probs%*%t(hap_Probs)))
	genotype<-MatrixVector(geno_matrix_3snps,c(hap_Probs%*%t(hap_Probs)),verbose=F)
	#return(cbind(hcube(rep(3,3))-1,genotype))
	return(genotype)
}

extractsnps_4snps<-function(haptype){
	hap_Probs<-haptype$totalProb
	names(hap_Probs)<-apply((haptype[,-5]-1),1,paste,collapse="")
	#add extra columns with probability 0 for the haplotypes which do not appear in our dataset
	hap_Probs<-hap_Probs[apply(hcube(rep(2,4))-1, 1,function(x){paste(x,collapse="")})]
	hap_Probs[is.na(hap_Probs)]<-0
	#compute genotype probabilities (may be some duplicates)
	#Probs<-c(hap_Probs%*%t(hap_Probs))
	#sum together hap pairs corresponding to the same genotype
	#genotype<-c(geno_matrix_4snps%*%c(hap_Probs%*%t(hap_Probs)))
	genotype<-MatrixVector(geno_matrix_4snps,c(hap_Probs%*%t(hap_Probs)),verbose=F)
	return(genotype)
}

extractsnps_5snps<-function(haptype){
	hap_Probs<-haptype$totalProb
	names(hap_Probs)<-apply((haptype[,-6]-1),1,paste,collapse="")
	#add extra columns with probability 0 for the haplotypes which do not appear in our dataset
	hap_Probs<-hap_Probs[apply(hcube(rep(2,5))-1, 1,function(x){paste(x,collapse="")})]
	hap_Probs[is.na(hap_Probs)]<-0
	#compute genotype probabilities (may be some duplicates)
	#Probs<-c(hap_Probs%*%t(hap_Probs))
	#sum together hap pairs corresponding to the same genotype
	#genotype<-c(geno_matrix_5snps%*%c(hap_Probs%*%t(hap_Probs)))
	genotype<-MatrixVector(geno_matrix_5snps,c(hap_Probs%*%t(hap_Probs)),verbose=F)
	return(genotype)
}

extractsnps_6snps<-function(haptype){
	hap_Probs<-haptype$totalProb
	names(hap_Probs)<-apply((haptype[,-7]-1),1,paste,collapse="")
	#add extra columns with probability 0 for the haplotypes which do not appear in our dataset
	hap_Probs<-hap_Probs[apply(hcube(rep(2,6))-1, 1,function(x){paste(x,collapse="")})]
	hap_Probs[is.na(hap_Probs)]<-0
	#compute genotype probabilities (may be some duplicates)
	#Probs<-c(hap_Probs%*%t(hap_Probs))
	#sum together hap pairs corresponding to the same genotype
	#genotype<-c(geno_matrix_6snps%*%c(hap_Probs%*%t(hap_Probs)))
	genotype<-MatrixVector(geno_matrix_6snps,c(hap_Probs%*%t(hap_Probs)),verbose=F)
	return(genotype)
}




