library(devtools)
library(annotSnpStats)
load_all("~/RP/simGWAS")

(load("~/scratch/test_data_2CV.RData"))

load_all("~/RP/simGWAS")

e.old <- extractsnps(snps[1],W,freq)
e.new <- fastextractsnps(snps[1],W,freq)
all.equal(do.call("cbind",e.old), do.call("cbind",e.new))

haplabs(4)

library(microbenchmark)
## microbenchmark(extractsnps(snps[6],W[1],freq),
##                fastextractsnps(snps[6],W[1],freq))

old <- make_GenoProbList(snps[1:5],W[1],freq)
new <- fastmake_GenoProbList(snps[1:5],W[1],freq)
all.equal(old,new)

microbenchmark(make_GenoProbList(snps[1:5],W[1],freq),
               fastmake_GenoProbList(snps[1:5],W[1],freq))
microbenchmark(make_GenoProbList(snps[1:5],W,freq),
               fastmake_GenoProbList(snps[1:5],W,freq))

GP <- fastmake_GenoProbList(snps,W,freq)
load_all("~/RP/simGWAS")
Z <- est_statistic(N0,N1,snps,W,rep(1.5,length(W)),freq,GP)
Zf <- fast_statistic(N0,N1,snps,W,rep(1.5,length(W)),freq,GP)
all.equal(Z,Zf)

microbenchmark(Z <- est_statistic(N0,N1,snps,W,rep(1.5,length(W)),freq,GP),
               Zf <- fast_statistic(N0,N1,snps,W,rep(1.5,length(W)),freq,GP))


names(Z) <- snps
Z[W]

GPold <- 	lapply(snps, extractsnps, W=W,freq=freq)
GPnew <- 	lapply(snps, fastextractsnps, W=W,freq=freq)
identical(GPold,GPnew)
identical(GenoProbList,GPnew)

simGWAS:::make_GenoProbList(snps=snps,W=W,freq=freq)

library(microbenchmark)
microbenchmark(happrobs(G,freq$Probability),
               extractsnps(X,W,freq))

################################################################################
## simulations
library(combinat)
library(devtools)
load_all("~/RP/simGWAS")
nsnps <- 100
nhaps <- 1000
maf <- runif(nsnps+1,0.1,0.5)
haps <- do.call("cbind", lapply(maf, function(f) rbinom(nhaps,1,f)))
r <- cor(haps)
summary(r[upper.tri(r)])
## add some autocorrelation
haps <- pmin(haps[,1:nsnps]+haps[,-1],1)
r <- cor(haps)
summary(r[upper.tri(r)])

snps <- colnames(haps) <- paste0("s",1:nsnps)
freq <- as.data.frame(haps+1)
freq$Probability <- 1/nrow(freq)
sum(freq$Probability)
N0 <- 2000
N1 <- 3000
N <- N1+N0
g1 <- 1.4

load_all("~/RP/simGWAS")
wrapper <- function() {
    g0 <- compute_gamma0(N0=N0,N1=N1,W=snps[CV],gamma.CV=g1,freq=freq)
    GP <- make_GenoProbList(snps=snps,W=snps[CV],freq=freq)
    est_statistic(N0,N1,snps,W=snps[CV],gamma=c(g0,g1),freq,GP)
}
fwrapper <- function() {
    FP <- fastmake_GenoProbList(snps=snps,W=snps[CV],freq=freq)
    fast_statistic(N0,N1,snps,W=snps[CV],gamma1=g1,freq,FP) ## NB gamma1 != gamma above
}

zold <- wrapper()
znew <- fwrapper()
all.equal(zold,znew)
plot(1:nsnps,zold); abline(v=CV,col="red"); points(1:nsnps,znew,col="blue"); abline(h=0)

microbenchmark::microbenchmark(wrapper(), fwrapper(), times=10)



## H <- haps[sample(1:nhaps,N0,replace=TRUE),] + haps[sample(1:nhaps,N0,replace=TRUE),]
## summary(colMeans(H)/2 - maf)

## PWgY0 <- p0 <- c(0.25,0.5,0.25)
## PWgY1 <- p1 <- c(0.16, 0.48, 0.36)
## CV <- which.min(abs(maf-0.5))
## df <- cbind(Y=0,H)
## D <- make_dataset(df,N0=N0,N1=N1,CV=CV,PWgY0=p0,PWgY1=p1)
## CV
## (g0 <- compute_gamma0(N1,N1,snps[CV],1.5,freq))
## g <- c(g0,1.4)

GPf <- fastmake_GenoProbList(snps=snps,W=snps[CV],freq=freq)
est_zscore(N0,N1,Ufactor,powerfactor,freq,GP[[1]])
est_zscore(N0,N1,Ufactor,powerfactor,freq,GP[[CV]])
Z <- make_zscore(N0,N1,snps[CV],1.5,freq,GP)
head(Z)


(g0 <- compute_gamma0(N1,N1,snps[1],1.5,freq))
gamma <- c(g0,1.5)
expeta<-exp(gamma[1]+rowSums(sweep((hcube(rep(3,length(gamma)-1))-1),MARGIN=2,gamma[-1],"*")))
#compute the constant factors we will multiply by
Ufactor<-N0*(N-1)*(N0*expeta-N1)/(N^2)
powerfactor<-N0*(expeta+1)/N

GP <- simGWAS:::make_GenoProbList(snps=snps,W=snps[CV],freq=freq)
est_zscore(N0,N1,Ufactor,powerfactor,freq,GP[[1]])
est_zscore(N0,N1,Ufactor,powerfactor,freq,GP[[CV]])
Z <- make_zscore(N0,N1,snps[CV],1.5,freq,GP)
head(Z)

library(Rcpp)
sourceCpp("~/RP/simGWAS/src/mean.cpp")


library(microbenchmark)
x <- runif(1e5)
microbenchmark(
  mean(x),
  meanC(x)
)

X=snps[1]
extractsnps(X,W,freq)

microbenchmark(extractsnps(snps[1],W,freq))



?strtoi

## THIS IS FASTER
x1=GenoProbList[[1]]
microbenchmark(
zold=sapply(1:nsnps, function(i)
    est_zscore(N0=N0,N1=N1,Ufactor=Ufactor,powerfactor=powerfactor,freq=freq,GenoProbXW=x1)),
znew=sapply(1:nsnps, function(i)
    fast_zscore(N0=N0,N1=N1,Ufactor=Ufactor,powerfactor=powerfactor,freq=freq,GenoProbXW=x1)))
all.equal(zold,znew)

## THIS IS NOT
nsnps <- length(GenoProbList)
microbenchmark(
zold <- sapply(1:nsnps, function(i)
    est_zscore(N0=N0,N1=N1,Ufactor=Ufactor,powerfactor=powerfactor,
               freq=freq,GenoProbXW=GenoProbList[[i]])),
znew <- sapply(1:nsnps, function(i)
    fast_zscore(N0=N0,N1=N1,Ufactor=Ufactor,powerfactor=powerfactor,
                freq=freq,GenoProbXW=GenoProbList[[i]])))
zold <- sapply(1:nsnps, function(i)
    est_zscore(N0=N0,N1=N1,Ufactor=Ufactor,powerfactor=powerfactor,
               freq=freq,GenoProbXW=GenoProbList[[i]]))
znew <- sapply(1:nsnps, function(i)
    fast_zscore(N0=N0,N1=N1,Ufactor=Ufactor,powerfactor=powerfactor,
                freq=freq,GenoProbXW=GenoProbList[[i]]))
all.equal(zold,znew)

## THIS IS NOT
microbenchmark(
zold <- sapply(1:nsnps, function(i) est_zscore(N0=N0,N1=N1,Ufactor=Ufactor,powerfactor=powerfactor,freq=freq,GenoProbXW=GenoProbList[[2]])),
znew <- sapply(1:nsnps, function(i) fast_zscore(N0=N0,N1=N1,Ufactor=Ufactor,powerfactor=powerfactor,freq=freq,GenoProbXW=GenoProbList[[2]])))
all.equal(zold,znew)

