#a function which takes in N0,N1, freq and the other values of gamma, and returns the value of gamma_0 to use
#note that gamma.CV is gamma_1...gamma_m only
#W should be the names of the SNPs, not the locations of the SNPs

compute_gamma0<-function(N0,N1,W,gamma.CV,freq){
  eta.CV<-rowSums(sweep((hcube(rep(3,length(W)))-1),MARGIN=2,gamma.CV,`*`))
  e.eta.CV<-exp(eta.CV)
  GenoProbW<-make_GenoProbList(W[1],W,freq)[[1]]
  PW<-GenoProbW[[1]]+GenoProbW[[2]]+GenoProbW[[3]]
  denom<-e.eta.CV*PW
  gamma0<-log(N1)-log(N0)-log(sum(denom))
  return(gamma0)
}


#for when we know PW
compute_gamma0_PW<-function(N0,N1,W,gamma.CV,PWgY0){
  e.eta.CV<-exp(rowSums(sweep((hcube(rep(3,length(W)))-1),MARGIN=2,gamma.CV,`*`)))
  gamma0<-log(N1)-log(N0)-log(sum(e.eta.CV*PWgY0))
  return(gamma0)
}
