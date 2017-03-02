#include<Rcpp.h>
#include<cstring>
#include<cstdlib>
#include<bitset>
#include <iostream> // for debugging
using namespace Rcpp;

// [[Rcpp::export]]
std::vector< std::string > haplabs(const int n) {
  const int n2=pow(2,n);
  std::vector< std::string > labs;
  for(int i=0; i<n2; i++) {
    std::bitset< sizeof(i)*CHAR_BIT > bits( i );
    std::string bitstr = bits.to_string();
    int j = bitstr.length();
    std::string tmp("0000000000",n);
    // std::cout << tmp << " -> ";
    for (int k=0; k<n; k++) {
      tmp[k] = bitstr[--j];
    }    
    // std::cout << tmp << std::endl;
    labs.push_back(tmp);
  }
  return labs;
}

// [[Rcpp::export]]
NumericVector happrobs(const NumericMatrix& G, const NumericVector& P) {
  const int m=G.ncol();
  const int n=G.nrow();
  const int m2=pow(2,m);
  // std::cout << "nsnps: " << m << std::endl;
  std::vector< std::string > labs = haplabs(m);
  NumericVector pr(m2);
  std::map<std::string,double> mpr;
  for(int i=0; i<m2; i++)
    mpr[ labs[i] ] = 0;
  for(int i=0; i<n; i++) {
    // label
    std::string line("000000000",m);
    for(int j=0; j<m; j++) {
      if( G(i,j) != 0)
	line[j]='1';
    }
    mpr[line] += P[i];
  }
  // sort into named NumericVector for return
  for(int i=0; i<m2; i++)
    pr[i] = mpr[ labs[i] ];
  pr.attr("names") = labs;
  return pr;
}
