#include <Rcpp.h>
using namespace Rcpp;

/* Compute the sum of squares between x and y, weighted by w */

// [[Rcpp::export]]

double wsumsq( const NumericVector& xx,  const NumericVector& yy,  const NumericVector& ww){
  	int k;
  	double sumw =0;
	double sumxy =0;
	double sumsq =0;
  	/*  get dimensions of data */
	const int n = xx.size();
        /* for loop to calculate sums, the for loop is needed as we have a check for pairwise missing: ISNAN  */
      for (k = 0; k < n; k++){
        if (!ISNAN(xx[k]) && !ISNAN(yy[k])){
          sumw += ww[k];
          sumxy += ww[k]*(xx[k]-yy[k])*(xx[k]-yy[k]);
        }
      }
      /*  final correlation */
      sumsq = sumxy/sumw;
  return sumsq;
}
