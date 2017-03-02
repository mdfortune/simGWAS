#include<Rcpp.h>
using namespace Rcpp;

/* Helper for est_zscore.R */

// [[Rcpp::export]]
double zscore( const double N0, const double N1, const NumericVector& Ufactor,  const NumericVector& powerfactor,
	       const NumericVector& PX1W, const NumericVector& PX2W){
  double N=N0+N1;
  double UfactorPX2W = sum(Ufactor*PX2W);
  double UfactorPX1W = sum(Ufactor*PX1W);
  double powerfactorPX2W = sum(powerfactor*PX2W);
  double powerfactorPX1W = sum(powerfactor*PX1W);
  double U = 2*UfactorPX2W+UfactorPX1W;
  // compute E(X^a) for various powers a
  double EX = 2*powerfactorPX2W+powerfactorPX1W;
  double EX2 = 4*powerfactorPX2W+powerfactorPX1W;
  double EX3 = 8*powerfactorPX2W+powerfactorPX1W;
  double EX4 = 16*powerfactorPX2W+powerfactorPX1W;
  //find expected value of VX
  double EVX = EX2-pow(EX,2.0);
  //find expected value of VX^2
  double Term1 = N*EX4+N*(N-1)*pow(EX2,2.0);
  double Term2 = EX4+2*(N-1)*EX3*EX+(N-1)*pow(EX2,2.0)+(N-1)*(N-2)*EX2*pow(EX,2.0);
  double Term3 = (EX4+4*(N-1)*EX3*EX+6*(N-1)*pow(EX2,2.0)+6*(N-1)*(N-2)*EX2*pow(EX,2.0)+(N-1)*(N-2)*(N-3)*pow(EX,4.0))/N;
  double EVX2 = pow(1/(N-1),2.0) * ( Term1 - (2*Term2) + (Term3));
  // distribution of V
  double a = (2*EVX2-EVX*EVX)/(EVX2-EVX*EVX);
  double b = (EVX2*EVX)/(EVX2-EVX*EVX);
  double EinvsqrtVX = pow(b,-0.5)*exp(lgamma((2*a+1)/2)-lgamma(a));
  //Value of Z
  double Z = U*EinvsqrtVX*pow(N/(N0*N1),0.5);
  return Z;
  // return sum(x*y);
  //   const int n = x.size();
  // for(int i=0; i<n; i++) {
  //   //    if (!ISNAN(x[i]) && !ISNAN(y[i]))
  //     ret =+ x[i]*y[i];
  // }
  // return ret;
}

// [[Rcpp::export]]
double psum( const NumericVector& x,  const NumericVector& y){
  double ret=0;
  return sum(x*y);
  //   const int n = x.size();
  // for(int i=0; i<n; i++) {
  //   //    if (!ISNAN(x[i]) && !ISNAN(y[i]))
  //     ret =+ x[i]*y[i];
  // }
  // return ret;
}
