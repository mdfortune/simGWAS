#include <Rcpp.h>
using namespace Rcpp;
#include <emmintrin.h>
#include <ctime>

#define ROUND_DOWN(x, s) ((x) & ~((s)-1))

// [[Rcpp::export]]
NumericVector MatrixVector(const NumericMatrix& matrix, const NumericVector& vector, const bool verbose = false){
    const int nrow = matrix.nrow();
    const int ncol = matrix.ncol();
    NumericVector ret(nrow);
    const double* ref;
    
    std::clock_t start;
    start = std::clock();
    
    for(int i = 0; i < ncol; ++i){
        ref = &matrix[nrow*i];
        // Load vector[i]
        __m128d P = _mm_set1_pd(vector[i]);
        
        int j = 0;
         for(; j < ROUND_DOWN(nrow, 4); j+=4){
            __m128d v0 = _mm_loadu_pd(&ref[j]); // Stupid R limitation
            __m128d v1 = _mm_loadu_pd(&ref[j+2]); // Stupid R limitation
            __m128d r = _mm_loadu_pd(&ret[j]); // Stupid R limitation
             __m128d r1 = _mm_loadu_pd(&ret[j+2]); // Stupid R limitation
            v0 = _mm_mul_pd(v0, P);
            v1 = _mm_mul_pd(v1, P);
            v0 = _mm_add_pd(v0, r);
            v1 = _mm_add_pd(v1, r1);
            _mm_storeu_pd(&ret[j], v0); // Stupid R limitation
             _mm_storeu_pd(&ret[j+2], v1); // Stupid R limitation
        }
        for(; j < nrow; ++j){
            ret[j] += ref[j]*vector[i];
        }
    }
    
    if(verbose){
        Rcout << "Clocks/cycle: " << CLOCKS_PER_SEC << std::endl;
        const double timediff = (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000);
        Rcout << "Operations: " << ncol*nrow*2 << " --> " << (ncol*nrow*2)/timediff/1000 << " MFLOPS" << std::endl;
         Rcpp::Rcout << "Time: " << timediff << " ms" << std::endl;
    }
    return(ret);
}
