#include "sum_vector.h"

SEXP sum_vector(SEXP sValue, SEXP sIdxVec, SEXP sNBin) {
    BEGIN_RCPP

    Rcpp::NumericVector value(sValue);
    Rcpp::IntegerVector idxVec(sIdxVec);
    // Rcpp uses int instead of size_t. Limitation of R.
    int nBin = Rcpp::as<int>(sNBin);
    int nPair = value.size();

    Rcpp::NumericVector result(nBin);

    //idxVec = idxVec - 1L;

    for (int i = 0; i < nPair; ++i) {
        result[idxVec[i]] += value[i];
    }
    return result;

    END_RCPP
}
