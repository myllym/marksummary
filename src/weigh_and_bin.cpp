#include "weigh_and_bin.h"

SEXP weighAndBin(SEXP sFValue, SEXP sWeight, SEXP sIdxVec, SEXP sNBin) {
    BEGIN_RCPP

    Rcpp::NumericVector fValue(sFValue);
    Rcpp::NumericVector weight(sWeight);
    Rcpp::IntegerVector idxVec(sIdxVec);
    // Rcpp uses int instead of size_t. Limitation of R.
    int nBin = Rcpp::as<int>(sNBin);
    int nPair = fValue.size();

    Rcpp::NumericVector result(nBin);

    for (int i = 0; i < nPair; ++i) {
        result[idxVec[i]] += fValue[i] * weight[i];
    }

    return result;

    END_RCPP
}
