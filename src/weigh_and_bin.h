#ifndef MARKSUMMARY_WEIGH_AND_BIN_H
#define MARKSUMMARY_WEIGH_AND_BIN_H

#include <Rcpp.h>

RcppExport SEXP weighAndBin(SEXP sFValue, SEXP sWeight, SEXP sIdxVec,
                            SEXP sNBin);

#endif
