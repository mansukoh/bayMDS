// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// distRcpp
Rcpp::NumericMatrix distRcpp(Rcpp::NumericMatrix X);
RcppExport SEXP _bayMDS_distRcpp(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(distRcpp(X));
    return rcpp_result_gen;
END_RCPP
}
// bmdsMCMC
extern "C" SEXP bmdsMCMC(SEXP DIST, SEXP p, int nwarm, int niter);
RcppExport SEXP _bayMDS_bmdsMCMC(SEXP DISTSEXP, SEXP pSEXP, SEXP nwarmSEXP, SEXP niterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type DIST(DISTSEXP);
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type nwarm(nwarmSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    rcpp_result_gen = Rcpp::wrap(bmdsMCMC(DIST, p, nwarm, niter));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bayMDS_distRcpp", (DL_FUNC) &_bayMDS_distRcpp, 1},
    {"_bayMDS_bmdsMCMC", (DL_FUNC) &_bayMDS_bmdsMCMC, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_bayMDS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
