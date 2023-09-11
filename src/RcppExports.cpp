// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// glsCpp
arma::mat glsCpp(arma::mat Y, arma::mat Z, arma::mat Z2, arma::mat H, arma::mat h, arma::mat A0, arma::mat B0, arma::mat C0, arma::cube Omega0, arma::vec break_pts, int r, int n_iter);
RcppExport SEXP _narwhalCA_glsCpp(SEXP YSEXP, SEXP ZSEXP, SEXP Z2SEXP, SEXP HSEXP, SEXP hSEXP, SEXP A0SEXP, SEXP B0SEXP, SEXP C0SEXP, SEXP Omega0SEXP, SEXP break_ptsSEXP, SEXP rSEXP, SEXP n_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z2(Z2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type H(HSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A0(A0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type C0(C0SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Omega0(Omega0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type break_pts(break_ptsSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type n_iter(n_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(glsCpp(Y, Z, Z2, H, h, A0, B0, C0, Omega0, break_pts, r, n_iter));
    return rcpp_result_gen;
END_RCPP
}
// glsConstCpp
arma::mat glsConstCpp(arma::mat Y, arma::mat Z, arma::mat Z2, arma::mat H, arma::mat h, arma::mat A0, arma::mat B0, arma::mat C0, arma::mat Omega0, arma::vec break_pts, int r, int n_iter);
RcppExport SEXP _narwhalCA_glsConstCpp(SEXP YSEXP, SEXP ZSEXP, SEXP Z2SEXP, SEXP HSEXP, SEXP hSEXP, SEXP A0SEXP, SEXP B0SEXP, SEXP C0SEXP, SEXP Omega0SEXP, SEXP break_ptsSEXP, SEXP rSEXP, SEXP n_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z2(Z2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type H(HSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A0(A0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type C0(C0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Omega0(Omega0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type break_pts(break_ptsSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type n_iter(n_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(glsConstCpp(Y, Z, Z2, H, h, A0, B0, C0, Omega0, break_pts, r, n_iter));
    return rcpp_result_gen;
END_RCPP
}
// vecm
arma::mat vecm(arma::mat Z0, arma::mat Z1, arma::mat Z2, int r, double dt, bool intercept);
RcppExport SEXP _narwhalCA_vecm(SEXP Z0SEXP, SEXP Z1SEXP, SEXP Z2SEXP, SEXP rSEXP, SEXP dtSEXP, SEXP interceptSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Z0(Z0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z1(Z1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z2(Z2SEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< bool >::type intercept(interceptSEXP);
    rcpp_result_gen = Rcpp::wrap(vecm(Z0, Z1, Z2, r, dt, intercept));
    return rcpp_result_gen;
END_RCPP
}
// johansenCpp
arma::mat johansenCpp(arma::mat Z0, arma::mat Z1, arma::mat Z2, int r, double dt, bool intercept);
RcppExport SEXP _narwhalCA_johansenCpp(SEXP Z0SEXP, SEXP Z1SEXP, SEXP Z2SEXP, SEXP rSEXP, SEXP dtSEXP, SEXP interceptSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Z0(Z0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z1(Z1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z2(Z2SEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< bool >::type intercept(interceptSEXP);
    rcpp_result_gen = Rcpp::wrap(johansenCpp(Z0, Z1, Z2, r, dt, intercept));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_narwhalCA_glsCpp", (DL_FUNC) &_narwhalCA_glsCpp, 12},
    {"_narwhalCA_glsConstCpp", (DL_FUNC) &_narwhalCA_glsConstCpp, 12},
    {"_narwhalCA_vecm", (DL_FUNC) &_narwhalCA_vecm, 6},
    {"_narwhalCA_johansenCpp", (DL_FUNC) &_narwhalCA_johansenCpp, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_narwhalCA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}