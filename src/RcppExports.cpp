// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// mvrnormArma
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);
RcppExport SEXP _Rcppocc_mvrnormArma(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrnormArma(n, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// mvrnormArma2
arma::mat mvrnormArma2(int n, arma::vec mu, arma::mat sigma);
RcppExport SEXP _Rcppocc_mvrnormArma2(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrnormArma2(n, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// mvrnormArma3
arma::mat mvrnormArma3(int n, arma::vec A, arma::mat invSigma);
RcppExport SEXP _Rcppocc_mvrnormArma3(SEXP nSEXP, SEXP ASEXP, SEXP invSigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type invSigma(invSigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrnormArma3(n, A, invSigma));
    return rcpp_result_gen;
END_RCPP
}
// mvrnormArma4
arma::mat mvrnormArma4(int n, arma::vec A, arma::mat invSigma);
RcppExport SEXP _Rcppocc_mvrnormArma4(SEXP nSEXP, SEXP ASEXP, SEXP invSigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type invSigma(invSigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrnormArma4(n, A, invSigma));
    return rcpp_result_gen;
END_RCPP
}
// posterior_r3
double posterior_r3(double x, arma::vec wr, arma::vec sr);
RcppExport SEXP _Rcppocc_posterior_r3(SEXP xSEXP, SEXP wrSEXP, SEXP srSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type wr(wrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sr(srSEXP);
    rcpp_result_gen = Rcpp::wrap(posterior_r3(x, wr, sr));
    return rcpp_result_gen;
END_RCPP
}
// invlogit
arma::vec invlogit(arma::mat lincomb);
RcppExport SEXP _Rcppocc_invlogit(SEXP lincombSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type lincomb(lincombSEXP);
    rcpp_result_gen = Rcpp::wrap(invlogit(lincomb));
    return rcpp_result_gen;
END_RCPP
}
// ln_invlogit
arma::vec ln_invlogit(arma::mat lincomb);
RcppExport SEXP _Rcppocc_ln_invlogit(SEXP lincombSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type lincomb(lincombSEXP);
    rcpp_result_gen = Rcpp::wrap(ln_invlogit(lincomb));
    return rcpp_result_gen;
END_RCPP
}
// log_not_prob
arma::vec log_not_prob(arma::mat lincomb);
RcppExport SEXP _Rcppocc_log_not_prob(SEXP lincombSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type lincomb(lincombSEXP);
    rcpp_result_gen = Rcpp::wrap(log_not_prob(lincomb));
    return rcpp_result_gen;
END_RCPP
}
// logitoccDA4
List logitoccDA4(arma::mat X, arma::mat Y, arma::mat W_vb, NumericVector siteids, int ndraws, arma::vec ysum, arma::vec z, arma::vec nvisits, arma::mat alpha_m, arma::mat beta_m, arma::mat sigma_inv_alpha_p, arma::mat sigma_inv_beta_p);
RcppExport SEXP _Rcppocc_logitoccDA4(SEXP XSEXP, SEXP YSEXP, SEXP W_vbSEXP, SEXP siteidsSEXP, SEXP ndrawsSEXP, SEXP ysumSEXP, SEXP zSEXP, SEXP nvisitsSEXP, SEXP alpha_mSEXP, SEXP beta_mSEXP, SEXP sigma_inv_alpha_pSEXP, SEXP sigma_inv_beta_pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W_vb(W_vbSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type siteids(siteidsSEXP);
    Rcpp::traits::input_parameter< int >::type ndraws(ndrawsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ysum(ysumSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nvisits(nvisitsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type alpha_m(alpha_mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_m(beta_mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_inv_alpha_p(sigma_inv_alpha_pSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_inv_beta_p(sigma_inv_beta_pSEXP);
    rcpp_result_gen = Rcpp::wrap(logitoccDA4(X, Y, W_vb, siteids, ndraws, ysum, z, nvisits, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p));
    return rcpp_result_gen;
END_RCPP
}
// rpg2
arma::vec rpg2(arma::mat scale);
RcppExport SEXP _Rcppocc_rpg2(SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(rpg2(scale));
    return rcpp_result_gen;
END_RCPP
}
// rpg5
arma::vec rpg5(arma::mat scale);
RcppExport SEXP _Rcppocc_rpg5(SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(rpg5(scale));
    return rcpp_result_gen;
END_RCPP
}
// parallelMatrixRpg
NumericMatrix parallelMatrixRpg(NumericMatrix x);
RcppExport SEXP _Rcppocc_parallelMatrixRpg(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(parallelMatrixRpg(x));
    return rcpp_result_gen;
END_RCPP
}
// logitoccPG3
List logitoccPG3(arma::mat X, arma::mat Y, arma::mat W_vb, NumericVector siteids, int ndraws, arma::vec ysum, arma::vec z, arma::vec nvisits, arma::mat alpha_m, arma::mat beta_m, arma::mat sigma_inv_alpha_p, arma::mat sigma_inv_beta_p, double percent_burn_in);
RcppExport SEXP _Rcppocc_logitoccPG3(SEXP XSEXP, SEXP YSEXP, SEXP W_vbSEXP, SEXP siteidsSEXP, SEXP ndrawsSEXP, SEXP ysumSEXP, SEXP zSEXP, SEXP nvisitsSEXP, SEXP alpha_mSEXP, SEXP beta_mSEXP, SEXP sigma_inv_alpha_pSEXP, SEXP sigma_inv_beta_pSEXP, SEXP percent_burn_inSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W_vb(W_vbSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type siteids(siteidsSEXP);
    Rcpp::traits::input_parameter< int >::type ndraws(ndrawsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ysum(ysumSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nvisits(nvisitsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type alpha_m(alpha_mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_m(beta_mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_inv_alpha_p(sigma_inv_alpha_pSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_inv_beta_p(sigma_inv_beta_pSEXP);
    Rcpp::traits::input_parameter< double >::type percent_burn_in(percent_burn_inSEXP);
    rcpp_result_gen = Rcpp::wrap(logitoccPG3(X, Y, W_vb, siteids, ndraws, ysum, z, nvisits, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, percent_burn_in));
    return rcpp_result_gen;
END_RCPP
}
// rgammadouble
double rgammadouble(int a, double b, double c);
RcppExport SEXP _Rcppocc_rgammadouble(SEXP aSEXP, SEXP bSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(rgammadouble(a, b, c));
    return rcpp_result_gen;
END_RCPP
}
// rnormdouble
double rnormdouble(double b, double c);
RcppExport SEXP _Rcppocc_rnormdouble(SEXP bSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(rnormdouble(b, c));
    return rcpp_result_gen;
END_RCPP
}
// quadform
arma::mat quadform(arma::mat X, arma::mat Xt, arma::vec dmat);
RcppExport SEXP _Rcppocc_quadform(SEXP XSEXP, SEXP XtSEXP, SEXP dmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xt(XtSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dmat(dmatSEXP);
    rcpp_result_gen = Rcpp::wrap(quadform(X, Xt, dmat));
    return rcpp_result_gen;
END_RCPP
}
// diagtimesX
arma::mat diagtimesX(arma::mat X, arma::vec dmat);
RcppExport SEXP _Rcppocc_diagtimesX(SEXP XSEXP, SEXP dmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dmat(dmatSEXP);
    rcpp_result_gen = Rcpp::wrap(diagtimesX(X, dmat));
    return rcpp_result_gen;
END_RCPP
}
// logitoccSPAT
List logitoccSPAT(arma::mat X, arma::mat W_vb, arma::mat Y, arma::mat z, arma::vec ysum, arma::vec nvisits, arma::mat K, arma::mat Minv, double n_obs, NumericVector siteids, arma::vec unsurveyed_ind, double tau_0, double a_tau, double b_tau, arma::mat alpha_m, arma::mat beta_m, arma::mat sigma_inv_alpha_p, arma::mat sigma_inv_beta_p, int ndraws, double percent_burn_in);
RcppExport SEXP _Rcppocc_logitoccSPAT(SEXP XSEXP, SEXP W_vbSEXP, SEXP YSEXP, SEXP zSEXP, SEXP ysumSEXP, SEXP nvisitsSEXP, SEXP KSEXP, SEXP MinvSEXP, SEXP n_obsSEXP, SEXP siteidsSEXP, SEXP unsurveyed_indSEXP, SEXP tau_0SEXP, SEXP a_tauSEXP, SEXP b_tauSEXP, SEXP alpha_mSEXP, SEXP beta_mSEXP, SEXP sigma_inv_alpha_pSEXP, SEXP sigma_inv_beta_pSEXP, SEXP ndrawsSEXP, SEXP percent_burn_inSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W_vb(W_vbSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ysum(ysumSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nvisits(nvisitsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Minv(MinvSEXP);
    Rcpp::traits::input_parameter< double >::type n_obs(n_obsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type siteids(siteidsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type unsurveyed_ind(unsurveyed_indSEXP);
    Rcpp::traits::input_parameter< double >::type tau_0(tau_0SEXP);
    Rcpp::traits::input_parameter< double >::type a_tau(a_tauSEXP);
    Rcpp::traits::input_parameter< double >::type b_tau(b_tauSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type alpha_m(alpha_mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_m(beta_mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_inv_alpha_p(sigma_inv_alpha_pSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_inv_beta_p(sigma_inv_beta_pSEXP);
    Rcpp::traits::input_parameter< int >::type ndraws(ndrawsSEXP);
    Rcpp::traits::input_parameter< double >::type percent_burn_in(percent_burn_inSEXP);
    rcpp_result_gen = Rcpp::wrap(logitoccSPAT(X, W_vb, Y, z, ysum, nvisits, K, Minv, n_obs, siteids, unsurveyed_ind, tau_0, a_tau, b_tau, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, ndraws, percent_burn_in));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Rcppocc_mvrnormArma", (DL_FUNC) &_Rcppocc_mvrnormArma, 3},
    {"_Rcppocc_mvrnormArma2", (DL_FUNC) &_Rcppocc_mvrnormArma2, 3},
    {"_Rcppocc_mvrnormArma3", (DL_FUNC) &_Rcppocc_mvrnormArma3, 3},
    {"_Rcppocc_mvrnormArma4", (DL_FUNC) &_Rcppocc_mvrnormArma4, 3},
    {"_Rcppocc_posterior_r3", (DL_FUNC) &_Rcppocc_posterior_r3, 3},
    {"_Rcppocc_invlogit", (DL_FUNC) &_Rcppocc_invlogit, 1},
    {"_Rcppocc_ln_invlogit", (DL_FUNC) &_Rcppocc_ln_invlogit, 1},
    {"_Rcppocc_log_not_prob", (DL_FUNC) &_Rcppocc_log_not_prob, 1},
    {"_Rcppocc_logitoccDA4", (DL_FUNC) &_Rcppocc_logitoccDA4, 12},
    {"_Rcppocc_rpg2", (DL_FUNC) &_Rcppocc_rpg2, 1},
    {"_Rcppocc_rpg5", (DL_FUNC) &_Rcppocc_rpg5, 1},
    {"_Rcppocc_parallelMatrixRpg", (DL_FUNC) &_Rcppocc_parallelMatrixRpg, 1},
    {"_Rcppocc_logitoccPG3", (DL_FUNC) &_Rcppocc_logitoccPG3, 13},
    {"_Rcppocc_rgammadouble", (DL_FUNC) &_Rcppocc_rgammadouble, 3},
    {"_Rcppocc_rnormdouble", (DL_FUNC) &_Rcppocc_rnormdouble, 2},
    {"_Rcppocc_quadform", (DL_FUNC) &_Rcppocc_quadform, 3},
    {"_Rcppocc_diagtimesX", (DL_FUNC) &_Rcppocc_diagtimesX, 2},
    {"_Rcppocc_logitoccSPAT", (DL_FUNC) &_Rcppocc_logitoccSPAT, 20},
    {NULL, NULL, 0}
};

RcppExport void R_init_Rcppocc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
