// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

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
// mvrnormArma5
arma::mat mvrnormArma5(int n, arma::mat invSigma);
RcppExport SEXP _Rcppocc_mvrnormArma5(SEXP nSEXP, SEXP invSigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type invSigma(invSigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrnormArma5(n, invSigma));
    return rcpp_result_gen;
END_RCPP
}
// mvrnormArma6
arma::mat mvrnormArma6(int n, arma::vec A, arma::mat invSigma);
RcppExport SEXP _Rcppocc_mvrnormArma6(SEXP nSEXP, SEXP ASEXP, SEXP invSigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type invSigma(invSigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrnormArma6(n, A, invSigma));
    return rcpp_result_gen;
END_RCPP
}
// mvrnormArma4_timers
void mvrnormArma4_timers(int n, arma::vec A, arma::mat invSigma);
RcppExport SEXP _Rcppocc_mvrnormArma4_timers(SEXP nSEXP, SEXP ASEXP, SEXP invSigmaSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type invSigma(invSigmaSEXP);
    mvrnormArma4_timers(n, A, invSigma);
    return R_NilValue;
END_RCPP
}
// setCubeslices
arma::cube setCubeslices(arma::mat K);
RcppExport SEXP _Rcppocc_setCubeslices(SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(setCubeslices(K));
    return rcpp_result_gen;
END_RCPP
}
// setCubeslices_p
arma::cube setCubeslices_p(arma::mat K, int ncores);
RcppExport SEXP _Rcppocc_setCubeslices_p(SEXP KSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(setCubeslices_p(K, ncores));
    return rcpp_result_gen;
END_RCPP
}
// cross_p
arma::mat cross_p(arma::cube S, arma::vec D, int ncores);
RcppExport SEXP _Rcppocc_cross_p(SEXP SSEXP, SEXP DSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type S(SSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(cross_p(S, D, ncores));
    return rcpp_result_gen;
END_RCPP
}
// multiplybyconstants
arma::mat multiplybyconstants(arma::mat S, arma::vec Dv);
RcppExport SEXP _Rcppocc_multiplybyconstants(SEXP SSEXP, SEXP DvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Dv(DvSEXP);
    rcpp_result_gen = Rcpp::wrap(multiplybyconstants(S, Dv));
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
// rpg5_openmp
arma::vec rpg5_openmp(arma::mat scale, int ncores);
RcppExport SEXP _Rcppocc_rpg5_openmp(SEXP scaleSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(rpg5_openmp(scale, ncores));
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
// logitoccPG3_z
List logitoccPG3_z(arma::mat X, arma::mat Y, arma::mat W_vb, NumericVector siteids, int ndraws, arma::vec ysum, arma::vec z, arma::vec nvisits, arma::mat alpha_m, arma::mat beta_m, arma::mat sigma_inv_alpha_p, arma::mat sigma_inv_beta_p, double percent_burn_in);
RcppExport SEXP _Rcppocc_logitoccPG3_z(SEXP XSEXP, SEXP YSEXP, SEXP W_vbSEXP, SEXP siteidsSEXP, SEXP ndrawsSEXP, SEXP ysumSEXP, SEXP zSEXP, SEXP nvisitsSEXP, SEXP alpha_mSEXP, SEXP beta_mSEXP, SEXP sigma_inv_alpha_pSEXP, SEXP sigma_inv_beta_pSEXP, SEXP percent_burn_inSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(logitoccPG3_z(X, Y, W_vb, siteids, ndraws, ysum, z, nvisits, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, percent_burn_in));
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
// seq_int
arma::vec seq_int(long int a, long int b);
RcppExport SEXP _Rcppocc_seq_int(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< long int >::type a(aSEXP);
    Rcpp::traits::input_parameter< long int >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(seq_int(a, b));
    return rcpp_result_gen;
END_RCPP
}
// quadform2
arma::mat quadform2(arma::mat X, arma::vec dmat);
RcppExport SEXP _Rcppocc_quadform2(SEXP XSEXP, SEXP dmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dmat(dmatSEXP);
    rcpp_result_gen = Rcpp::wrap(quadform2(X, dmat));
    return rcpp_result_gen;
END_RCPP
}
// matrix_multiplication
arma::mat matrix_multiplication(arma::mat x, arma::mat y, int ncores);
RcppExport SEXP _Rcppocc_matrix_multiplication(SEXP xSEXP, SEXP ySEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_multiplication(x, y, ncores));
    return rcpp_result_gen;
END_RCPP
}
// quadform3
arma::mat quadform3(arma::mat X, arma::vec dmat, int ncores);
RcppExport SEXP _Rcppocc_quadform3(SEXP XSEXP, SEXP dmatSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dmat(dmatSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(quadform3(X, dmat, ncores));
    return rcpp_result_gen;
END_RCPP
}
// xtx
arma::mat xtx(arma::mat x, int ncores);
RcppExport SEXP _Rcppocc_xtx(SEXP xSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(xtx(x, ncores));
    return rcpp_result_gen;
END_RCPP
}
// quadform4
arma::mat quadform4(arma::mat X, arma::vec dmat, int ncores);
RcppExport SEXP _Rcppocc_quadform4(SEXP XSEXP, SEXP dmatSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dmat(dmatSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(quadform4(X, dmat, ncores));
    return rcpp_result_gen;
END_RCPP
}
// rbinom2
int rbinom2(arma::vec prob);
RcppExport SEXP _Rcppocc_rbinom2(SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(rbinom2(prob));
    return rcpp_result_gen;
END_RCPP
}
// rbinom3
int rbinom3(double prob);
RcppExport SEXP _Rcppocc_rbinom3(SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(rbinom3(prob));
    return rcpp_result_gen;
END_RCPP
}
// matrix_multiplication3
arma::mat matrix_multiplication3(arma::mat x, arma::mat y);
RcppExport SEXP _Rcppocc_matrix_multiplication3(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_multiplication3(x, y));
    return rcpp_result_gen;
END_RCPP
}
// matrix_multiplication4
arma::mat matrix_multiplication4(arma::mat x);
RcppExport SEXP _Rcppocc_matrix_multiplication4(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_multiplication4(x));
    return rcpp_result_gen;
END_RCPP
}
// transpose_sq
arma::mat transpose_sq(arma::mat x);
RcppExport SEXP _Rcppocc_transpose_sq(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(transpose_sq(x));
    return rcpp_result_gen;
END_RCPP
}
// transpose_notsq
arma::mat transpose_notsq(arma::mat x, int ncores);
RcppExport SEXP _Rcppocc_transpose_notsq(SEXP xSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(transpose_notsq(x, ncores));
    return rcpp_result_gen;
END_RCPP
}
// Transpose2
arma::mat Transpose2(arma::mat x);
RcppExport SEXP _Rcppocc_Transpose2(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(Transpose2(x));
    return rcpp_result_gen;
END_RCPP
}
// rmatmult2
arma::mat rmatmult2(arma::mat X, arma::vec Dv, int nsamps, int ncores);
RcppExport SEXP _Rcppocc_rmatmult2(SEXP XSEXP, SEXP DvSEXP, SEXP nsampsSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Dv(DvSEXP);
    Rcpp::traits::input_parameter< int >::type nsamps(nsampsSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(rmatmult2(X, Dv, nsamps, ncores));
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
// logitoccSPAT2
List logitoccSPAT2(arma::mat X, arma::mat W_vb, arma::mat Y, arma::mat z, arma::vec ysum, arma::vec nvisits, arma::mat K, arma::mat Minv, double n_obs, NumericVector siteids, arma::vec unsurveyed_ind, double tau_0, double a_tau, double b_tau, arma::mat alpha_m, arma::mat beta_m, arma::mat sigma_inv_alpha_p, arma::mat sigma_inv_beta_p, int ndraws, double percent_burn_in, int ncores);
RcppExport SEXP _Rcppocc_logitoccSPAT2(SEXP XSEXP, SEXP W_vbSEXP, SEXP YSEXP, SEXP zSEXP, SEXP ysumSEXP, SEXP nvisitsSEXP, SEXP KSEXP, SEXP MinvSEXP, SEXP n_obsSEXP, SEXP siteidsSEXP, SEXP unsurveyed_indSEXP, SEXP tau_0SEXP, SEXP a_tauSEXP, SEXP b_tauSEXP, SEXP alpha_mSEXP, SEXP beta_mSEXP, SEXP sigma_inv_alpha_pSEXP, SEXP sigma_inv_beta_pSEXP, SEXP ndrawsSEXP, SEXP percent_burn_inSEXP, SEXP ncoresSEXP) {
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
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(logitoccSPAT2(X, W_vb, Y, z, ysum, nvisits, K, Minv, n_obs, siteids, unsurveyed_ind, tau_0, a_tau, b_tau, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, ndraws, percent_burn_in, ncores));
    return rcpp_result_gen;
END_RCPP
}
// logitoccSPAT3
List logitoccSPAT3(arma::mat X, arma::mat W_vb, arma::mat Y, arma::mat z, arma::vec ysum, arma::vec nvisits, arma::mat K, arma::mat Minv, double n_obs, NumericVector siteids, arma::vec unsurveyed_ind, double tau_0, double a_tau, double b_tau, arma::mat alpha_m, arma::mat beta_m, arma::mat sigma_inv_alpha_p, arma::mat sigma_inv_beta_p, int ndraws, double percent_burn_in);
RcppExport SEXP _Rcppocc_logitoccSPAT3(SEXP XSEXP, SEXP W_vbSEXP, SEXP YSEXP, SEXP zSEXP, SEXP ysumSEXP, SEXP nvisitsSEXP, SEXP KSEXP, SEXP MinvSEXP, SEXP n_obsSEXP, SEXP siteidsSEXP, SEXP unsurveyed_indSEXP, SEXP tau_0SEXP, SEXP a_tauSEXP, SEXP b_tauSEXP, SEXP alpha_mSEXP, SEXP beta_mSEXP, SEXP sigma_inv_alpha_pSEXP, SEXP sigma_inv_beta_pSEXP, SEXP ndrawsSEXP, SEXP percent_burn_inSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(logitoccSPAT3(X, W_vb, Y, z, ysum, nvisits, K, Minv, n_obs, siteids, unsurveyed_ind, tau_0, a_tau, b_tau, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, ndraws, percent_burn_in));
    return rcpp_result_gen;
END_RCPP
}
// logitoccSPATsplit
List logitoccSPATsplit(arma::mat X, arma::mat W_vb, arma::mat Y, arma::mat z, arma::vec ysum, arma::vec nvisits, arma::mat K, arma::mat Minv, double n_obs, NumericVector siteids, arma::vec unsurveyed_ind, double tau_0, double a_tau, double b_tau, arma::mat alpha_m, arma::mat beta_m, arma::mat sigma_inv_alpha_p, arma::mat sigma_inv_beta_p, int ndraws, double percent_burn_in);
RcppExport SEXP _Rcppocc_logitoccSPATsplit(SEXP XSEXP, SEXP W_vbSEXP, SEXP YSEXP, SEXP zSEXP, SEXP ysumSEXP, SEXP nvisitsSEXP, SEXP KSEXP, SEXP MinvSEXP, SEXP n_obsSEXP, SEXP siteidsSEXP, SEXP unsurveyed_indSEXP, SEXP tau_0SEXP, SEXP a_tauSEXP, SEXP b_tauSEXP, SEXP alpha_mSEXP, SEXP beta_mSEXP, SEXP sigma_inv_alpha_pSEXP, SEXP sigma_inv_beta_pSEXP, SEXP ndrawsSEXP, SEXP percent_burn_inSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(logitoccSPATsplit(X, W_vb, Y, z, ysum, nvisits, K, Minv, n_obs, siteids, unsurveyed_ind, tau_0, a_tau, b_tau, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, ndraws, percent_burn_in));
    return rcpp_result_gen;
END_RCPP
}
// rpg6
double rpg6(int n, double scale);
RcppExport SEXP _Rcppocc_rpg6(SEXP nSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(rpg6(n, scale));
    return rcpp_result_gen;
END_RCPP
}
// rpg7
arma::vec rpg7(arma::vec n, arma::mat scale);
RcppExport SEXP _Rcppocc_rpg7(SEXP nSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(rpg7(n, scale));
    return rcpp_result_gen;
END_RCPP
}
// rnorm6
double rnorm6(double a, double c);
RcppExport SEXP _Rcppocc_rnorm6(SEXP aSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(rnorm6(a, c));
    return rcpp_result_gen;
END_RCPP
}
// rnorm_1
double rnorm_1(double mu, double sd);
RcppExport SEXP _Rcppocc_rnorm_1(SEXP muSEXP, SEXP sdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    rcpp_result_gen = Rcpp::wrap(rnorm_1(mu, sd));
    return rcpp_result_gen;
END_RCPP
}
// rnorm7
arma::vec rnorm7(arma::vec n, arma::mat scale);
RcppExport SEXP _Rcppocc_rnorm7(SEXP nSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(rnorm7(n, scale));
    return rcpp_result_gen;
END_RCPP
}
// logitoccSPATBINOM
List logitoccSPATBINOM(arma::mat X, arma::mat W_vb, arma::mat Y, arma::mat z, arma::vec ysum, arma::vec nvisits, arma::mat K, arma::mat Minv, double n_obs, NumericVector siteids, arma::vec unsurveyed_ind, double tau_0, double a_tau, double b_tau, arma::mat alpha_m, arma::mat beta_m, arma::mat sigma_inv_alpha_p, arma::mat sigma_inv_beta_p, int ndraws, double percent_burn_in);
RcppExport SEXP _Rcppocc_logitoccSPATBINOM(SEXP XSEXP, SEXP W_vbSEXP, SEXP YSEXP, SEXP zSEXP, SEXP ysumSEXP, SEXP nvisitsSEXP, SEXP KSEXP, SEXP MinvSEXP, SEXP n_obsSEXP, SEXP siteidsSEXP, SEXP unsurveyed_indSEXP, SEXP tau_0SEXP, SEXP a_tauSEXP, SEXP b_tauSEXP, SEXP alpha_mSEXP, SEXP beta_mSEXP, SEXP sigma_inv_alpha_pSEXP, SEXP sigma_inv_beta_pSEXP, SEXP ndrawsSEXP, SEXP percent_burn_inSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(logitoccSPATBINOM(X, W_vb, Y, z, ysum, nvisits, K, Minv, n_obs, siteids, unsurveyed_ind, tau_0, a_tau, b_tau, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, ndraws, percent_burn_in));
    return rcpp_result_gen;
END_RCPP
}
// logitoccSPATBINOMPG
List logitoccSPATBINOMPG(arma::mat X, arma::mat W_vb, arma::mat Y, arma::mat z, arma::vec ysum, arma::vec nvisits, arma::mat K, arma::mat Minv, double n_obs, NumericVector siteids, arma::vec unsurveyed_ind, double tau_0, double a_tau, double b_tau, arma::mat alpha_m, arma::mat beta_m, arma::mat sigma_inv_alpha_p, arma::mat sigma_inv_beta_p, int ndraws, double percent_burn_in);
RcppExport SEXP _Rcppocc_logitoccSPATBINOMPG(SEXP XSEXP, SEXP W_vbSEXP, SEXP YSEXP, SEXP zSEXP, SEXP ysumSEXP, SEXP nvisitsSEXP, SEXP KSEXP, SEXP MinvSEXP, SEXP n_obsSEXP, SEXP siteidsSEXP, SEXP unsurveyed_indSEXP, SEXP tau_0SEXP, SEXP a_tauSEXP, SEXP b_tauSEXP, SEXP alpha_mSEXP, SEXP beta_mSEXP, SEXP sigma_inv_alpha_pSEXP, SEXP sigma_inv_beta_pSEXP, SEXP ndrawsSEXP, SEXP percent_burn_inSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(logitoccSPATBINOMPG(X, W_vb, Y, z, ysum, nvisits, K, Minv, n_obs, siteids, unsurveyed_ind, tau_0, a_tau, b_tau, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, ndraws, percent_burn_in));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Rcppocc_mvrnormArma4", (DL_FUNC) &_Rcppocc_mvrnormArma4, 3},
    {"_Rcppocc_mvrnormArma5", (DL_FUNC) &_Rcppocc_mvrnormArma5, 2},
    {"_Rcppocc_mvrnormArma6", (DL_FUNC) &_Rcppocc_mvrnormArma6, 3},
    {"_Rcppocc_mvrnormArma4_timers", (DL_FUNC) &_Rcppocc_mvrnormArma4_timers, 3},
    {"_Rcppocc_setCubeslices", (DL_FUNC) &_Rcppocc_setCubeslices, 1},
    {"_Rcppocc_setCubeslices_p", (DL_FUNC) &_Rcppocc_setCubeslices_p, 2},
    {"_Rcppocc_cross_p", (DL_FUNC) &_Rcppocc_cross_p, 3},
    {"_Rcppocc_multiplybyconstants", (DL_FUNC) &_Rcppocc_multiplybyconstants, 2},
    {"_Rcppocc_logitoccDA4", (DL_FUNC) &_Rcppocc_logitoccDA4, 12},
    {"_Rcppocc_rpg5", (DL_FUNC) &_Rcppocc_rpg5, 1},
    {"_Rcppocc_rpg5_openmp", (DL_FUNC) &_Rcppocc_rpg5_openmp, 2},
    {"_Rcppocc_logitoccPG3", (DL_FUNC) &_Rcppocc_logitoccPG3, 13},
    {"_Rcppocc_logitoccPG3_z", (DL_FUNC) &_Rcppocc_logitoccPG3_z, 13},
    {"_Rcppocc_quadform", (DL_FUNC) &_Rcppocc_quadform, 3},
    {"_Rcppocc_seq_int", (DL_FUNC) &_Rcppocc_seq_int, 2},
    {"_Rcppocc_quadform2", (DL_FUNC) &_Rcppocc_quadform2, 2},
    {"_Rcppocc_matrix_multiplication", (DL_FUNC) &_Rcppocc_matrix_multiplication, 3},
    {"_Rcppocc_quadform3", (DL_FUNC) &_Rcppocc_quadform3, 3},
    {"_Rcppocc_xtx", (DL_FUNC) &_Rcppocc_xtx, 2},
    {"_Rcppocc_quadform4", (DL_FUNC) &_Rcppocc_quadform4, 3},
    {"_Rcppocc_rbinom2", (DL_FUNC) &_Rcppocc_rbinom2, 1},
    {"_Rcppocc_rbinom3", (DL_FUNC) &_Rcppocc_rbinom3, 1},
    {"_Rcppocc_matrix_multiplication3", (DL_FUNC) &_Rcppocc_matrix_multiplication3, 2},
    {"_Rcppocc_matrix_multiplication4", (DL_FUNC) &_Rcppocc_matrix_multiplication4, 1},
    {"_Rcppocc_transpose_sq", (DL_FUNC) &_Rcppocc_transpose_sq, 1},
    {"_Rcppocc_transpose_notsq", (DL_FUNC) &_Rcppocc_transpose_notsq, 2},
    {"_Rcppocc_Transpose2", (DL_FUNC) &_Rcppocc_Transpose2, 1},
    {"_Rcppocc_rmatmult2", (DL_FUNC) &_Rcppocc_rmatmult2, 4},
    {"_Rcppocc_logitoccSPAT", (DL_FUNC) &_Rcppocc_logitoccSPAT, 20},
    {"_Rcppocc_logitoccSPAT2", (DL_FUNC) &_Rcppocc_logitoccSPAT2, 21},
    {"_Rcppocc_logitoccSPAT3", (DL_FUNC) &_Rcppocc_logitoccSPAT3, 20},
    {"_Rcppocc_logitoccSPATsplit", (DL_FUNC) &_Rcppocc_logitoccSPATsplit, 20},
    {"_Rcppocc_rpg6", (DL_FUNC) &_Rcppocc_rpg6, 2},
    {"_Rcppocc_rpg7", (DL_FUNC) &_Rcppocc_rpg7, 2},
    {"_Rcppocc_rnorm6", (DL_FUNC) &_Rcppocc_rnorm6, 2},
    {"_Rcppocc_rnorm_1", (DL_FUNC) &_Rcppocc_rnorm_1, 2},
    {"_Rcppocc_rnorm7", (DL_FUNC) &_Rcppocc_rnorm7, 2},
    {"_Rcppocc_logitoccSPATBINOM", (DL_FUNC) &_Rcppocc_logitoccSPATBINOM, 20},
    {"_Rcppocc_logitoccSPATBINOMPG", (DL_FUNC) &_Rcppocc_logitoccSPATBINOMPG, 20},
    {NULL, NULL, 0}
};

RcppExport void R_init_Rcppocc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
