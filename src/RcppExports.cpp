// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Q_theta_cpp
double Q_theta_cpp(NumericVector para_vec_ori, NumericVector para_vec_0_ori, NumericVector yi_vec, NumericVector m_star_vec, NumericVector x_i_vec, NumericVector l_i_vec, NumericMatrix confound_mat, bool x4_inter, bool x5_inter);
RcppExport SEXP _MarZIC_Q_theta_cpp(SEXP para_vec_oriSEXP, SEXP para_vec_0_oriSEXP, SEXP yi_vecSEXP, SEXP m_star_vecSEXP, SEXP x_i_vecSEXP, SEXP l_i_vecSEXP, SEXP confound_matSEXP, SEXP x4_interSEXP, SEXP x5_interSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type para_vec_ori(para_vec_oriSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type para_vec_0_ori(para_vec_0_oriSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yi_vec(yi_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type m_star_vec(m_star_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x_i_vec(x_i_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type l_i_vec(l_i_vecSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type confound_mat(confound_matSEXP);
    Rcpp::traits::input_parameter< bool >::type x4_inter(x4_interSEXP);
    Rcpp::traits::input_parameter< bool >::type x5_inter(x5_interSEXP);
    rcpp_result_gen = Rcpp::wrap(Q_theta_cpp(para_vec_ori, para_vec_0_ori, yi_vec, m_star_vec, x_i_vec, l_i_vec, confound_mat, x4_inter, x5_inter));
    return rcpp_result_gen;
END_RCPP
}
// Q_theta_cpp_nz
double Q_theta_cpp_nz(NumericVector para_vec_ori, NumericVector para_vec_0_ori, NumericVector yi_vec, NumericVector m_star_vec, NumericVector x_i_vec, NumericVector l_i_vec, NumericMatrix confound_mat, bool x4_inter, bool x5_inter);
RcppExport SEXP _MarZIC_Q_theta_cpp_nz(SEXP para_vec_oriSEXP, SEXP para_vec_0_oriSEXP, SEXP yi_vecSEXP, SEXP m_star_vecSEXP, SEXP x_i_vecSEXP, SEXP l_i_vecSEXP, SEXP confound_matSEXP, SEXP x4_interSEXP, SEXP x5_interSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type para_vec_ori(para_vec_oriSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type para_vec_0_ori(para_vec_0_oriSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yi_vec(yi_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type m_star_vec(m_star_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x_i_vec(x_i_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type l_i_vec(l_i_vecSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type confound_mat(confound_matSEXP);
    Rcpp::traits::input_parameter< bool >::type x4_inter(x4_interSEXP);
    Rcpp::traits::input_parameter< bool >::type x5_inter(x5_interSEXP);
    rcpp_result_gen = Rcpp::wrap(Q_theta_cpp_nz(para_vec_ori, para_vec_0_ori, yi_vec, m_star_vec, x_i_vec, l_i_vec, confound_mat, x4_inter, x5_inter));
    return rcpp_result_gen;
END_RCPP
}
// Q_theta_cpp_nomix
double Q_theta_cpp_nomix(NumericVector para_vec_ori, NumericVector para_vec_0_ori, NumericVector yi_vec, NumericVector m_star_vec, NumericVector x_i_vec, NumericVector l_i_vec, NumericMatrix confound_mat, bool x4_inter, bool x5_inter);
RcppExport SEXP _MarZIC_Q_theta_cpp_nomix(SEXP para_vec_oriSEXP, SEXP para_vec_0_oriSEXP, SEXP yi_vecSEXP, SEXP m_star_vecSEXP, SEXP x_i_vecSEXP, SEXP l_i_vecSEXP, SEXP confound_matSEXP, SEXP x4_interSEXP, SEXP x5_interSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type para_vec_ori(para_vec_oriSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type para_vec_0_ori(para_vec_0_oriSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yi_vec(yi_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type m_star_vec(m_star_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x_i_vec(x_i_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type l_i_vec(l_i_vecSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type confound_mat(confound_matSEXP);
    Rcpp::traits::input_parameter< bool >::type x4_inter(x4_interSEXP);
    Rcpp::traits::input_parameter< bool >::type x5_inter(x5_interSEXP);
    rcpp_result_gen = Rcpp::wrap(Q_theta_cpp_nomix(para_vec_ori, para_vec_0_ori, yi_vec, m_star_vec, x_i_vec, l_i_vec, confound_mat, x4_inter, x5_inter));
    return rcpp_result_gen;
END_RCPP
}
// Q_theta_cpp_nomix_nz
double Q_theta_cpp_nomix_nz(NumericVector para_vec_ori, NumericVector para_vec_0_ori, NumericVector yi_vec, NumericVector m_star_vec, NumericVector x_i_vec, NumericVector l_i_vec, NumericMatrix confound_mat, bool x4_inter, bool x5_inter);
RcppExport SEXP _MarZIC_Q_theta_cpp_nomix_nz(SEXP para_vec_oriSEXP, SEXP para_vec_0_oriSEXP, SEXP yi_vecSEXP, SEXP m_star_vecSEXP, SEXP x_i_vecSEXP, SEXP l_i_vecSEXP, SEXP confound_matSEXP, SEXP x4_interSEXP, SEXP x5_interSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type para_vec_ori(para_vec_oriSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type para_vec_0_ori(para_vec_0_oriSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yi_vec(yi_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type m_star_vec(m_star_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x_i_vec(x_i_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type l_i_vec(l_i_vecSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type confound_mat(confound_matSEXP);
    Rcpp::traits::input_parameter< bool >::type x4_inter(x4_interSEXP);
    Rcpp::traits::input_parameter< bool >::type x5_inter(x5_interSEXP);
    rcpp_result_gen = Rcpp::wrap(Q_theta_cpp_nomix_nz(para_vec_ori, para_vec_0_ori, yi_vec, m_star_vec, x_i_vec, l_i_vec, confound_mat, x4_inter, x5_inter));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MarZIC_Q_theta_cpp", (DL_FUNC) &_MarZIC_Q_theta_cpp, 9},
    {"_MarZIC_Q_theta_cpp_nz", (DL_FUNC) &_MarZIC_Q_theta_cpp_nz, 9},
    {"_MarZIC_Q_theta_cpp_nomix", (DL_FUNC) &_MarZIC_Q_theta_cpp_nomix, 9},
    {"_MarZIC_Q_theta_cpp_nomix_nz", (DL_FUNC) &_MarZIC_Q_theta_cpp_nomix_nz, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_MarZIC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
