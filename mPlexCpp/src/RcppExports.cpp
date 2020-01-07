// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// run_mPlex
void run_mPlex(const uint_least32_t& seed_, const uint_least32_t& numReps_, const uint_least32_t& numThreads_, const Rcpp::List& networkParameters_, const Rcpp::List& reproductionReference_, const Rcpp::List& patchReleases_, const Rcpp::NumericMatrix& migrationMale_, const Rcpp::NumericMatrix& migrationFemale_, const Rcpp::List& migrationBatch_, const std::string& reproductionType_, const std::string& outputDirectory_, const bool& verbose_);
RcppExport SEXP _mPlexCpp_run_mPlex(SEXP seed_SEXP, SEXP numReps_SEXP, SEXP numThreads_SEXP, SEXP networkParameters_SEXP, SEXP reproductionReference_SEXP, SEXP patchReleases_SEXP, SEXP migrationMale_SEXP, SEXP migrationFemale_SEXP, SEXP migrationBatch_SEXP, SEXP reproductionType_SEXP, SEXP outputDirectory_SEXP, SEXP verbose_SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const uint_least32_t& >::type seed_(seed_SEXP);
    Rcpp::traits::input_parameter< const uint_least32_t& >::type numReps_(numReps_SEXP);
    Rcpp::traits::input_parameter< const uint_least32_t& >::type numThreads_(numThreads_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type networkParameters_(networkParameters_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type reproductionReference_(reproductionReference_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type patchReleases_(patchReleases_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type migrationMale_(migrationMale_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type migrationFemale_(migrationFemale_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type migrationBatch_(migrationBatch_SEXP);
    Rcpp::traits::input_parameter< const std::string& >::type reproductionType_(reproductionType_SEXP);
    Rcpp::traits::input_parameter< const std::string& >::type outputDirectory_(outputDirectory_SEXP);
    Rcpp::traits::input_parameter< const bool& >::type verbose_(verbose_SEXP);
    run_mPlex(seed_, numReps_, numThreads_, networkParameters_, reproductionReference_, patchReleases_, migrationMale_, migrationFemale_, migrationBatch_, reproductionType_, outputDirectory_, verbose_);
    return R_NilValue;
END_RCPP
}
// calcCos
Rcpp::NumericMatrix calcCos(const Rcpp::NumericMatrix& latLongs, const double& r);
RcppExport SEXP _mPlexCpp_calcCos(SEXP latLongsSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type latLongs(latLongsSEXP);
    Rcpp::traits::input_parameter< const double& >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(calcCos(latLongs, r));
    return rcpp_result_gen;
END_RCPP
}
// calcHaversine
Rcpp::NumericMatrix calcHaversine(const Rcpp::NumericMatrix& latLongs, const double& r);
RcppExport SEXP _mPlexCpp_calcHaversine(SEXP latLongsSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type latLongs(latLongsSEXP);
    Rcpp::traits::input_parameter< const double& >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(calcHaversine(latLongs, r));
    return rcpp_result_gen;
END_RCPP
}
// calcVinSph
Rcpp::NumericMatrix calcVinSph(const Rcpp::NumericMatrix& latLongs, const double& r);
RcppExport SEXP _mPlexCpp_calcVinSph(SEXP latLongsSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type latLongs(latLongsSEXP);
    Rcpp::traits::input_parameter< const double& >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(calcVinSph(latLongs, r));
    return rcpp_result_gen;
END_RCPP
}
// calcVinEll
Rcpp::NumericMatrix calcVinEll(const Rcpp::NumericMatrix& latLongs, const double& a, const double& b, const double& f, const double& eps, const double& iter);
RcppExport SEXP _mPlexCpp_calcVinEll(SEXP latLongsSEXP, SEXP aSEXP, SEXP bSEXP, SEXP fSEXP, SEXP epsSEXP, SEXP iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type latLongs(latLongsSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double& >::type f(fSEXP);
    Rcpp::traits::input_parameter< const double& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const double& >::type iter(iterSEXP);
    rcpp_result_gen = Rcpp::wrap(calcVinEll(latLongs, a, b, f, eps, iter));
    return rcpp_result_gen;
END_RCPP
}
// calcLognormalKernel
Rcpp::NumericMatrix calcLognormalKernel(const Rcpp::NumericMatrix& distMat, const double& meanlog, const double& sdlog);
RcppExport SEXP _mPlexCpp_calcLognormalKernel(SEXP distMatSEXP, SEXP meanlogSEXP, SEXP sdlogSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type distMat(distMatSEXP);
    Rcpp::traits::input_parameter< const double& >::type meanlog(meanlogSEXP);
    Rcpp::traits::input_parameter< const double& >::type sdlog(sdlogSEXP);
    rcpp_result_gen = Rcpp::wrap(calcLognormalKernel(distMat, meanlog, sdlog));
    return rcpp_result_gen;
END_RCPP
}
// calcGammaKernel
Rcpp::NumericMatrix calcGammaKernel(const Rcpp::NumericMatrix& distMat, const double& shape, const double& rate);
RcppExport SEXP _mPlexCpp_calcGammaKernel(SEXP distMatSEXP, SEXP shapeSEXP, SEXP rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type distMat(distMatSEXP);
    Rcpp::traits::input_parameter< const double& >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< const double& >::type rate(rateSEXP);
    rcpp_result_gen = Rcpp::wrap(calcGammaKernel(distMat, shape, rate));
    return rcpp_result_gen;
END_RCPP
}
// calcExpKernel
Rcpp::NumericMatrix calcExpKernel(const Rcpp::NumericMatrix& distMat, const double& rate);
RcppExport SEXP _mPlexCpp_calcExpKernel(SEXP distMatSEXP, SEXP rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type distMat(distMatSEXP);
    Rcpp::traits::input_parameter< const double& >::type rate(rateSEXP);
    rcpp_result_gen = Rcpp::wrap(calcExpKernel(distMat, rate));
    return rcpp_result_gen;
END_RCPP
}
// calcHurdleExpKernel
Rcpp::NumericMatrix calcHurdleExpKernel(const Rcpp::NumericMatrix& distMat, double rate, double p0);
RcppExport SEXP _mPlexCpp_calcHurdleExpKernel(SEXP distMatSEXP, SEXP rateSEXP, SEXP p0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type distMat(distMatSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    Rcpp::traits::input_parameter< double >::type p0(p0SEXP);
    rcpp_result_gen = Rcpp::wrap(calcHurdleExpKernel(distMat, rate, p0));
    return rcpp_result_gen;
END_RCPP
}
// simAgg
void simAgg(const std::vector<std::vector<std::string> >& readFiles_, const std::vector<std::vector<std::string> >& writeFiles_, const std::string& largeFile_, const int& simTime_, const Rcpp::List& genKey_);
RcppExport SEXP _mPlexCpp_simAgg(SEXP readFiles_SEXP, SEXP writeFiles_SEXP, SEXP largeFile_SEXP, SEXP simTime_SEXP, SEXP genKey_SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::vector<std::string> >& >::type readFiles_(readFiles_SEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<std::string> >& >::type writeFiles_(writeFiles_SEXP);
    Rcpp::traits::input_parameter< const std::string& >::type largeFile_(largeFile_SEXP);
    Rcpp::traits::input_parameter< const int& >::type simTime_(simTime_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type genKey_(genKey_SEXP);
    simAgg(readFiles_, writeFiles_, largeFile_, simTime_, genKey_);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mPlexCpp_run_mPlex", (DL_FUNC) &_mPlexCpp_run_mPlex, 12},
    {"_mPlexCpp_calcCos", (DL_FUNC) &_mPlexCpp_calcCos, 2},
    {"_mPlexCpp_calcHaversine", (DL_FUNC) &_mPlexCpp_calcHaversine, 2},
    {"_mPlexCpp_calcVinSph", (DL_FUNC) &_mPlexCpp_calcVinSph, 2},
    {"_mPlexCpp_calcVinEll", (DL_FUNC) &_mPlexCpp_calcVinEll, 6},
    {"_mPlexCpp_calcLognormalKernel", (DL_FUNC) &_mPlexCpp_calcLognormalKernel, 3},
    {"_mPlexCpp_calcGammaKernel", (DL_FUNC) &_mPlexCpp_calcGammaKernel, 3},
    {"_mPlexCpp_calcExpKernel", (DL_FUNC) &_mPlexCpp_calcExpKernel, 2},
    {"_mPlexCpp_calcHurdleExpKernel", (DL_FUNC) &_mPlexCpp_calcHurdleExpKernel, 3},
    {"_mPlexCpp_simAgg", (DL_FUNC) &_mPlexCpp_simAgg, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_mPlexCpp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
