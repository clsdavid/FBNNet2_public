// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// extractGeneStateFromTimeSeriesCube
Rcpp::NumericMatrix extractGeneStateFromTimeSeriesCube(Rcpp::List& timeSeriesCube, Rcpp::IntegerVector temporal);
RcppExport SEXP _FBNNet_extractGeneStateFromTimeSeriesCube(SEXP timeSeriesCubeSEXP, SEXP temporalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type timeSeriesCube(timeSeriesCubeSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type temporal(temporalSEXP);
    rcpp_result_gen = Rcpp::wrap(extractGeneStateFromTimeSeriesCube(timeSeriesCube, temporal));
    return rcpp_result_gen;
END_RCPP
}
// getGenePrababilities_basic
Rcpp::List getGenePrababilities_basic(Rcpp::Environment& main_parameters_in_ref, Rcpp::Nullable<Rcpp::List>& fixedgenestate, Rcpp::CharacterVector& target_gene, Rcpp::CharacterVector& new_conditional_gene, Rcpp::IntegerVector temporal, Rcpp::Nullable<Rcpp::List> targetCounts);
RcppExport SEXP _FBNNet_getGenePrababilities_basic(SEXP main_parameters_in_refSEXP, SEXP fixedgenestateSEXP, SEXP target_geneSEXP, SEXP new_conditional_geneSEXP, SEXP temporalSEXP, SEXP targetCountsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Environment& >::type main_parameters_in_ref(main_parameters_in_refSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List>& >::type fixedgenestate(fixedgenestateSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector& >::type target_gene(target_geneSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector& >::type new_conditional_gene(new_conditional_geneSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type temporal(temporalSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type targetCounts(targetCountsSEXP);
    rcpp_result_gen = Rcpp::wrap(getGenePrababilities_basic(main_parameters_in_ref, fixedgenestate, target_gene, new_conditional_gene, temporal, targetCounts));
    return rcpp_result_gen;
END_RCPP
}
// getGenePrababilities_advanced
Rcpp::List getGenePrababilities_advanced(Rcpp::List& getGenePrababilities_basic);
RcppExport SEXP _FBNNet_getGenePrababilities_advanced(SEXP getGenePrababilities_basicSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type getGenePrababilities_basic(getGenePrababilities_basicSEXP);
    rcpp_result_gen = Rcpp::wrap(getGenePrababilities_advanced(getGenePrababilities_basic));
    return rcpp_result_gen;
END_RCPP
}
// getGenePrababilities
Rcpp::List getGenePrababilities(Rcpp::Environment& main_parameters_in_ref, Rcpp::Nullable<Rcpp::List>& fixedgenestate, Rcpp::CharacterVector& target_gene, Rcpp::CharacterVector& new_conditional_gene, Rcpp::IntegerVector temporal, Rcpp::Nullable<Rcpp::List> targetCounts);
RcppExport SEXP _FBNNet_getGenePrababilities(SEXP main_parameters_in_refSEXP, SEXP fixedgenestateSEXP, SEXP target_geneSEXP, SEXP new_conditional_geneSEXP, SEXP temporalSEXP, SEXP targetCountsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Environment& >::type main_parameters_in_ref(main_parameters_in_refSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List>& >::type fixedgenestate(fixedgenestateSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector& >::type target_gene(target_geneSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector& >::type new_conditional_gene(new_conditional_geneSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type temporal(temporalSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type targetCounts(targetCountsSEXP);
    rcpp_result_gen = Rcpp::wrap(getGenePrababilities(main_parameters_in_ref, fixedgenestate, target_gene, new_conditional_gene, temporal, targetCounts));
    return rcpp_result_gen;
END_RCPP
}
// networkFiltering
Rcpp::List networkFiltering(Rcpp::List& res);
RcppExport SEXP _FBNNet_networkFiltering(SEXP resSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type res(resSEXP);
    rcpp_result_gen = Rcpp::wrap(networkFiltering(res));
    return rcpp_result_gen;
END_RCPP
}
// getGenePrababilities_measurements
Rcpp::List getGenePrababilities_measurements(Rcpp::CharacterVector& targetGene, Rcpp::Environment& mainParameters, Rcpp::CharacterVector& genes, Rcpp::Nullable<Rcpp::List>& matchedgenes, Rcpp::IntegerVector temporal, Nullable<Rcpp::List> targetCounts);
RcppExport SEXP _FBNNet_getGenePrababilities_measurements(SEXP targetGeneSEXP, SEXP mainParametersSEXP, SEXP genesSEXP, SEXP matchedgenesSEXP, SEXP temporalSEXP, SEXP targetCountsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::CharacterVector& >::type targetGene(targetGeneSEXP);
    Rcpp::traits::input_parameter< Rcpp::Environment& >::type mainParameters(mainParametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector& >::type genes(genesSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List>& >::type matchedgenes(matchedgenesSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type temporal(temporalSEXP);
    Rcpp::traits::input_parameter< Nullable<Rcpp::List> >::type targetCounts(targetCountsSEXP);
    rcpp_result_gen = Rcpp::wrap(getGenePrababilities_measurements(targetGene, mainParameters, genes, matchedgenes, temporal, targetCounts));
    return rcpp_result_gen;
END_RCPP
}
// buildProbabilityTreeOnTargetGene
Rcpp::List buildProbabilityTreeOnTargetGene(Rcpp::CharacterVector& targetGene, Rcpp::Environment& mainParameters, Rcpp::CharacterVector& genes, Rcpp::Nullable<Rcpp::List> matchedgenes, Rcpp::Nullable<Rcpp::CharacterVector> matchedexpression, Rcpp::IntegerVector maxK, Rcpp::IntegerVector temporal, Nullable<Rcpp::List> targetCounts, bool findPositiveRegulate, bool findNegativeRegulate);
RcppExport SEXP _FBNNet_buildProbabilityTreeOnTargetGene(SEXP targetGeneSEXP, SEXP mainParametersSEXP, SEXP genesSEXP, SEXP matchedgenesSEXP, SEXP matchedexpressionSEXP, SEXP maxKSEXP, SEXP temporalSEXP, SEXP targetCountsSEXP, SEXP findPositiveRegulateSEXP, SEXP findNegativeRegulateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::CharacterVector& >::type targetGene(targetGeneSEXP);
    Rcpp::traits::input_parameter< Rcpp::Environment& >::type mainParameters(mainParametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector& >::type genes(genesSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type matchedgenes(matchedgenesSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::CharacterVector> >::type matchedexpression(matchedexpressionSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type maxK(maxKSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type temporal(temporalSEXP);
    Rcpp::traits::input_parameter< Nullable<Rcpp::List> >::type targetCounts(targetCountsSEXP);
    Rcpp::traits::input_parameter< bool >::type findPositiveRegulate(findPositiveRegulateSEXP);
    Rcpp::traits::input_parameter< bool >::type findNegativeRegulate(findNegativeRegulateSEXP);
    rcpp_result_gen = Rcpp::wrap(buildProbabilityTreeOnTargetGene(targetGene, mainParameters, genes, matchedgenes, matchedexpression, maxK, temporal, targetCounts, findPositiveRegulate, findNegativeRegulate));
    return rcpp_result_gen;
END_RCPP
}
// process_cube_algorithm
Rcpp::List process_cube_algorithm(Rcpp::CharacterVector& target_gene, Rcpp::CharacterVector& conditional_genes, Rcpp::IntegerVector maxK, Rcpp::IntegerVector temporal, Rcpp::Environment& mainParameters);
RcppExport SEXP _FBNNet_process_cube_algorithm(SEXP target_geneSEXP, SEXP conditional_genesSEXP, SEXP maxKSEXP, SEXP temporalSEXP, SEXP mainParametersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::CharacterVector& >::type target_gene(target_geneSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector& >::type conditional_genes(conditional_genesSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type maxK(maxKSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type temporal(temporalSEXP);
    Rcpp::traits::input_parameter< Rcpp::Environment& >::type mainParameters(mainParametersSEXP);
    rcpp_result_gen = Rcpp::wrap(process_cube_algorithm(target_gene, conditional_genes, maxK, temporal, mainParameters));
    return rcpp_result_gen;
END_RCPP
}
// splitExpression
Rcpp::CharacterVector splitExpression(Rcpp::CharacterVector& expression, int outputType, bool lowerCase);
RcppExport SEXP _FBNNet_splitExpression(SEXP expressionSEXP, SEXP outputTypeSEXP, SEXP lowerCaseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::CharacterVector& >::type expression(expressionSEXP);
    Rcpp::traits::input_parameter< int >::type outputType(outputTypeSEXP);
    Rcpp::traits::input_parameter< bool >::type lowerCase(lowerCaseSEXP);
    rcpp_result_gen = Rcpp::wrap(splitExpression(expression, outputType, lowerCase));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FBNNet_extractGeneStateFromTimeSeriesCube", (DL_FUNC) &_FBNNet_extractGeneStateFromTimeSeriesCube, 2},
    {"_FBNNet_getGenePrababilities_basic", (DL_FUNC) &_FBNNet_getGenePrababilities_basic, 6},
    {"_FBNNet_getGenePrababilities_advanced", (DL_FUNC) &_FBNNet_getGenePrababilities_advanced, 1},
    {"_FBNNet_getGenePrababilities", (DL_FUNC) &_FBNNet_getGenePrababilities, 6},
    {"_FBNNet_networkFiltering", (DL_FUNC) &_FBNNet_networkFiltering, 1},
    {"_FBNNet_getGenePrababilities_measurements", (DL_FUNC) &_FBNNet_getGenePrababilities_measurements, 6},
    {"_FBNNet_buildProbabilityTreeOnTargetGene", (DL_FUNC) &_FBNNet_buildProbabilityTreeOnTargetGene, 10},
    {"_FBNNet_process_cube_algorithm", (DL_FUNC) &_FBNNet_process_cube_algorithm, 5},
    {"_FBNNet_splitExpression", (DL_FUNC) &_FBNNet_splitExpression, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_FBNNet(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
