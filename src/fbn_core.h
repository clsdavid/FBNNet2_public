#ifndef FBN_CORE_H
#define FBN_CORE_H

#include <Rcpp.h>

Rcpp::NumericMatrix extractGeneStateFromTimeSeriesCube(Rcpp::List& timeSeriesCube, Rcpp::IntegerVector temporal=1);
Rcpp::NumericMatrix extractGeneStates(Rcpp::NumericMatrix& stateMatrix, Rcpp::CharacterVector& targetgenes);
Rcpp::List generateTemporalGeneStates(Rcpp::Environment& mainParameters, Rcpp::CharacterVector& targetgene, Rcpp::CharacterVector& conditional_genes, Rcpp::IntegerVector temporal=1);
Rcpp::List getBasicMeasures(
    Rcpp::NumericVector& stateTCond,
    Rcpp::NumericMatrix& m,
    Rcpp::NumericMatrix& mc,
    Rcpp::NumericVector& cond_T_target_T_state,
    Rcpp::NumericVector& cond_F_target_T_state,
    Rcpp::NumericVector& cond_T_target_F_state,
    Rcpp::NumericVector& cond_F_target_F_state,
    Rcpp::NumericVector& cond_T_target_T_state_c,
    Rcpp::NumericVector& cond_F_target_T_state_c,
    Rcpp::NumericVector& cond_T_target_F_state_c,
    Rcpp::NumericVector& cond_F_target_F_state_c,
    bool recount_target
);
Rcpp::List getAdvancedMeasures(Rcpp::List& basic_measures);
Rcpp::List getGenePrababilities_advanced(Rcpp::List& getGenePrababilities_basic);
Rcpp::List getGenePrababilities_basic(Rcpp::Environment& main_parameters_in_ref, Rcpp::Nullable<Rcpp::List>& fixedgenestate, Rcpp::CharacterVector& target_gene, Rcpp::CharacterVector& new_conditional_gene, Rcpp::IntegerVector temporal=1, Rcpp::Nullable<Rcpp::List> targetCounts = R_NilValue);
Rcpp::List getGenePrababilities(Rcpp::Environment& main_parameters_in_ref, Rcpp::Nullable<Rcpp::List>& fixedgenestate, Rcpp::CharacterVector& target_gene, Rcpp::CharacterVector& new_conditional_gene, Rcpp::IntegerVector temporal=1, Rcpp::Nullable<Rcpp::List> targetCounts = R_NilValue);
#endif
