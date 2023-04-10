#ifndef FBN_UTILITY_H
#define FBN_UTILITY_H

#include <Rcpp.h>

std::string to_string(double val);
Rcpp::String mpaste(Rcpp::CharacterVector& x, std::string sep = "");
Rcpp::CharacterVector concatenator(Rcpp::CharacterVector& a, Rcpp::CharacterVector& b);
Rcpp::IntegerVector concatenatorI(Rcpp::IntegerVector& a, Rcpp::IntegerVector& b);
Rcpp::NumericVector concatenatorN(Rcpp::NumericVector& a, Rcpp::NumericVector& b);
Rcpp::NumericMatrix mrbind(Rcpp::NumericMatrix& a, Rcpp::NumericMatrix& b);
Rcpp::NumericMatrix mcbind(Rcpp::NumericMatrix& a, Rcpp::NumericMatrix& b);
bool isReallyNA(double val);
int countZeros(Rcpp::NumericVector& v);
Rcpp::List fisher_test_cpp(const Rcpp::NumericMatrix& x, double conf_level = 0.95);
Rcpp::NumericMatrix substractM(Rcpp::NumericMatrix& m, Rcpp::NumericVector& v);
int matchCount(Rcpp::NumericMatrix& m, Rcpp::NumericVector& v);
double dround(double val, int decimal);
Rcpp::LogicalVector a_in_b(Rcpp::CharacterVector& names1, Rcpp::CharacterVector& names2);
Rcpp::IntegerVector a_in_b_index(Rcpp::CharacterVector& names1, Rcpp::CharacterVector& names2);
Rcpp::List resizel( const Rcpp::List& x, int n );
Rcpp::List orderByname( const Rcpp::List& x, Rcpp::CharacterVector& names);
Rcpp::LogicalVector a_not_in_b(Rcpp::CharacterVector& names1, Rcpp::CharacterVector& names2);
Rcpp::IntegerVector a_not_in_b_index(Rcpp::CharacterVector& names1, Rcpp::CharacterVector& names2);
Rcpp::List removeEmptyElement(Rcpp::List& x);
Rcpp::CharacterVector char_sort(Rcpp::CharacterVector& x, bool dsc=true);
Rcpp::IntegerVector int_sort(Rcpp::IntegerVector& x, bool dsc=true);
Rcpp::NumericVector num_sort(Rcpp::NumericVector& x, bool dsc=true);
Rcpp::CharacterVector convertStringIntoVector(std::string value, int outputType = 1, bool lowerCase = false);
Rcpp::CharacterVector subCPP(Rcpp::CharacterVector& pattern, Rcpp::CharacterVector& replacement, Rcpp::CharacterVector& x);
Rcpp::CharacterVector splitExpression(Rcpp::CharacterVector& expression, int outputType = 1, bool lowerCase = false);
#endif
