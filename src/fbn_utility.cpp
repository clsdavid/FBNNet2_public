#include "fbn_utility.h"
#include <Rcpp.h>
using namespace Rcpp;

/////////// external concatenate methods
// At the C-level, all R objects are stored in a common datatype, the SEXP, or S-expression. All R objects are S-expressions so every C function that you create must return a SEXP as output and take SEXPs as inputs. (Technically, this is a pointer to a structure with typedef SEXPREC.) A SEXP is a variant type, with subtypes for all R’s data structures. The most important types are:
//   
//   REALSXP: numeric vector
//   INTSXP: integer vector
//   LGLSXP: logical vector
//   STRSXP: character vector
//   VECSXP: list
//   CLOSXP: function (closure)
//   ENVSXP: environment

// [[Rcpp::export]]
std::string to_string(double val){
  std::ostringstream stm ;
  stm << val ;
  return stm.str() ;
}

// [[Rcpp::export]]
Rcpp::String mpaste(Rcpp::CharacterVector x, std::string sep){
  Rcpp::String res;
  for(int i=0;i<x.length();i++){
    if(i<x.length()-1){
      
      res += String((std::string)x[i]+sep);
    }else{
      res += String(x[i]);
    }
  }
  return (res);
}

// [[Rcpp::export]]
Rcpp::CharacterVector concatenator(Rcpp::CharacterVector a, Rcpp::CharacterVector b) {
  int alen = a.length();
  int blen = b.length();
  Rcpp::CharacterVector out(alen+blen);
  
  for (int j = 0; j < alen+blen; j++) {
    if (j < alen) {
      out(j) = a(j);
    } else {
      out(j) = b(j - alen);
    }
  }
  return(out);
}

// [[Rcpp::export]]
Rcpp::IntegerVector concatenatorI(Rcpp::IntegerVector a, Rcpp::IntegerVector b) {
  int alen = a.length();
  int blen = b.length();
  Rcpp::IntegerVector out(alen+blen);
  
  for (int j = 0; j < alen+blen; j++) {
    if (j < alen) {
      out(j) = a(j);
    } else {
      out(j) = b(j - alen);
    }
  }
  return(out);
}

// [[Rcpp::export]]
Rcpp::NumericVector concatenatorN(Rcpp::NumericVector a, Rcpp::NumericVector b) {
  int alen = a.length();
  int blen = b.length();
  Rcpp::NumericVector out(alen+blen);
  
  for (int j = 0; j < alen+blen; j++) {
    if (j < alen) {
      out(j) = a(j);
    } else {
      out(j) = b(j - alen);
    }
  }
  return(out);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix mcbind(Rcpp::NumericMatrix a, Rcpp::NumericMatrix b) {
  int acoln = a.ncol();
  int bcoln = b.ncol();
  int arown = a.nrow();
  int brown = b.nrow();
  if(arown != brown){
    Rcpp::String msg("The two matrixes must have the same number of rows: nrow(a)=");
    msg += (Rcpp::String)acoln;
    msg += ", ncol(b)=";
    msg += (Rcpp::String)bcoln;
    throw std::runtime_error(msg);
  }
  
  Rcpp::NumericMatrix out = Rcpp::NumericMatrix(a.nrow(), acoln + bcoln);
  for (int j = 0; j < acoln + bcoln; j++) {
    if (j < acoln) {
      out(_, j) = a(_, j);
    } else {
      out(_, j) = b(_, j - acoln);
    }
  }
  rownames(out) = rownames(a);
  if(!Rf_isNull(colnames(a))&&!Rf_isNull(colnames(b)))
  {
    Rcpp::CharacterVector new_colName =  concatenator(colnames(a),colnames(b));
    colnames(out) = new_colName;
  }
  
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix mrbind(Rcpp::NumericMatrix a, Rcpp::NumericMatrix b) {
  int acoln = a.ncol();
  int bcoln = b.ncol();
  int arown = a.nrow();
  int brown = b.nrow();
  if(acoln != bcoln){
    
    Rcpp::String msg("The two matrixes must have the same number of cols: ncol(a)=");
    msg += (Rcpp::String)acoln;
    msg += ", ncol(b)=";
    msg += (Rcpp::String)bcoln;
    throw std::runtime_error(msg);
  }
  
  Rcpp::NumericMatrix out = Rcpp::NumericMatrix(arown + brown, a.ncol());
  for (int j = 0; j < arown + brown; j++) {
    if (j < arown) {
      out(j,_) = a(j,_);
    } else {
      out(j,_) = b(j - arown,_);
    }
  }
  colnames(out) = colnames(a);
  if(!Rf_isNull(rownames(a))&&!Rf_isNull(rownames(b)))
  {
    Rcpp::CharacterVector new_rowName =  concatenator(rownames(a),rownames(b));
    rownames(out) = new_rowName;
  }
  
  return out;
}

// [[Rcpp::export]]
bool isReallyNA(double val) {
  NumericVector val2 = NumericVector::create(val);
  return Rcpp::all(Rcpp::is_na(val2));
}


// [[Rcpp::export]]
int countZeros(Rcpp::NumericVector v){
  int c = 0;
  for(int i=0; i<v.length(); ++i){
    if(v[i]==0)c++;
  }
  return(c);
}


//' @title Accessing R's fisher.test function from Rcpp
// [[Rcpp::export]]
Rcpp::List fisher_test_cpp(const Rcpp::NumericMatrix& x, double conf_level){
  
  // Obtain environment containing function
  Rcpp::Environment base("package:stats"); 
  
  // Make function callable from C++
  Rcpp::Function fisher_test = base["fisher.test"];    
  
  // Call the function and receive its list output
  Rcpp::List test_out = fisher_test(Rcpp::_["x"] = x,
                                    Rcpp::_["conf.level"]  = conf_level);
  
  // Return test object in list structure
  return test_out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix substractM(Rcpp::NumericMatrix m, Rcpp::NumericVector v)
{
  Rcpp::NumericMatrix res =NumericMatrix(m.nrow(),m.ncol());// clone(m);
  const int ncol=m.ncol();
  for(int i=0; i<ncol; ++i){
    res(_,i)= abs(m(_,i)-v);//abs is Rcpp type abs
  }
  return(res);
}

// [[Rcpp::export]]
int matchCount(Rcpp::NumericMatrix m, Rcpp::NumericVector v){
  // Rcpp::NumericMatrix mid=substractM(m,v);
  // return(countZeros(Rcpp::colSums(mid)));
  Rcpp::NumericMatrix mid=substractM(m,v);
  Rcpp::NumericVector col_sumed = Rcpp::colSums(mid);
  return(std::count(col_sumed.begin(),col_sumed.end(),0));
}

// [[Rcpp::export]]
double dround(double val, int decimal){
  double dec = pow(10.0,decimal);
  return(round( val * dec) / dec);
}

// [[Rcpp::export]]
Rcpp::LogicalVector a_in_b(Rcpp::CharacterVector names1, Rcpp::CharacterVector names2){
  Rcpp::LogicalVector res(names1.length());
  for(int i=0;i<names1.length();i++){
    for(int j=0;j<names2.length();j++){
      res[i]=false;
      if(names1[i]==names2[j]){
        res[i]=true;
        break;
      }
    }
  }
  return(res);
}

// [[Rcpp::export]]
Rcpp::IntegerVector a_in_b_index(Rcpp::CharacterVector names1, Rcpp::CharacterVector names2){
  Rcpp::IntegerVector res;
  for(int i=0;i<names1.length();i++){
    for(int j=0;j<names2.length();j++){
      if(names1[i]==names2[j]){
        res.push_back(j);
        break;
      }
    }
  }
  return(res);
}

// [[Rcpp::export]]
Rcpp::LogicalVector a_not_in_b(Rcpp::CharacterVector names1, Rcpp::CharacterVector names2){
  Rcpp::LogicalVector res(names1.length());
  for(int i=0;i<names1.length();i++){
    for(int j=0;j<names2.length();j++){
      res[i]=true;
      if(names1[i]==names2[j]){
        res[i]=false;
        break;
      }
    }
  }
  return(res);
}

// [[Rcpp::export]]
Rcpp::IntegerVector a_not_in_b_index(Rcpp::CharacterVector names1, Rcpp::CharacterVector names2){
  Rcpp::IntegerVector res;
  for(int i=0;i<names1.length();i++){
    bool match = false;
    for(int j=0;j<names2.length();j++){
      if(names1[i]==names2[j]){
        match = true;
        break;
      }
    }
    if(!match){
      res.push_back(i);
    }
  }
  return(res);
}

// [[Rcpp::export]]
Rcpp::List resizel( const Rcpp::List& x, int n ){
  int oldsize = x.size() ;
  List y(n) ;
  for( int i=0; i<oldsize; i++) y[i] = x[i] ;
  return y ;
}

// [[Rcpp::export]]
Rcpp::List orderByname( const Rcpp::List& x, Rcpp::CharacterVector names){
  int len = names.length() ;
  List y(len) ;
  for( int i=0; i<len; i++) {
    std::string name = (std::string)names[i];
    y[i] = x[name] ;
  }
  y.attr("names") = Rcpp::wrap(names);
  return y ;
}

// [[Rcpp::export]]
Rcpp::List removeEmptyElement(Rcpp::List x)
{
  IntegerVector idx;
  for(int i=0;i<x.length();i++){
    if(!Rf_isNull(x[i])){
      idx.push_back(i);
    }
  }
  return(x[idx]);
}

// [[Rcpp::export]]
Rcpp::CharacterVector char_sort(Rcpp::CharacterVector x, bool dsc) {
  Rcpp::CharacterVector res = Rcpp::clone(x);
  res.sort(dsc);
  return res;
}

// [[Rcpp::export]]
Rcpp::IntegerVector int_sort(Rcpp::IntegerVector x, bool dsc) {
  Rcpp::IntegerVector res = Rcpp::clone(x);
  res.sort(dsc);
  return res;
}

// [[Rcpp::export]]
Rcpp::NumericVector num_sort(Rcpp::NumericVector x, bool dsc) {
  Rcpp::NumericVector res = Rcpp::clone(x);
  res.sort(dsc);
  return res;
}

// [[Rcpp::export]]
Rcpp::CharacterVector convertStringIntoVector(std::string value, int outputType, bool lowerCase) {
  int len = value.length();
  Rcpp::CharacterVector res;
  std::string t;
  for(int i=0;i<len;i++) {
    char v = value[i];
    if(lowerCase) {
      v = std::tolower(v);
    }
    if(v == ' ') {
      if(t.length()>0) {
        res.push_back(t);
        t.clear();
      }
      continue;
    }
    
    if(outputType == 1 || outputType ==0) {
      if(v == '&' || v == ',' || v == '!') {
        if(t.length()>0) {
          res.push_back(t);
          t.clear();
        }
        
        if(outputType==1) {
          t.push_back(v);
        } else {
          if(v!='!') {
            t.push_back(v);;
          }
        }
        res.push_back(t);
        t.clear();
        continue;
      }
      
      if(outputType==1) {
        t.push_back(v);
      } else {
        if(v!='!') {
          t.push_back(v);;
        }
      }
    } else {
      if(v == '&' || v == ',') {
        if(t.length()>0) {
          res.push_back(t);
          t.clear();
        }
        t.push_back(v);
        res.push_back(t);
        t.clear();
        continue;
      }
      t.push_back(v);
    }
  }
  if(t.length()>0) {
    res.push_back(t);
  }
  
  return(res);
}


// [[Rcpp::export]]
Rcpp::CharacterVector subCPP(Rcpp::CharacterVector pattern, Rcpp::CharacterVector replacement, Rcpp::CharacterVector x) {
  int len = x.size();
  CharacterVector y(len);
  int patlen = pattern.size();
  int replen = replacement.size();
  if (patlen != replen)
    Rcout<<"Error: Pattern and replacement length do not match";
  for(int i = 0; i < patlen; ++i) {
    if (*(char*)x[i] == *(char*)pattern[i])
      y[x[i] == pattern[i]] = replacement[i];
  }
  return y;
} 

// [[Rcpp::export]]
Rcpp::CharacterVector splitExpression(Rcpp::CharacterVector expression, int outputType, bool lowerCase) {
  std::string strExpression = (std::string)expression[0];

  Rcpp::CharacterVector res;
  if (strExpression =="1" || strExpression == "0") {
    res.push_back(strExpression);
  } else{
    if(outputType == 1) {
      res = convertStringIntoVector(strExpression, 1, lowerCase);
    } else {
      
      res = convertStringIntoVector(strExpression, 2, lowerCase);
    }
  }
  return(res);
}