#include "fbn_utility.h"
#include <Rcpp.h>
using namespace Rcpp;

/////////// external concatenate methods
// At the C-level, all R objects are stored in a common datatype, the SEXP,
// or S-expression. All R objects are S-expressions so every C function that you create
// must return a SEXP as output and take SEXPs as inputs. (Technically, this is a pointer
// to a structure with typedef SEXPREC.) A SEXP is a variant type, with subtypes for all
// R’s data structures. The most important types are:
//
//   REALSXP: numeric vector
//   INTSXP: integer vector
//   LGLSXP: logical vector
//   STRSXP: character vector
//   VECSXP: list
//   CLOSXP: function (closure)
//   ENVSXP: environment

std::string to_string(double val){
   std::ostringstream stm ;
   stm << val ;
   return stm.str() ;
}

Rcpp::String mpaste(Rcpp::CharacterVector& x,
                    std::string sep){
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

// concatenate two char vectors
Rcpp::CharacterVector concatenator(Rcpp::CharacterVector& a,
                                   Rcpp::CharacterVector& b) {
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

// concatenate two integer vectors
Rcpp::IntegerVector concatenatorI(Rcpp::IntegerVector& a,
                                  Rcpp::IntegerVector& b) {
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

// concatenate two numeric vectors
Rcpp::NumericVector concatenatorN(Rcpp::NumericVector& a,
                                  Rcpp::NumericVector& b) {
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

// concatenate two numeric matrix by columns
Rcpp::NumericMatrix mcbind(Rcpp::NumericMatrix& a,
                           Rcpp::NumericMatrix& b) {
   int acoln = a.ncol();
   int bcoln = b.ncol();
   int arown = a.nrow();
   int brown = b.nrow();
   CharacterVector a_names = colnames(a);
   CharacterVector b_names = colnames(b);
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
   if(!Rf_isNull(a_names)&&!Rf_isNull(b_names))
   {
      Rcpp::CharacterVector new_colName =  concatenator(a_names, b_names);
      colnames(out) = new_colName;
   }

   return out;
}

// concatenate two numeric matrix by rows
Rcpp::NumericMatrix mrbind(Rcpp::NumericMatrix& a,
                           Rcpp::NumericMatrix& b) {
   int acoln = a.ncol();
   int bcoln = b.ncol();
   int arown = a.nrow();
   int brown = b.nrow();
   CharacterVector a_names = rownames(a);
   CharacterVector b_names = rownames(b);
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
      Rcpp::CharacterVector new_rowName =  concatenator(a_names, b_names);
      rownames(out) = new_rowName;
   }

   return out;
}

bool isReallyNA(double val) {
  NumericVector val2 = NumericVector::create(val);
  return Rcpp::all(Rcpp::is_na(val2));
}


int countZeros(Rcpp::NumericVector& v){
   int c = 0;
   for(int i=0; i<v.length(); ++i){
      if(v[i]==0)c++;
   }
   return(c);
}


Rcpp::List fisher_test_cpp(const Rcpp::NumericMatrix& x,
                           double conf_level){

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

Rcpp::NumericMatrix substractM(Rcpp::NumericMatrix& m,
                               Rcpp::NumericVector& v)
{
   Rcpp::NumericMatrix res =NumericMatrix(m.nrow(),m.ncol());// clone(m);
   const int ncol=m.ncol();
   for(int i=0; i<ncol; ++i){
      res(_,i)= abs(m(_,i)-v);//abs is Rcpp type abs
   }
   return(res);
}

int matchCount(Rcpp::NumericMatrix& m, Rcpp::NumericVector& v){
   Rcpp::NumericMatrix mid=substractM(m,v);
   Rcpp::NumericVector col_sumed = Rcpp::colSums(mid);
   return(std::count(col_sumed.begin(),col_sumed.end(),0));
}

double dround(double val, int decimal){
  double dec = pow(10.0,decimal);
  return(round( val * dec) / dec);
}

Rcpp::LogicalVector a_in_b(Rcpp::CharacterVector& names1, Rcpp::CharacterVector& names2){
   Rcpp::LogicalVector res(names1.length());
   if(names2.length() == 0) {
      return(Rcpp::LogicalVector(0));
   }


   if(names1.length() == 0) {
      return(res);
   }
   Rcpp::IntegerVector v = Rcpp::seq(0, names2.size()-1);
   Rcpp::IntegerVector v_matched = match(names1,names2);
   res[!Rcpp::is_na(v_matched)] = true;
   return(res);
}



Rcpp::IntegerVector a_in_b_index(Rcpp::CharacterVector& names1, Rcpp::CharacterVector& names2){
   if(names2.length() == 0) {
      return(Rcpp::IntegerVector(0));
   }

   Rcpp::IntegerVector v = Rcpp::seq(0, names2.size()-1);
   if(names1.length() == 0) {
      return(v);
   }

   Rcpp::IntegerVector v_matched = match(names2,names1);
   return(v[!Rcpp::is_na(v_matched)]);
}

Rcpp::LogicalVector a_not_in_b(Rcpp::CharacterVector& names1, Rcpp::CharacterVector& names2){
   Rcpp::LogicalVector res(names1.length());
   if(names2.length() == 0) {
      return(Rcpp::LogicalVector(0));
   }


   if(names1.length() == 0) {
      return(res);
   }
   Rcpp::IntegerVector v = Rcpp::seq(0, names2.size()-1);
   Rcpp::IntegerVector v_matched = match(names1,names2);
   res[Rcpp::is_na(v_matched)] = true;
   return(res);
}


Rcpp::IntegerVector a_not_in_b_index(Rcpp::CharacterVector& names1, Rcpp::CharacterVector& names2){
   if(names1.length() == 0) {
      return(Rcpp::IntegerVector(0));
   }

   Rcpp::IntegerVector v = Rcpp::seq(0, names1.size()-1);
   if(names2.length() == 0) {
      return(v);
   }

   Rcpp::IntegerVector v_matched = match(names1, names2);
   return(v[Rcpp::is_na(v_matched)]);
}

Rcpp::List resizel( const Rcpp::List& x, int n ){
   int oldsize = x.size() ;
   List y(n) ;
   for( int i=0; i<oldsize; i++) y[i] = x[i] ;
   return y ;
}


Rcpp::List orderByname( const Rcpp::List& x, Rcpp::CharacterVector& names){
   int len = names.length() ;
   List y(len) ;
   for( int i=0; i<len; i++) {
      std::string name = (std::string)names[i];
      y[i] = x[name] ;
   }
   y.attr("names") = Rcpp::wrap(names);
   return y ;
}


Rcpp::List removeEmptyElement(Rcpp::List& x)
{
   IntegerVector idx;
   for(int i=0;i<x.length();i++){
      if(!Rf_isNull(x[i])){
         idx.push_back(i);
      }
   }
   return(x[idx]);
}


Rcpp::CharacterVector char_sort(Rcpp::CharacterVector& x, bool dsc) {
   Rcpp::CharacterVector res = Rcpp::clone(x);
   res.sort(dsc);
   return res;
}


Rcpp::IntegerVector int_sort(Rcpp::IntegerVector& x, bool dsc) {
   Rcpp::IntegerVector res = Rcpp::clone(x);
   res.sort(dsc);
   return res;
}


Rcpp::NumericVector num_sort(Rcpp::NumericVector& x, bool dsc) {
    Rcpp::NumericVector res = Rcpp::clone(x);
    res.sort(dsc);
    return res;
}


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



Rcpp::CharacterVector subCPP(Rcpp::CharacterVector& pattern, Rcpp::CharacterVector& replacement, Rcpp::CharacterVector& x) {
   CharacterVector y(x.size());
   int patlen = pattern.size();
   int replen = replacement.size();
   if (patlen != replen)
      Rcout<<"Error: Pattern and replacement length do not match";
   for(int i = 0; i < patlen; ++i) {
      if (*(char*)x[i] == *(char*) pattern[i])
      y[x[i] == pattern[i]] = replacement[i];
   }
   return y;
}

//' A function to split an expression into a vector of input
//' @name splitExpression
//' @param expression The expression fo a FBN network connection
//' @param outputType The type of output.
//' @param lowerCase Optional, if TRUE convert them to lower case.
// [[Rcpp::export]]
Rcpp::CharacterVector splitExpression(Rcpp::CharacterVector& expression,
                                      int outputType,
                                      bool lowerCase) {
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
