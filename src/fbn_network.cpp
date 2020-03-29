#include <Rcpp.h>
#include "fbn_utility.h"
using namespace Rcpp;

//' A function to filter networks.
//' @name networkFiltering
//' @param res A list of named network interactions
// [[Rcpp::export]]
Rcpp::List networkFiltering(Rcpp::List res) {
  
  //finalFilteredlist <- list()
  Rcpp::CharacterVector targetgenes = res.attr("names");
  int len = targetgenes.length();
  Rcpp::List filteredres(len);
  Rcpp::List filtered =  Rcpp::clone(res);
  for(int i=0;i<len;i++) {
    std::string target = (std::string)targetgenes[i];
    Rcpp::List sub_rule_set = res[target];
    Rcpp::IntegerVector filteredIndex;
    //filteredres = sub_rule_set;
    Rcpp::IntegerVector indexes;
    int cond_len = sub_rule_set.length();
    for(int j=0;j<cond_len;j++) {
      indexes.push_back(j);
      if(Rf_isNull(sub_rule_set[j])) {
        filteredIndex.push_back(j);
        continue;
      }
      
      Rcpp::CharacterVector rule = sub_rule_set[j];
      for(int k=0;k<cond_len;k++) {
        if(Rf_isNull(sub_rule_set[k])) {
          filteredIndex.push_back(k);
          continue;
        }
        
        Rcpp::CharacterVector rule2 = sub_rule_set[k];
        std::string str_numOfInput_rule1 = (std::string)rule["numOfInput"];
        std::string str_numOfInput_rule2 = (std::string)rule2["numOfInput"];

        int numOfInput_rule1 = std::atoi(str_numOfInput_rule1.c_str());
        int numOfInput_rule2 = std::atoi(str_numOfInput_rule2.c_str());
        
        std::string str_type_rule1 = (std::string)rule["type"];
        std::string str_type_rule2 = (std::string)rule2["type"];

        std::string str_timestep_rule1 = (std::string)rule["timestep"];
        std::string str_timestep_rule2 = (std::string)rule2["timestep"];

        std::string str_input_rule1 = (std::string)rule["input"];
        std::string str_input_rule2 = (std::string)rule2["input"];

        Rcpp::CharacterVector input_rule1 = splitExpression(str_input_rule1, 2, false);
        Rcpp::CharacterVector input_rule2 = splitExpression(str_input_rule2, 2, false);
        if (is_true(all(a_in_b(input_rule1,input_rule2))))
        {
          if(numOfInput_rule1< numOfInput_rule2 && str_type_rule1 == str_type_rule2 && str_timestep_rule1 == str_timestep_rule2) {
            filteredIndex.push_back(k);
          }
        }
      }
    }
    
    Rcpp::IntegerVector diffset = setdiff(indexes, filteredIndex);
    filtered[target] = sub_rule_set[diffset];
  }
  
  return(filtered);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

