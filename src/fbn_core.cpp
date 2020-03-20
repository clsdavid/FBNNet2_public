#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <String.h>
#include<limits.h>
#include <math.h>
#include <stdexcept>
#include "fbn_utility.h"
#include "fbn_core.h"
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericMatrix extractGeneStateFromTimeSeriesCube(
    Rcpp::List timeSeriesCube, 
    Rcpp::IntegerVector temporal)
{
  int numOfmats = timeSeriesCube.length();
  if(numOfmats == 1){
    return(timeSeriesCube[0]);
  }
  Rcpp::CharacterVector genes = Rcpp::rownames(timeSeriesCube[0]);
  int thisTemporal = temporal[0];
  Rcpp::NumericMatrix nine_matrix(genes.length(),thisTemporal);
  std::fill(nine_matrix.begin(), nine_matrix.end(), 9);

  int final_len = 0;
  for(int i =0;i<numOfmats;i++){
    NumericMatrix temp = timeSeriesCube[i];
    final_len += temp.ncol();
  }
  
  Rcpp::NumericMatrix res(genes.length(),final_len+thisTemporal*numOfmats-thisTemporal);
  int offset=0;
  for(int i=0;i<numOfmats;i++){
    Rcpp::NumericMatrix cur_mat = timeSeriesCube[i];
    int cur_col_len = cur_mat.ncol();
    
    for(int j=0;j<cur_col_len;j++){
      res(_,offset) = cur_mat(_,j);
      offset +=1;
    }
    if(i<(numOfmats-1))
    {
      for(int j=0;j<thisTemporal;j++){
        res(_,offset) = nine_matrix(_,j);
        offset +=1;
      }
    }

  }
  
  Rcpp::rownames(res) = genes;
  return (res);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix extractGeneStates(Rcpp::NumericMatrix stateMatrix,
                                      Rcpp::CharacterVector targetgenes)
{
  Rcpp::CharacterVector names = Rcpp::rownames(stateMatrix);
  IntegerVector r_index = a_in_b_index(targetgenes,names);
  r_index = int_sort(r_index,false);
  Rcpp::NumericMatrix sub(targetgenes.length(),stateMatrix.ncol());
  for(int i=0; i<r_index.length(); i++)
     sub.row(i) = stateMatrix.row((int)r_index[i]);
  
  rownames(sub) = names[r_index];
  return (sub);
}

// [[Rcpp::export]]
Rcpp::List generateTemporalGeneStates(Rcpp::Environment mainParameters,
                                      Rcpp::CharacterVector targetgene, 
                                      Rcpp::CharacterVector conditional_genes, 
                                      Rcpp::IntegerVector temporal)
{
    Rcpp::List getCurrentStates  =  mainParameters["currentStates"];
    Rcpp::List getpreviousStates  =  mainParameters["previousStates"];
    Rcpp::List getCurrentStates_c  =  mainParameters["currentStates_c"];
    Rcpp::List getpreviousStates_c  =  mainParameters["previousStates_c"];
    int thisTemporal = temporal[0];
    
    if(thisTemporal<1)
        thisTemporal =1;
     
    int cur_len = getCurrentStates.length();
    if(cur_len < thisTemporal || cur_len == 0)
    {
      throw std::runtime_error("getCurrentStates subscript out of bounds");
    }

    int pre_len =  getpreviousStates.length(); 
    if(pre_len < thisTemporal || pre_len == 0)
    {
      throw std::runtime_error("getPreviousStates counter subscript out of bounds");
    }

    int cur_len_c = getCurrentStates_c.length();    
    if(cur_len_c < thisTemporal || cur_len_c == 0)
    {
      throw std::runtime_error("getCurrentStates_c subscript out of bounds");
    }

    int pre_len_c =  getpreviousStates_c.length();      
    if(pre_len_c < thisTemporal || pre_len_c == 0)
    {
      throw std::runtime_error("getPreviousStates counter subscript out of bounds");
    }
    
    std::vector<std::string> names;
    List result(thisTemporal);
    for(int i=0;i<thisTemporal;i++){
      Rcpp::NumericMatrix getCurrentState =  getCurrentStates[i];
      Rcpp::NumericMatrix getpreviousState  =  getpreviousStates[i];
      Rcpp::NumericMatrix getCurrentState_c  =  getCurrentStates_c[i];
      Rcpp::NumericMatrix getpreviousState_c  =  getpreviousStates_c[i];
      int n_state = getpreviousState.ncol();
      
      if((n_state-thisTemporal)<2){
        throw std::runtime_error("No enough states for this temporal");
      }
      n_state = n_state-1;
      Rcpp::NumericMatrix t_getCurrentState = getCurrentState(_, Range(i+1,n_state)); //remove begin
      rownames(t_getCurrentState)= rownames(getCurrentState);
      
      Rcpp::NumericMatrix t_getpreviousState = getpreviousState(_, Range(0, n_state-(i+1)));//remove last part
      rownames(t_getpreviousState)= rownames(getpreviousState);
      
      Rcpp::NumericMatrix t_getCurrentState_c = getCurrentState_c(_, Range(i+1,n_state)); //remove begin
      rownames(t_getCurrentState_c)= rownames(getCurrentState_c);
      
      Rcpp::NumericMatrix t_getpreviousState_c = getpreviousState_c(_, Range(0, n_state-(i+1)));//remove last part
      rownames(t_getpreviousState_c)= rownames(getpreviousState_c);
      
      
      Rcpp::NumericMatrix extracedCondition = extractGeneStates(t_getpreviousState, conditional_genes);
      Rcpp::NumericMatrix extracedTargetGeneState = extractGeneStates(t_getCurrentState,targetgene);
        
      Rcpp::NumericMatrix extracedCondition_c = extractGeneStates(t_getpreviousState_c, targetgene);
      Rcpp::NumericMatrix extracedTargetGeneState_c = extractGeneStates(t_getCurrentState_c, conditional_genes);
        
      Rcpp::NumericMatrix concatenatedMatrix = mrbind(extracedCondition,extracedTargetGeneState);
      Rcpp::NumericMatrix concatenatedMatrix_c = mrbind(extracedCondition_c,extracedTargetGeneState_c);
        
      std::vector<std::string> subnames(3);
      Rcpp::List subresult(3);
      subresult[0] = concatenatedMatrix; 
      subnames[0] = "computation_Matrix";
      subresult[1]= concatenatedMatrix_c; 
      subnames[1] = "computation_Matrix_c";
      subresult[2]= i+1; 
      subnames[2] = "timeStep";
      subresult.attr("names") = Rcpp::wrap(subnames); 
        
      result[i]=subresult;
    }
   return(result);
}

// [[Rcpp::export]]
Rcpp::List getBasicMeasures(
    Rcpp::NumericVector stateTCond,
    Rcpp::NumericMatrix m,
    Rcpp::NumericMatrix mc,
    Rcpp::NumericVector cond_T_target_T_state,
    Rcpp::NumericVector cond_F_target_T_state,
    Rcpp::NumericVector cond_T_target_F_state,
    Rcpp::NumericVector cond_F_target_F_state,
    Rcpp::NumericVector cond_T_target_T_state_c,
    Rcpp::NumericVector cond_F_target_T_state_c,
    Rcpp::NumericVector cond_T_target_F_state_c,
    Rcpp::NumericVector cond_F_target_F_state_c,
    bool recount_target
  ) {
  int target_T_count=0;
  int target_F_count=0;

//condition#####################################################################################
  int lenTT = matchCount(m, cond_T_target_T_state);
  int lenTF = matchCount(m, cond_F_target_T_state);
  int lenFT = matchCount(m, cond_T_target_F_state);
  int lenFF = matchCount(m, cond_F_target_F_state);
  if(recount_target){
    if(stateTCond.length()>1){
      Rcpp::NumericMatrix m2 = m(Range(m.nrow()-2,m.nrow()-1),_);
      int sresTT=matchCount(m2, NumericVector::create(1,1));
      int sresTF=matchCount(m2, NumericVector::create(0,1));
      int sresFT=matchCount(m2, NumericVector::create(1,0));
      int sresFF=matchCount(m2, NumericVector::create(0,0));
      
      target_T_count = sresTT+sresTF;
      target_F_count = sresFT+sresFF;
    }else{
      target_T_count = lenTT + lenTF;
      target_F_count = lenFT + lenFF;
    }
  }

  int lenTT_c = matchCount(mc,cond_T_target_T_state_c);
  int lenTF_c = matchCount(mc,cond_F_target_T_state_c);
  int lenFT_c = matchCount(mc,cond_T_target_F_state_c);
  int lenFF_c = matchCount(mc,cond_F_target_F_state_c);

  List res = List::create(_("target_T_count") = target_T_count,
                          _("target_F_count") = target_F_count,
                          _("cond_T_count") = lenTT+lenFT,
                          _("cond_F_count") = lenTF+lenFF,
                          _("cond_T_count_c") = lenTT_c + lenFT_c,
                          _("cond_F_count_c") = lenTF_c + lenFF_c,
                          _("lenTT") = lenTT,
                          _("lenTF") = lenTF,
                          _("lenFT") = lenFT,
                          _("lenFF") = lenFF,
                          _("lenTT_c") = lenTT_c,
                          _("lenTF_c") = lenTF_c,
                          _("lenFT_c") = lenFT_c,
                          _("lenFF_c") = lenFF_c);
  return res;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

// [[Rcpp::export]]
Rcpp::List getGenePrababilities_basic(Rcpp::Environment main_parameters_in_ref,
                                      Rcpp::Nullable<Rcpp::List> fixedgenestate, 
                                      Rcpp::CharacterVector target_gene, 
                                      Rcpp::CharacterVector new_conditional_gene,
                                      Rcpp::IntegerVector temporal,
                                      Rcpp::Nullable<Rcpp::List> targetCounts)
{
  
  //P(B|A)=P(A & B)/P(A)=(frq(A&B)/n)/(frq(A)/n)=frq(A&B)/frq(A)=P(A&D)/(P(A&D)+P(!A & D)) Bayers rule
  //order of target_genes is very important
  int total_samples = main_parameters_in_ref["total_samples"];//ToDo, should pass from up stream as it is static
  CharacterVector all_gene_names = main_parameters_in_ref["all_gene_names"];
  int n_timepoints = main_parameters_in_ref["total_timepoints"];;
  CharacterVector conditional_genes;
  
  List cur_fixed_state;
  
  if(fixedgenestate.isNull())
  {
      conditional_genes = new_conditional_gene;
  } else
  {
      Rcpp::List temp_state(fixedgenestate);
      cur_fixed_state = temp_state;
      //sub level
      conditional_genes = cur_fixed_state.names();
    
    
      if(!is_true(all(a_in_b(conditional_genes,all_gene_names))))
      {
          throw std::runtime_error("All or some part of the conditional genes are not founded in the timeseries cube");
      }
    
      if(is_true(any(a_in_b(new_conditional_gene,conditional_genes))))
      {
          IntegerVector indexs = a_in_b_index(new_conditional_gene, conditional_genes);
          for(int i=0;i<indexs.length();i++){
              conditional_genes.erase(i);
          }
        
          IntegerVector indexs2 = a_in_b_index(new_conditional_gene, cur_fixed_state.names());
          for(int i=0;i<indexs2.length();i++){
              cur_fixed_state.erase(i);
          }
      }else
      {
          conditional_genes = concatenator(conditional_genes, new_conditional_gene);
      }
    
  }
  
  List cond_gene_T_states = cur_fixed_state;
  List cond_gene_F_states = cur_fixed_state;
  
  int s_len = cond_gene_T_states.length();
  cond_gene_T_states = resizel(cond_gene_T_states,s_len+1);
  cond_gene_F_states = resizel(cond_gene_F_states,s_len+1);
  
  //#Add the new conditional gene to the current fixedgenestate
  cond_gene_T_states[s_len]=1;
  cond_gene_F_states[s_len]=0;
  cond_gene_T_states.attr("names") = Rcpp::wrap(conditional_genes);
  cond_gene_F_states.attr("names") = Rcpp::wrap(conditional_genes);
  
  //#get states in order, the order is very important
  IntegerVector indexs3 = a_in_b_index(conditional_genes,all_gene_names).sort();
  CharacterVector uniqued_conditional_genes = all_gene_names[indexs3];
  
  int num_of_conditional_genes = uniqued_conditional_genes.length();
  
  cond_gene_T_states = orderByname(cond_gene_T_states, uniqued_conditional_genes);
  cond_gene_F_states = orderByname(cond_gene_F_states, uniqued_conditional_genes);
  

  NumericVector stateTCond;
  NumericVector stateFCond;
  for(int i=0;i < cond_gene_T_states.length(); i++){
      stateTCond.push_back(cond_gene_T_states[i], (std::string)uniqued_conditional_genes[i]);
  }
  
  for(int i=0;i < cond_gene_F_states.length(); i++){
      stateFCond.push_back(cond_gene_F_states[i], (std::string)uniqued_conditional_genes[i]);
  }
  
  NumericVector mTRUE = NumericVector::create(1);
  NumericVector mFALSE = NumericVector::create(0);
  
  //#prepare, the last bit is the target, the variable name TT, the front is target, the other is conditional
  NumericVector cond_T_target_T_state = concatenatorN(stateTCond, mTRUE);
  NumericVector cond_F_target_T_state = concatenatorN(stateFCond, mTRUE);
  NumericVector cond_T_target_F_state = concatenatorN(stateTCond, mFALSE);
  NumericVector cond_F_target_F_state = concatenatorN(stateFCond, mFALSE);
  //#prepare counter
  NumericVector cond_T_target_T_state_c = concatenatorN(mTRUE, stateTCond);
  NumericVector cond_F_target_T_state_c = concatenatorN(mFALSE, stateTCond);
  NumericVector cond_T_target_F_state_c = concatenatorN(mTRUE, stateFCond);
  NumericVector cond_F_target_F_state_c = concatenatorN(mFALSE, stateFCond);
  
  
  //#get all combination of temporal timeserise
  //#for example if temporal =1, the index is 1, then the data will be column 0 and column 1
  //# if temporal =2, the index is 2, then data will be 1) column 0 abd column 2,;2) column 1 and colun 2
  Rcpp::List getAllTemporalStates = generateTemporalGeneStates(main_parameters_in_ref,
                                                               target_gene,
                                                               conditional_genes,
                                                               temporal);
  
  List resultGroup(getAllTemporalStates.length());
  bool recount_target = false;
  Rcpp::List new_targetCounts(getAllTemporalStates.length());
  if(targetCounts.isNull()){
    recount_target= true;
  }else{
    new_targetCounts = targetCounts.as();
  }
  
  for(int i=0; i < getAllTemporalStates.length(); i++){
    List temporalState = getAllTemporalStates[i];
    int time_step = temporalState["timeStep"];
    int total_calculated_timepoints = 0;
    total_calculated_timepoints = n_timepoints-(total_samples * time_step);
    List result = getBasicMeasures(
      stateTCond,
      temporalState["computation_Matrix"],
                   temporalState["computation_Matrix_c"],
                                cond_T_target_T_state,
                                cond_F_target_T_state,
                                cond_T_target_F_state,
                                cond_F_target_F_state,
                                cond_T_target_T_state_c,
                                cond_F_target_T_state_c,
                                cond_T_target_F_state_c,
                                cond_F_target_F_state_c,
                                recount_target
    );
    
    if(recount_target){
      IntegerVector targets = IntegerVector::create(_["target_T_count"]=0, _["target_F_count"]=0);
      targets["target_T_count"] = result["target_T_count"];
      targets["target_F_count"] = result["target_F_count"];
      new_targetCounts[i] = targets;
    }else{
      IntegerVector targets = new_targetCounts[i];
      result["target_T_count"] = targets["target_T_count"];
      result["target_F_count"] = targets["target_F_count"];
    }
    
    int len = result.length();
    CharacterVector res_names = result.names();
    
    result = resizel(result,len+3);
    result[len] = total_calculated_timepoints;
    res_names.push_back("total_calculated_timepoints");
    
    result[len + 1] = num_of_conditional_genes;
    res_names.push_back("num_of_conditional_genes");
    
    result[len + 2] = time_step;
    res_names.push_back("timestep");
    
    result.attr("names") = Rcpp::wrap(res_names);
    resultGroup[i] = result;
  }
  resultGroup.attr("class") = "FBNBasicMeasure";
  return( resultGroup);
}

// [[Rcpp::export]]
Rcpp::List getAdvancedMeasures(const Rcpp::List basic_measures){
    //extract information
    double cond_T_count  =  basic_measures["cond_T_count"];
    double cond_F_count  =  basic_measures["cond_F_count"];
    double target_T_count = basic_measures["target_T_count"];
    double target_F_count = basic_measures["target_F_count"];
    double cond_T_count_c  = basic_measures["cond_T_count_c"];
    double cond_F_count_c  = basic_measures["cond_F_count_c"];
  
    int time_step = basic_measures["timestep"];
    //#condition#####################################################################################333333
    double lenTT = basic_measures["lenTT"];
    double lenTF = basic_measures["lenTF"];
    double lenFT = basic_measures["lenFT"];
    double lenFF = basic_measures["lenFF"];
     
    double total_calculated_timepoints = basic_measures["total_calculated_timepoints"];
    //#condition counter################################################################################
    double lenTT_c = basic_measures["lenTT_c"];
    double lenTF_c = basic_measures["lenTF_c"];
    double lenFT_c = basic_measures["lenFT_c"];
    double lenFF_c = basic_measures["lenFF_c"];
   
    double condition_T_support = dround(cond_T_count/total_calculated_timepoints,5);
    double condition_F_support = dround(cond_F_count/total_calculated_timepoints,5);
    double target_T_support = dround(target_T_count/total_calculated_timepoints,5);
    double target_F_support = dround(target_F_count/total_calculated_timepoints,5);
  
    //#final p(B if A)=p(B and A)/p(A), A=conditions, B is the target
    double confidence_TT = 0;
    double confidence_FT = 0;
    if(cond_T_count > 0){
      confidence_TT = dround(lenTT/cond_T_count,5);
      confidence_FT = dround(lenFT/cond_T_count,5);
    }
  
    double confidence_TF = 0;
    double confidence_FF = 0;
    if(cond_F_count > 0){
      confidence_TF = dround(lenTF/cond_F_count,5);
      confidence_FF = dround(lenFF/cond_F_count,5);
    }
  
  
    //#final p(A if B)=p(B and A)/p(B), A=conditions, B is the target
    double counter_confidence_TT = 0;
    double counter_confidence_FT = 0;
    double counter_confidence_TF = 0;
    double counter_confidence_FF = 0;
    if(cond_T_count_c > 0){
      counter_confidence_TT = dround(lenTT_c/cond_T_count_c,5);
      counter_confidence_FT = dround(lenFT_c/cond_T_count_c,5);
    }
  
    if(cond_F_count_c > 0){
      counter_confidence_TF = dround(lenTF_c/cond_F_count_c,5);
      counter_confidence_FF = dround(lenFF_c/cond_F_count_c,5);
    }
  
     //###############################end counter####################################################################
  
    Rcpp::NumericVector vectors = NumericVector::create(lenTT + 1,lenFT + 1,lenTF + 1,lenFF + 1);
    vectors.attr("dim") = Dimension(2, 2);
    Rcpp::NumericMatrix pTable = as<NumericMatrix>(vectors);
    //df = (r-1)(c-1) where r is the number of rows and c is the number of columns.
    //chiSQ = chisq.test(pTable,correct = FALSE,simulate.p.value = TRUE)
    //NumericVector chiSQ = Rcpp::qchisq(v,1);
    Rcpp::List pTest =fisher_test_cpp(pTable);
    //#from the book Data mining concepts and techniques
    //#isNegativeCorrelated means the inhibited of the conditional genes regulate the target gene
    //#isPossitiveCorrelated means the activated of the conditional genes regulate the target gene
  
    bool isNegativeCorrelated=false;
    bool isPossitiveCorrelated = false;
    double test1 = (lenTF/total_calculated_timepoints)*(lenFT/total_calculated_timepoints);
    double test2 = (lenTT/total_calculated_timepoints)*(lenFF/total_calculated_timepoints);
    if(total_calculated_timepoints > 0 && test1 > test2){
      isNegativeCorrelated = true;
    }
  
    if(total_calculated_timepoints > 0  && test1 < test2){
      isPossitiveCorrelated = true;
    }
  
  
  //#at moment use total_calculated_timepoints receive the best result with 0.00013
    double supportTT = 0;
    if(total_calculated_timepoints > 0){
      supportTT = dround(lenTT/total_calculated_timepoints,5);
    }
    double supportFT = 0;
    if(total_calculated_timepoints > 0){
      supportFT = dround(lenFT/total_calculated_timepoints,5);
    }
    double supportTF = 0;
    if(total_calculated_timepoints > 0){
      supportTF = dround(lenTF/total_calculated_timepoints,5);
    }
    double supportFF = 0;
    if(total_calculated_timepoints > 0){
      supportFF = dround(lenFF/total_calculated_timepoints,5);
    }
  
  //#calculate all confidence and max confidence
    double max_confidence_TT = dround(fmax(confidence_TT,counter_confidence_TT),5);
    double all_confidence_TT = dround(fmin(confidence_TT,counter_confidence_TT),5);
    double max_confidence_TF = dround(fmax(confidence_TF,counter_confidence_TF),5);
    double all_confidence_TF = dround(fmin(confidence_TF,counter_confidence_TF),5);
    double max_confidence_FT = dround(fmax(confidence_FT,counter_confidence_FT),5);
    double all_confidence_FT = dround(fmin(confidence_FT,counter_confidence_FT),5);
    double max_confidence_FF = dround(fmax(confidence_FF,counter_confidence_FF),5);
    double all_confidence_FF = dround(fmin(confidence_FF,counter_confidence_FF),5);
  
  //#conditional causality test: very important and effective measures, >=1 means the cause relationship is from the condition to target
  
    double causality_test_TT = 99999;
    if(counter_confidence_TT != 0){
      causality_test_TT = dround(confidence_TT/counter_confidence_TT,2);
    }
  
    double causality_test_TF =  99999;
    if(counter_confidence_TF!=0){
      causality_test_TF = dround(confidence_TF/counter_confidence_TF,2);
    }
  
    double causality_test_FT = 99999;
    if(counter_confidence_FT!=0){
      causality_test_FT = dround(confidence_FT/counter_confidence_FT,2);
    }
  
    double causality_test_FF =  99999;
    if(counter_confidence_FF!=0){
      causality_test_FF = dround(confidence_FF/counter_confidence_FF,2);
    }
  
    double signal_activator = 0;
    double signal_inhibitor = 0;
    double error_activator = 0;
    double error_inhibitor = 0;
    double pickT_support = 0;
    double pickT_causality_test = 0;
    double pickT_confidenceCounter = 0;
    double pickT_all_confidence = 0;
    double pickT_max_confidence = 0 ;
    std::string signal_sign_T ="";
    double pickF_support = 0;
    double pickF_causality_test = 0;
    double pickF_confidenceCounter = 0;
    double pickF_all_confidence = 0;
    double pickF_max_confidence = 0 ;
    std::string signal_sign_F ="";

    if (confidence_TT >= confidence_TF)
    {
      signal_activator = confidence_TT;
      error_activator = 1 - confidence_TT;
      pickT_support = supportTT;
      pickT_causality_test = causality_test_TT;
      pickT_confidenceCounter = counter_confidence_TT;
      pickT_all_confidence = all_confidence_TT;
      pickT_max_confidence = max_confidence_TT;
      signal_sign_T = "TT";

    }else
    {
      signal_activator = confidence_TF;
      error_activator = 1 - confidence_TF;
      pickT_support = supportTF;
      pickT_causality_test = causality_test_TF;
      pickT_confidenceCounter = counter_confidence_TF;
      pickT_all_confidence = all_confidence_TF;
      pickT_max_confidence = max_confidence_TF;
      signal_sign_T = "TF";
    }

    if(confidence_FT >= confidence_FF)
    {
      signal_inhibitor = confidence_FT;
      error_inhibitor = 1 - confidence_FT;
      pickF_support = supportFT;
      pickF_causality_test = causality_test_FT;
      pickF_confidenceCounter = counter_confidence_FT;
      pickF_all_confidence = all_confidence_FT;
      pickF_max_confidence = max_confidence_FT;
      signal_sign_F = "FT";
    }else
    {
      signal_inhibitor = confidence_FF;
      error_inhibitor = 1 - confidence_FF;
      pickF_support = supportFF;
      pickF_causality_test = causality_test_FF;
      pickF_confidenceCounter = counter_confidence_FF;
      pickF_all_confidence = all_confidence_FF;
      pickF_max_confidence = max_confidence_FF;
      signal_sign_F = "FF";
    }

        bool is_Essential = true;
        if (signal_activator == 1 && confidence_TT == confidence_TF)
          is_Essential = false;
  
        if (signal_inhibitor == 1 && confidence_FT == confidence_FF)
          is_Essential = false;
  
        if (!isNegativeCorrelated && !isPossitiveCorrelated)
          is_Essential = false;
  
        double p_value = (double)pTest["p.value"];
        if (p_value > 0.05)
          is_Essential = false;

        int essential = 0;
        if(is_Essential)essential=1;

        int causality_test_T=0;
        if(pickT_causality_test>=1)causality_test_T=1;

        int causality_test_F = 0;
        if(pickF_causality_test >= 1)causality_test_F=1;

        //key values for temporal FBN
        double bestFitP = sqrt(pow((pickT_max_confidence - signal_activator),2) + pow((pickT_all_confidence - pickT_confidenceCounter),2) + pow((signal_activator - 1),2) + pow((causality_test_T - 1),2) + pow((essential - 1),2));
        double bestFitN = sqrt(pow((pickF_max_confidence - signal_inhibitor),2) + pow((pickF_all_confidence - pickF_confidenceCounter),2) + pow((signal_inhibitor - 1),2) + pow((causality_test_F - 1),2) + pow((essential - 1),2));


        if(std::isinf(bestFitP) || isReallyNA(bestFitP))
          bestFitP = 99999;

        if(std::isinf(bestFitN) || isReallyNA(bestFitN))
          bestFitN = 99999;
        
        std::vector<std::string> names;
        
        List result(42);//42
        result[0]=confidence_TT;
        names.push_back("TT");
        result[1]=confidence_TF;
        names.push_back("TF");
        result[2]=confidence_FT;
        names.push_back("FT");
        result[3]=confidence_FF;
        names.push_back("FF");
        result[4]=counter_confidence_TT;
        names.push_back("TT_c");
        result[5]=counter_confidence_TF;
        names.push_back("TF_c");
        result[6]=counter_confidence_FT;
        names.push_back("FT_c");
        result[7]=counter_confidence_FF;
        names.push_back("FF_c");
        result[8]=condition_T_support;
        names.push_back("conditionT");
        result[9]=condition_F_support;
        names.push_back("conditionF");
        result[10]=target_T_support;
        names.push_back("targetT");
        result[11]=target_F_support;
        names.push_back("targetF");
        result[12]=isNegativeCorrelated;
        names.push_back("isNegativeCorrelated");
        result[13]=isPossitiveCorrelated;
        names.push_back("isPossitiveCorrelated");
        result[14]=supportTT;
        names.push_back("supportTT");
        result[15]=supportFT;
        names.push_back("supportFT");
        result[16]=supportTF;
        names.push_back("supportTF");
        result[17]=supportFF;
        names.push_back("supportFF");
        result[18]=signal_sign_T;
        names.push_back("signal_sign_T");
        result[19]=signal_sign_F;
        names.push_back("signal_sign_F");
        result[20]=time_step;
        names.push_back("timestep");
        result[21]=bestFitP;
        names.push_back("bestFitP");
        result[22]=bestFitN;
        names.push_back("bestFitN");
        result[23]=is_Essential;
        names.push_back("is_essential_gene");  
        result[24]=error_activator;
        names.push_back("Noise_P");   
        result[25]=error_inhibitor;
        names.push_back("Noise_N");  
        result[26]=signal_activator;
        names.push_back("Signal_P");    
        result[27]=signal_inhibitor;
        names.push_back("Signal_N");  
        result[28]=pickT_support;
        names.push_back("pickT_support");
        result[29]=pickF_support;
        names.push_back("pickF_support");
        result[30]=pickT_causality_test;
        names.push_back("pickT_causality_test");
        result[31]=pickF_causality_test;
        names.push_back("pickF_causality_test");
        result[32]=pickT_confidenceCounter;
        names.push_back("pickT_confidenceCounter");
        result[33]=pickF_confidenceCounter;
        names.push_back("pickF_confidenceCounter");
        result[34]=pickT_all_confidence;
        names.push_back("pickT_all_confidence");
        result[35]=pickF_all_confidence;
        names.push_back("pickF_all_confidence");
        result[36]=pickT_max_confidence;
        names.push_back("pickT_max_confidence");
        result[37]=pickF_max_confidence;
        names.push_back("pickF_max_confidence");
        result[38]=basic_measures;
        names.push_back("basic_measures");
        result[39]=p_value;
        names.push_back("p_value");
        
        result.attr("names") = Rcpp::wrap(names);
       return(result);
}




// [[Rcpp::export]]
Rcpp::List getGenePrababilities_advanced(const Rcpp::List getGenePrababilities_basic)
{
  
  //P(B|A)=P(A & B)/P(A)=(frq(A&B)/n)/(frq(A)/n)=frq(A&B)/frq(A)=P(A&D)/(P(A&D)+P(!A & D)) Bayers rule
  //order of target_genes is very important
  int len = getGenePrababilities_basic.length();
  List resultGroup(len);
  List targetCounts(len);
  //data cleaning
  List bestFitP;
  List bestFitN;
  for(int j=0;j<len;j++){
    List basic = getGenePrababilities_basic[j];
    resultGroup[j] = getAdvancedMeasures(basic);
    IntegerVector targets = IntegerVector::create(_["target_T_count"]=0, _["target_F_count"]=0);
    List temp = resultGroup[j];
    int timestep = temp["timestep"];
    targets["target_T_count"] = basic["target_T_count"];
    targets["target_F_count"] = basic["target_F_count"];
    if(j==0){
      bestFitP = temp;
      bestFitN = temp;
    }else{
      double this_bestFitP = (double)bestFitP["bestFitP"];
      double this_bestFitN = (double)bestFitN["bestFitN"];
      double this_temp_bestFitP = (double)temp["bestFitP"];
      double this_temp_bestFitN = (double)temp["bestFitN"];
      
      if(this_bestFitP > this_temp_bestFitP){
        bestFitP = temp;
      }
      
      if(this_bestFitN > this_temp_bestFitN){
        bestFitN = temp;
      }
      
      if(this_bestFitP == this_temp_bestFitP){
        if((int)bestFitP["timestep"]>timestep){
          bestFitP = temp;
        }
      }
      
      if(this_bestFitN == this_temp_bestFitN){
        if((int)bestFitN["timestep"]>timestep){
          bestFitN = temp;
        }
      }
    }
    
    targetCounts[j] = targets;
  }

  //Named("result_group") = resultGroup,
  return (List::create( _["getBestFitP"] = bestFitP, 
                        _["getBestFitN"] = bestFitN,
                        _["targetCounts"] =targetCounts));
}

//' The main function to main FBN probabilities from time series data
//' 
//' @param main_parameters_in_ref The environment that contains the required data.
//' @param fixedgenestate The up stream fixed gene state.
//' @param target_gene The target gene.
//' @param new_conditional_gene The current stream conditional gene.
//' @param temporal The temporal time step.
//' @param targetCounts A list of pre-calculated targe genes.
// [[Rcpp::export]]
Rcpp::List getGenePrababilities(Rcpp::Environment main_parameters_in_ref,
                                Rcpp::Nullable<Rcpp::List> fixedgenestate, 
                                Rcpp::CharacterVector target_gene, 
                                Rcpp::CharacterVector new_conditional_gene,
                                Rcpp::IntegerVector temporal, 
                                Rcpp::Nullable<Rcpp::List> targetCounts)
{
  
    Rcpp::List basic_measures = getGenePrababilities_basic(main_parameters_in_ref,
                                                           fixedgenestate,target_gene,
                                                           new_conditional_gene,temporal, 
                                                           targetCounts);
    Rcpp::List probability = getGenePrababilities_advanced(basic_measures);

    if(Rf_isNull(probability))
      return (R_NilValue);
    
    Rcpp::List new_targetCounts =  probability["targetCounts"];
    return (List::create(_["getBestFitP"] = probability["getBestFitP"], 
                         _["getBestFitN"] = probability["getBestFitN"], 
                         _["targetCounts"] = new_targetCounts));
}