#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <String.h>
#include<string>
#include<limits.h>
#include <math.h>
#include <stdexcept>
#include "fbn_utility.h"
#include "fbn_core.h"

using namespace Rcpp;

//' Get the main measurements based on the input data
//' 
//' @param targetGene The target gene
//' @param mainParameters An environment variable holds all input data
//' @param genes All conditional genes
//' @param matchedgenes processed genes
//' @param temporal The temporal time steps
//' @param targetCounts 
//' 
// [[Rcpp::export]]
Rcpp::List getGenePrababilities_measurements(
      Rcpp::CharacterVector targetGene,
      Rcpp::Environment mainParameters,
      Rcpp::CharacterVector genes,
      Rcpp::Nullable<Rcpp::List> matchedgenes,
      Rcpp::IntegerVector temporal=1,
      Nullable<Rcpp::List> targetCounts = R_NilValue)
{

   int len = genes.length();
   List res(len);

   for(int i = 0;i < len;i++){
      std::string gene = as<std::string>(genes[i]);
      List probabilityOfFourCombines = getGenePrababilities(mainParameters,
                                                            matchedgenes,
                                                            targetGene,
                                                            gene,
                                                            temporal,
                                                            targetCounts);
      if(Rf_isNull(probabilityOfFourCombines))
      {
         continue;
      }
      List probabilityOfFourCombines_P;
      List probabilityOfFourCombines_N;
      Rcpp::List new_targetCounts = probabilityOfFourCombines["targetCounts"];
      probabilityOfFourCombines_P = probabilityOfFourCombines["getBestFitP"];
      probabilityOfFourCombines_N = probabilityOfFourCombines["getBestFitN"];

      if(!probabilityOfFourCombines_P["is_essential_gene"] 
            && !probabilityOfFourCombines_N["is_essential_gene"]){
         continue;
      }
      
      if(!probabilityOfFourCombines_P["is_essential_gene"])
      {
         probabilityOfFourCombines_P=probabilityOfFourCombines_N;
      }
      
      if(!probabilityOfFourCombines_N["is_essential_gene"]){
         probabilityOfFourCombines_N=probabilityOfFourCombines_P;
      }
      
      
      if(!probabilityOfFourCombines_P["isPossitiveCorrelated"] 
            && !probabilityOfFourCombines_P["isNegativeCorrelated"]){
         probabilityOfFourCombines_P=probabilityOfFourCombines_N;
      }
      
      if(!probabilityOfFourCombines_N["isPossitiveCorrelated"] 
            && !probabilityOfFourCombines_N["isNegativeCorrelated"])
      {
         if(is_true(Rcpp::all(probabilityOfFourCombines_P==probabilityOfFourCombines_N)))
            continue;
         
         probabilityOfFourCombines_N=probabilityOfFourCombines_P;
      }
      
      List in_res(3);
      in_res[0]= probabilityOfFourCombines_P;
      in_res[1]= probabilityOfFourCombines_N;
      in_res[2]=new_targetCounts;
      CharacterVector in_res_names = CharacterVector::create("probabilityOfFourCombines_P",
                                                             "probabilityOfFourCombines_N",
                                                             "targetCounts");

      in_res.attr("names") = in_res_names;
      res[i] = in_res;
   }
   res.attr("names") = genes;
   res = removeEmptyElement(res);
   return(res);
}



// [[Rcpp::export]]
Rcpp::List buildProbabilityTreeOnTargetGene(
    Rcpp::CharacterVector targetGene,
    Rcpp::Environment mainParameters,
    Rcpp::CharacterVector genes,
    Rcpp::Nullable<Rcpp::List> matchedgenes,
    Rcpp::Nullable<Rcpp::CharacterVector> matchedexpression,
    Rcpp::IntegerVector maxK=4,
    Rcpp::IntegerVector temporal=1,
    Nullable<Rcpp::List> targetCounts = R_NilValue,
    bool findPositiveRegulate = false,
    bool findNegativeRegulate = false)
{
 //tree mining stop when 1 or 0
   Rcpp::List measuements = getGenePrababilities_measurements(targetGene, mainParameters, genes, matchedgenes, temporal, targetCounts);
   CharacterVector new_genes = measuements.names();
   CharacterVector unprocessedGenes = clone(new_genes);
   int len = new_genes.length();

   CharacterVector processedGenes;
   List res(len);
   maxK = ifelse(maxK>unprocessedGenes.length(),unprocessedGenes.length(),maxK);

  for(int i=0;i<len;i++){
     int pmaxK = (int)maxK[0];

    std::string gene = (std::string)new_genes[i];
    IntegerVector unprocessed_index = a_not_in_b_index(unprocessedGenes,CharacterVector::create(new_genes[i]));
    unprocessedGenes = unprocessedGenes[unprocessed_index];
    Rcpp::List newMatchedGenesT;
    Rcpp::List newMatchedGenesF;
    CharacterVector preprocessed;
    if(!matchedgenes.isNull())
    {
      List temp_matchedgenes(matchedgenes);
      newMatchedGenesT=temp_matchedgenes;
      newMatchedGenesF=temp_matchedgenes;
      preprocessed= temp_matchedgenes.names();
    }

    CharacterVector expressionT;
    CharacterVector expressionF;

  if(!matchedexpression.isNull())
  {
     Rcpp::CharacterVector pmatchedexpression(matchedexpression);
     expressionT = CharacterVector::create("&", gene);
     expressionF = CharacterVector::create("&","!",gene);
     expressionT = concatenator(pmatchedexpression,expressionT);
     expressionF = concatenator(pmatchedexpression,expressionF);

     if(!newMatchedGenesT.containsElementNamed(gene.c_str()))
     {
       int index=newMatchedGenesT.length()+1;
       CharacterVector temp_name = newMatchedGenesT.names();
       temp_name.push_back(gene);
       newMatchedGenesT = resizel(newMatchedGenesT,index);
       newMatchedGenesT[index-1]=1;
       newMatchedGenesT.attr("names") = Rcpp::wrap(temp_name);
     }


     if(!newMatchedGenesF.containsElementNamed(gene.c_str()))
     {
       int index=newMatchedGenesF.length()+1;
       CharacterVector temp_name = newMatchedGenesF.names();
       temp_name.push_back(gene);
       newMatchedGenesF = resizel(newMatchedGenesF,index);
       newMatchedGenesF[index-1]=0;
       newMatchedGenesF.attr("names") = Rcpp::wrap(temp_name);
     }

   }else
   {

     newMatchedGenesT=  List::create(Named(gene) = 1);
     newMatchedGenesF=  List::create(Named(gene) = 0);
     expressionT=CharacterVector::create(gene);
     expressionF=CharacterVector::create("!",gene);
   }

   String expT=mpaste(expressionT);
   String expF=mpaste(expressionF);

   Rcpp::CharacterVector inputgenes=newMatchedGenesT.names();
   inputgenes = char_sort(inputgenes,true);


   if(preprocessed.length()>0){
     if(is_true(any(a_in_b(gene,preprocessed)))){
       continue;
     }
   }

   List probabilityOfFourCombines_P;
   List probabilityOfFourCombines_N;
   List measuement = measuements[i];
   probabilityOfFourCombines_P = measuement["probabilityOfFourCombines_P"];
   probabilityOfFourCombines_N = measuement["probabilityOfFourCombines_N"];
   Rcpp::List new_targetCounts = measuement["targetCounts"];

   double bestFit_P=probabilityOfFourCombines_P["bestFitP"];
   bool is_essential_gene_P=probabilityOfFourCombines_P["is_essential_gene"];
   double bestFit_N=probabilityOfFourCombines_N["bestFitN"];
   bool is_essential_gene_N=probabilityOfFourCombines_N["is_essential_gene"];
   String sign_P=probabilityOfFourCombines_P["signal_sign_T"];
   String sign_N=probabilityOfFourCombines_N["signal_sign_F"];

   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ///find next branch
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   List subresultT;
   List subresultF;

   CharacterVector exlcudedSubgenes;
   CharacterVector nextGenes_T;
   CharacterVector nextGenes_F;
   if(pmaxK>1 && !(findPositiveRegulate && findNegativeRegulate))
   {
      findPositiveRegulate = findPositiveRegulate || bestFit_P == 0;
      findNegativeRegulate = findNegativeRegulate || bestFit_N == 0;
      pmaxK=pmaxK-1;
      exlcudedSubgenes=newMatchedGenesT.names();
      IntegerVector unprocessed_index2 = a_not_in_b_index(unprocessedGenes, exlcudedSubgenes);
      nextGenes_T=unprocessedGenes[unprocessed_index2];

      if(nextGenes_T.length()>0 && is_essential_gene_P && (bestFit_P>0 || bestFit_N>0))
      {
         subresultT=buildProbabilityTreeOnTargetGene(targetGene,mainParameters,nextGenes_T,newMatchedGenesT,CharacterVector::create(expT),IntegerVector::create(pmaxK),temporal, new_targetCounts, findPositiveRegulate, findNegativeRegulate);
      }

      exlcudedSubgenes=newMatchedGenesF.names();
      unprocessed_index2 = a_not_in_b_index(unprocessedGenes, exlcudedSubgenes);
      nextGenes_F=unprocessedGenes[unprocessed_index2];
      if(nextGenes_F.length()>0 && is_essential_gene_N && (bestFit_P>0 || bestFit_N>0))
      {
         subresultF=buildProbabilityTreeOnTargetGene(targetGene,mainParameters,nextGenes_F,newMatchedGenesF,CharacterVector::create(expF),IntegerVector::create(pmaxK),temporal, new_targetCounts, findPositiveRegulate, findNegativeRegulate);
      }
   }

   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //positive regulation
   double pickT_support=probabilityOfFourCombines_P["pickT_support"];
   double pickT_causality_test=probabilityOfFourCombines_P["pickT_causality_test"];
   double pickT_value=probabilityOfFourCombines_P["Signal_P"];
   double pickT_noise=probabilityOfFourCombines_P["Noise_P"];
   double pickT_confidenceCounter=probabilityOfFourCombines_P["pickT_confidenceCounter"];

   double pickT_all_confidence=probabilityOfFourCombines_P["pickT_all_confidence"];
   double pickT_max_confidence=probabilityOfFourCombines_P["pickT_max_confidence"];
   bool isNegativeCorrelated_T=probabilityOfFourCombines_P["isNegativeCorrelated"];
   bool isPossitiveCorrelated_T=probabilityOfFourCombines_P["isPossitiveCorrelated"];
   int timestep_T=probabilityOfFourCombines_P["timestep"];
   double bestFitP = probabilityOfFourCombines_P["bestFitP"];
   //#negative regulation
   double pickF_support=probabilityOfFourCombines_N["pickF_support"];
   double pickF_causality_test=probabilityOfFourCombines_N["pickF_causality_test"];
   double pickF_value=probabilityOfFourCombines_N["Signal_N"];
   double pickF_noise=probabilityOfFourCombines_N["Noise_N"];
   double pickF_ConfidenceCounter=probabilityOfFourCombines_N["pickF_confidenceCounter"];

   double pickF_all_confidence=probabilityOfFourCombines_N["pickF_all_confidence"];
   double pickF_max_confidence=probabilityOfFourCombines_N["pickF_max_confidence"];
   bool isNegativeCorrelated_F=probabilityOfFourCombines_N["isNegativeCorrelated"];
   bool isPossitiveCorrelated_F=probabilityOfFourCombines_N["isPossitiveCorrelated"];
   int timestep_F=probabilityOfFourCombines_N["timestep"];
   double bestFitN = probabilityOfFourCombines_N["bestFitN"];
      
   String pick_expT;
   String pick_expF;
   String identity_T;
   String identity_F;
   CharacterVector pattern;
   if(sign_P=="TT")
   {
     for(int i=0;i<inputgenes.length();i++){
       std::string s_gene = (std::string)inputgenes[i];
       int gene_state = newMatchedGenesT[s_gene];
       std::string this_gene = (std::string)inputgenes[i] + "$" + to_string(gene_state);
       pattern.push_back(this_gene);
     }
     identity_T=mpaste(pattern, "_");
     pick_expT = expT;

   }else if(sign_P=="TF")
   {
     for(int i=0;i<inputgenes.length();i++){
       std::string s_gene = (std::string)inputgenes[i];
       int gene_state = newMatchedGenesF[s_gene];
       std::string this_gene = (std::string)inputgenes[i] + "$" + to_string(gene_state);
       pattern.push_back(this_gene);
     }
     identity_T=mpaste(pattern,"_");
     pick_expT = expF;
   }
   pattern = CharacterVector::create((std::string)identity_T);
   pattern.push_back("Activator_of");
   pattern.push_back(targetGene[0]);
   identity_T=mpaste(pattern,"_");

   CharacterVector activator = CharacterVector::create(_["factor"]=pick_expT, 
                                                       _["Confidence"]=(String)pickT_value, 
                                                       _["ConfidenceCounter"]=(String)pickT_confidenceCounter,
                                                       _["all_confidence"]=(String)pickT_all_confidence, 
                                                       _["max_confidence"]=(String)pickT_max_confidence, 
                                                       _["support"]=(String)pickT_support,
                                                       _["causality_test"]=(String)pickT_causality_test,
                                                       _["Noise"]=(String)pickT_noise, _["Identity"]=identity_T,
                                                       _["type"]=sign_P, _["timestep"]=(String)timestep_T,
                                                       _["isNegativeCorrelated"]=(String)isNegativeCorrelated_T,
                                                       _["isPossitiveCorrelated"]=(String)isPossitiveCorrelated_T, 
                                                       _["bestFitP"]=(String)bestFitP);


   pattern = CharacterVector::create();
   if(sign_N=="FT")
   {
     for(int i=0;i<inputgenes.length();i++){
       std::string s_gene = (std::string)inputgenes[i];
       int gene_state = newMatchedGenesT[s_gene];
       std::string this_gene = (std::string)inputgenes[i] + "$" + to_string(gene_state);
       pattern.push_back(this_gene);
     }
     identity_F=mpaste(pattern,"_");
     pick_expF = expT;

   }else if (sign_N=="FF")
   {
     for(int i=0;i<inputgenes.length();i++){
       std::string s_gene = (std::string)inputgenes[i];
       int gene_state = newMatchedGenesF[s_gene];
       std::string this_gene = (std::string)inputgenes[i] + "$" + to_string(gene_state);
       pattern.push_back(this_gene);
     }
     identity_F=mpaste(pattern,"_");
     pick_expF = expF;
   }
   pattern = CharacterVector::create((std::string)identity_F);
   pattern.push_back("Inhibitor_of");
   pattern.push_back(targetGene[0]);
   identity_F=mpaste(pattern,"_");

   CharacterVector inhibitor = CharacterVector::create(_["factor"]=pick_expF, 
                                                       _["Confidence"]=(String)pickF_value, 
                                                       _["ConfidenceCounter"]=(String)pickF_ConfidenceCounter,
                                                       _["all_confidence"]=(String)pickF_all_confidence, 
                                                       _["max_confidence"]=(String)pickF_max_confidence, 
                                                       _["support"]=(String)pickF_support,
                                                       _["causality_test"]=(String)pickF_causality_test,
                                                       _["Noise"]=(String)pickF_noise,_["Identity"]=identity_F,
                                                       _["type"]=sign_N, _["timestep"]=(String)timestep_F,
                                                       _["isNegativeCorrelated"]=(String)isNegativeCorrelated_F, 
                                                       _["isPossitiveCorrelated"]=(String)isPossitiveCorrelated_F,
                                                       _["bestFitN"]=(String)bestFitN);


   List a_i_tor(2);
   a_i_tor[0] = activator;
   a_i_tor[1] = inhibitor;
   a_i_tor.attr("names") = CharacterVector::create("Activator", "Inhibitor");

   List in_res(2);
   in_res[0]= a_i_tor;
   in_res[1]= inputgenes;
   CharacterVector in_res_names = CharacterVector::create("ActivatorAndInhibitor", "Input");

   if(subresultT.length()>0){
      int oindex = in_res.length();
      in_res = resizel(in_res,oindex+1);
      in_res[oindex]=subresultT;
      in_res_names.push_back("SubGenesT");
   }

   if(subresultF.length()>0){
      int oindex = in_res.length();
      in_res = resizel(in_res,oindex+1);
      in_res[oindex]=subresultF;
      in_res_names.push_back("SubGenesF");
   }
   in_res.attr("names") = in_res_names;
   res[i] = in_res;
 }
 res.attr("names") = new_genes;
 res = removeEmptyElement(res);
  return(res);
}



//' Get the main measurements based on the input data
//' 
//' @param target_gene The target gene
//' @param conditional_genes conditional genes
//' @param maxK The maximum under ground levels
//' @param temporal The temporal time steps
//' @param mainParameters An environment variable holds all input data
//' 
// [[Rcpp::export]]
Rcpp::List process_cube_algorithm(Rcpp::CharacterVector target_gene, 
                                Rcpp::CharacterVector conditional_genes, 
                                Rcpp::IntegerVector maxK, 
                                Rcpp::IntegerVector temporal, 
                                Rcpp::Environment mainParameters) {
   List res(1);
   List sub_res(1);
   sub_res[0] = buildProbabilityTreeOnTargetGene(target_gene, mainParameters, conditional_genes, R_NilValue, R_NilValue, maxK, temporal, R_NilValue);
   sub_res.attr("names") = CharacterVector::create("SubGenes");
   res[0] = sub_res;
   res.attr("names") = target_gene;
   return(res);
}