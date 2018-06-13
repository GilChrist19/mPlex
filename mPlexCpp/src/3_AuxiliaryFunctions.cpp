/*                  ____  __
 *       ____ ___  / __ \/ /__  _  __
 *      / __ `__ \/ /_/ / / _ \| |/_/
 *     / / / / / / ____/ /  __/>  <
 *    /_/ /_/ /_/_/   /_/\___/_/|_|
 *
 *    Marshall Lab
 *    Mosquito class
 *    January 2018
 */

#include "3_DriveDefinitions.hpp"

/******************************************************************************
 * Mosquito Creation Functions
 ******************************************************************************/
// This one is for Daisy and multiLocus
void CreateMosquitoes2Allele(int numMos, int minAge, dVec ageDist,
                             Rcpp::ListOf<Rcpp::List> aTypes, popVec returnPop){
  
  // This only good for things that have 2 alleles per locus
  
  // number of loci is the number of sublists in aTypes
  size_t numLoci = aTypes.size();
  
  // holders for things in loop
  std::vector<double> probs;
  std::vector<std::string> alleles;
  std::vector<std::string> genHold(2);
  std::string genotype;
  int age;
  
  
  // loop over number of mosquitoes to create
  for(size_t count=0; count < numMos; ++count){
    
    
    // loop over each locus
    for(size_t loci=0; loci < numLoci; ++ loci){
      
      probs = Rcpp::as<std::vector<double> >( aTypes[loci]["probs"] );
      alleles = Rcpp::as<std::vector<std::string> >( aTypes[loci]["alleles"] );
      
      // generate each allele
      for(size_t gen=0; gen < 2; gen++){
        genHold[gen] = alleles[prng::instance().get_oneSample(probs)];
      }
      
      // sort alleles
      std::sort(genHold.begin(), genHold.end());
      
      // combine alleles into final genotype
      for(auto const& s : genHold){ genotype += s; }
    }
    
    // set age of mosquito
    age = minAge + prng::instance().get_oneSample(ageDist);
    
    // add new mosquito to population
    returnPop.push_back(Mosquito(age, genotype));
    
    // clear genotype for next one
    genotype.clear();
  } // end loop over mosquitoes
  
}

// This one is for oneLocus
void CreateMosquitoes2Loci(int numMos, int minAge, dVec ageDist,
                           Rcpp::ListOf<Rcpp::List> aTypes, popVec returnPop){
  
  // This only good for things that have 2 loci, but multiple alleles
  
  // number of loci is the number of sublists in aTypes
  size_t numLoci = aTypes.size();
  
  // holders for things in loop
  std::vector<double> probs;
  std::vector<std::string> alleles;
  
  std::vector<std::string> locus1(numLoci), locus2(numLoci);
  std::string genotype;
  int age;
  
  
  // loop over number of mosquitoes to create
  for(size_t count=0; count < numMos; ++count){
    
    // loop over each locus
    for(size_t loci=0; loci < numLoci; ++ loci){
      
      probs = Rcpp::as<std::vector<double> >( aTypes[loci]["probs"] );
      alleles = Rcpp::as<std::vector<std::string> >( aTypes[loci]["alleles"] );
      
      // generate each allele
      locus1[loci] = alleles[prng::instance().get_oneSample(probs)];
      locus2[loci] = alleles[prng::instance().get_oneSample(probs)];
    }
    
    // combine alleles into final genotype
    for(auto const& s : locus1){ genotype += s; }
    for(auto const& s : locus2){ genotype += s; }
    
    
    // set age of mosquito
    age = minAge + prng::instance().get_oneSample(ageDist);
    
    // add new mosquito to population
    returnPop.push_back(Mosquito(age, genotype));
    
    // clear genotype for next one
    genotype.clear();
  } // end loop over mosquitoes
  
}

/******************************************************************************
 * Mosquito Generating Functions
 ******************************************************************************/

/**************************************
 * Daisy Generating Function
 **************************************/
void DaisyOffspring(const std::string& fGen, const std::string& mGen){
  
  
  //// list of objects I am currently creating ////
  
  // int numAlleles
  
  // vector<bool> fScore, mScore
  // int index
  
  // dMat fProbs, mProbs
  // sMat fAllele, mAllele
  
  // std::string holdAllele;
  // std::unordered_map<std::string, double>   duplicates;
  // std::unordered_map<std::string, double>::iterator value;
  
  
  // sVec finalGenotypes;
  // dVec finalProbs;
  
  // sVec holdGens;
  // dVec holdProbs;
  
  
  
  
  
  
  
  
  
  
  
  // get number of alleles
  int numAlleles = fGen.size()/2;
  
  /*****************************************************************************/
  // Score Each Allele
  /*****************************************************************************/
  std::vector<bool> fScore(numAlleles+1, false), mScore(numAlleles+1, false);
  
  //loop over loci, separate alleles and score
  int index=1;
  for(size_t i=0; i<fGen.size(); i+=2, ++index){
    // female score
    if( (fGen[i] == 'H') || (fGen[i+1] == 'H')) {fScore[index] = true;}
    // male score
    if( (mGen[i] == 'H') || (mGen[i+1] == 'H')) {mScore[index] = true;}
    
  } // end scoring loop
  /*****************************************************************************/
  //End Split and Score
  /*****************************************************************************/
  
  /*****************************************************************************/
  // Determine Next-Gen alleles
  /*****************************************************************************/
  dMat fProbs(numAlleles), mProbs(numAlleles);
  sMat fAllele(numAlleles), mAllele(numAlleles);
  
  // loop over all loci
  index=0;
  for(size_t i=0; i<numAlleles; ++i, index+=2){
    
    // FEMALES
    if( (!fScore[i] && !fScore[i+1]) || (!fScore[i] && fScore[i+1]) ){
      //FF  or FT case
      // loop over alleles, all loci diploid
      for(size_t j=0; j<2; ++j){
        // fill allele and probs
        if(fGen[index+j] == 'W'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_mendelian_allele(i,0).begin(),
                            reference::instance().get_mendelian_allele(i,0).end());
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_mendelian_probs(i,0).begin(),
                           reference::instance().get_mendelian_probs(i,0).end());
        } else if(fGen[index+j] == 'H'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_mendelian_allele(i,1).begin(),
                            reference::instance().get_mendelian_allele(i,1).end());
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_mendelian_probs(i,1).begin(),
                           reference::instance().get_mendelian_probs(i,1).end());
        } else if(fGen[index+j] == 'R'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_mendelian_allele(i,2).begin(),
                            reference::instance().get_mendelian_allele(i,2).end());
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_mendelian_probs(i,2).begin(),
                           reference::instance().get_mendelian_probs(i,2).end());
        } else if(fGen[index+j] == 'S'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_mendelian_allele(i,3).begin(),
                            reference::instance().get_mendelian_allele(i,3).end());
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_mendelian_probs(i,3).begin(),
                           reference::instance().get_mendelian_probs(i,3).end());
        }
      } // end allele loop
      
    } else if(fScore[i] && !fScore[i+1]){
      //TF case
      // loop over alleles, all loci diploid
      for(size_t j=0; j<2; ++j){
        // fill allele and probs
        if(fGen[index+j] == 'W'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_cutting_allele(i,0).begin(),
                            reference::instance().get_cutting_allele(i,0).end());
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_cutting_probs(i,0).begin(),
                           reference::instance().get_cutting_probs(i,0).end());
        } else if(fGen[index+j] == 'H'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_cutting_allele(i,1).begin(),
                            reference::instance().get_cutting_allele(i,1).end());
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_cutting_probs(i,1).begin(),
                           reference::instance().get_cutting_probs(i,1).end());
        } else if(fGen[index+j] == 'R'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_cutting_allele(i,2).begin(),
                            reference::instance().get_cutting_allele(i,2).end());
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_cutting_probs(i,2).begin(),
                           reference::instance().get_cutting_probs(i,2).end());
        } else if(fGen[index+j] == 'S'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_cutting_allele(i,3).begin(),
                            reference::instance().get_cutting_allele(i,3).end());
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_cutting_probs(i,3).begin(),
                           reference::instance().get_cutting_probs(i,3).end());
        }
      } // end allele loop
      
    } else if(fScore[i] && fScore[i+1]){
      //TT case
      // loop over alleles, all loci diploid
      for(size_t j=0; j<2; ++j){
        // fill allele and probs
        if(fGen[index+j] == 'W'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_homing_allele(i,0).begin(),
                            reference::instance().get_homing_allele(i,0).end());
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_homing_probs(i,0).begin(),
                           reference::instance().get_homing_probs(i,0).end());
        } else if(fGen[index+j] == 'H'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_homing_allele(i,1).begin(),
                            reference::instance().get_homing_allele(i,1).end());
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_homing_probs(i,1).begin(),
                           reference::instance().get_homing_probs(i,1).end());
        } else if(fGen[index+j] == 'R'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_homing_allele(i,2).begin(),
                            reference::instance().get_homing_allele(i,2).end());
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_homing_probs(i,2).begin(),
                           reference::instance().get_homing_probs(i,2).end());
        } else if(fGen[index+j] == 'S'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_homing_allele(i,3).begin(),
                            reference::instance().get_homing_allele(i,3).end());
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_homing_probs(i,3).begin(),
                           reference::instance().get_homing_probs(i,3).end());
        }
      } // end allele loop 
    } // end females
    
    // MALES
    if( (!mScore[i] && !mScore[i+1]) || (!mScore[i] && mScore[i+1]) ){
      //FF  or FT case
      // loop over alleles, all loci diploid
      for(size_t j=0; j<2; ++j){
        // fill allele and probs
        if(mGen[index+j] == 'W'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_mendelian_allele(i,0).begin(),
                            reference::instance().get_mendelian_allele(i,0).end());
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_mendelian_probs(i,0).begin(),
                           reference::instance().get_mendelian_probs(i,0).end());
        } else if(mGen[index+j] == 'H'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_mendelian_allele(i,1).begin(),
                            reference::instance().get_mendelian_allele(i,1).end());
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_mendelian_probs(i,1).begin(),
                           reference::instance().get_mendelian_probs(i,1).end());
        } else if(mGen[index+j] == 'R'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_mendelian_allele(i,2).begin(),
                            reference::instance().get_mendelian_allele(i,2).end());
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_mendelian_probs(i,2).begin(),
                           reference::instance().get_mendelian_probs(i,2).end());
        } else if(mGen[index+j] == 'S'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_mendelian_allele(i,3).begin(),
                            reference::instance().get_mendelian_allele(i,3).end());
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_mendelian_probs(i,3).begin(),
                           reference::instance().get_mendelian_probs(i,3).end());
        }
      } // end allele loop
      
    } else if(mScore[i] && !mScore[i+1]){
      //TF case
      // loop over alleles, all loci diploid
      for(size_t j=0; j<2; ++j){
        // fill allele and probs
        if(mGen[index+j] == 'W'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_cutting_allele(i,0).begin(),
                            reference::instance().get_cutting_allele(i,0).end());
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_cutting_probs(i,0).begin(),
                           reference::instance().get_cutting_probs(i,0).end());
        } else if(mGen[index+j] == 'H'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_cutting_allele(i,1).begin(),
                            reference::instance().get_cutting_allele(i,1).end());
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_cutting_probs(i,1).begin(),
                           reference::instance().get_cutting_probs(i,1).end());
        } else if(mGen[index+j] == 'R'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_cutting_allele(i,2).begin(),
                            reference::instance().get_cutting_allele(i,2).end());
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_cutting_probs(i,2).begin(),
                           reference::instance().get_cutting_probs(i,2).end());
        } else if(mGen[index+j] == 'S'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_cutting_allele(i,3).begin(),
                            reference::instance().get_cutting_allele(i,3).end());
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_cutting_probs(i,3).begin(),
                           reference::instance().get_cutting_probs(i,3).end());
        }
      } // end allele loop
      
      
    } else if(mScore[i] && mScore[i+1]){
      //TT case
      // loop over alleles, all loci diploid
      for(size_t j=0; j<2; ++j){
        // fill allele and probs
        if(mGen[index+j] == 'W'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_homing_allele(i,0).begin(),
                            reference::instance().get_homing_allele(i,0).end());
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_homing_probs(i,0).begin(),
                           reference::instance().get_homing_probs(i,0).end());
        } else if(mGen[index+j] == 'H'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_homing_allele(i,1).begin(),
                            reference::instance().get_homing_allele(i,1).end());
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_homing_probs(i,1).begin(),
                           reference::instance().get_homing_probs(i,1).end());
        } else if(mGen[index+j] == 'R'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_homing_allele(i,2).begin(),
                            reference::instance().get_homing_allele(i,2).end());
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_homing_probs(i,2).begin(),
                           reference::instance().get_homing_probs(i,2).end());
        } else if(mGen[index+j] == 'S'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_homing_allele(i,3).begin(),
                            reference::instance().get_homing_allele(i,3).end());
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_homing_probs(i,3).begin(),
                           reference::instance().get_homing_probs(i,3).end());
        }
      } // end allele loop
    } // end males
    
  } // end loci loop
  
  /*****************************************************************************/
  // End Next-Gen alleles
  /*****************************************************************************/
  
  /*****************************************************************************/
  // All Combinations of Male/Female for each Loci
  /*****************************************************************************/
  std::string holdAllele;
  
  // Aggregate duplicates before storing
  std::unordered_map<std::string, double>   duplicates;
  std::unordered_map<std::string, double>::iterator value;
  
  
  // loop over alleles
  for(size_t i=0; i<numAlleles; ++i){
    
    // all combinations of female/male alleles, and probs
    for(size_t fem=0; fem<fAllele[i].size(); ++fem){
      for(size_t mal=0; mal<mAllele[i].size(); ++mal){
        // combine and sort allele
        holdAllele = fAllele[i][fem] + mAllele[i][mal];
        std::sort(holdAllele.begin(), holdAllele.end() );
        
        // Here we aggregate non-unique values
        // check if it is in map
        value = duplicates.find(holdAllele);
        if(value == duplicates.end()){
          // not in map, add it
          duplicates.insert(std::make_pair(holdAllele, fProbs[i][fem]*mProbs[i][mal]));
        } else {
          // is in map, combine with existing value
          value->second += fProbs[i][fem]*mProbs[i][mal];
        }
        
      } // end male loop
    } // end female loop
    
    // clear fAllele and fProbs for re-use
    fAllele[i].clear();
    fProbs[i].clear();
    
    // Take unique values in map, return to matrix
    for(auto elem : duplicates){
      fAllele[i].push_back(elem.first);
      fProbs[i].push_back(elem.second);
    }
    
    // clear map
    duplicates.clear();
  } // end 
  
  /*****************************************************************************/
  // End All Combinations of MAle/Female for Each Loci
  /*****************************************************************************/
  
  /*****************************************************************************/
  // Cartesian Product of All Loci
  /*****************************************************************************/
  
  sVec finalGenotypes;
  dVec finalProbs;
  
  sVec holdGens;
  dVec holdProbs;
  
  // fill first index
  for(index=0; index<fAllele[0].size(); ++index){
    finalGenotypes.push_back(fAllele[0][index]);
    finalProbs.push_back(fProbs[0][index]);
  }
  
  // loop over rest of them
  for(index=1; index<numAlleles; ++index){
    // loop over current elements in return, and elements in next vector
    for(size_t current=0; current<finalGenotypes.size(); ++current){
      for(size_t newOne=0; newOne<fAllele[index].size(); ++newOne){
        // combine strings and probabilities
        holdGens.push_back(finalGenotypes[current] + fAllele[index][newOne]);
        holdProbs.push_back(finalProbs[current] * fProbs[index][newOne]);
      } // end loop over new things
    } // end loop over current things
    
    // set returns
    finalGenotypes = holdGens;
    finalProbs = holdProbs; 
    
    // clear holders
    holdGens.clear();
    holdProbs.clear();
  } // end loop over loci
  
  /*****************************************************************************/
  // End Cartesian Product of All Loci
  /*****************************************************************************/
  
} // end function











/**************************************
 * multiLocus Generating Function
 **************************************/
void MultiplexOffspring_mLoci(const std::string& fGen, const std::string& mGen){
  
  
  
  //// list of objects I am currently creating ////
  
  // int numAlleles
  
  // vector<bool> fScore, mScore
  // int index
  
  // dMat fProbs, mProbs
  // sMat fAllele, mAllele
  
  // std::string holdAllele;
  // std::unordered_map<std::string, double>   duplicates;
  // std::unordered_map<std::string, double>::iterator value;
  
  
  // sVec finalGenotypes;
  // dVec finalProbs;
  
  // sVec holdGens;
  // dVec holdProbs;
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  // get number of alleles
  int numAlleles = fGen.size()/2;
  
  /*****************************************************************************/
  // Score Each Allele
  /*****************************************************************************/
  std::vector<bool> fScore(numAlleles, false), mScore(numAlleles, false);
  
  //loop over loci, separate alleles and score
  int index=0;
  for(size_t i=0; i<fGen.size(); i+=2, ++index){
    // female score
    if( (fGen[i] == 'H') || (fGen[i+1] == 'H')) {fScore[index] = true;}
    // male score
    if( (mGen[i] == 'H') || (mGen[i+1] == 'H')) {mScore[index] = true;}
    
  } // end scoring loop
  /*****************************************************************************/
  //End Score Each Allele
  /*****************************************************************************/
  
  /*****************************************************************************/
  // Determine Next-Gen alleles
  /*****************************************************************************/
  sMat fAllele(numAlleles), mAllele(numAlleles);
  dMat fProbs(numAlleles), mProbs(numAlleles);
  
  
  // loop over all loci
  index=0;
  for(size_t i=0; i<numAlleles; ++i, index+=2){
    
    // FEMALES
    if( fScore[i] ){
      // homing
      // loop over alleles, all loci diploid
      for(size_t j=0; j<2; ++j){
        // fill allele and probs
        if(fGen[index+j] == 'W'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_homing_allele(i,0).begin(),
                            reference::instance().get_homing_allele(i,0).end());
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_homing_probs(i,0).begin(),
                           reference::instance().get_homing_probs(i,0).end());
        } else if(fGen[index+j] == 'H'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_homing_allele(i,1).begin(),
                            reference::instance().get_homing_allele(i,1).end());
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_homing_probs(i,1).begin(),
                           reference::instance().get_homing_probs(i,1).end());
        } else if(fGen[index+j] == 'R'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_homing_allele(i,2).begin(),
                            reference::instance().get_homing_allele(i,2).end());
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_homing_probs(i,2).begin(),
                           reference::instance().get_homing_probs(i,2).end());
        } else if(fGen[index+j] == 'S'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_homing_allele(i,3).begin(),
                            reference::instance().get_homing_allele(i,3).end());
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_homing_probs(i,3).begin(),
                           reference::instance().get_homing_probs(i,3).end());
        }
        
      } // end allele loop
    } else {
      // no homing
      // loop over alleles, all loci diploid
      for(size_t j=0; j<2; ++j){
        // fill allele and probs
        if(fGen[index+j] == 'W'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_mendelian_allele(i,0).begin(),
                            reference::instance().get_mendelian_allele(i,0).end());
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_mendelian_probs(i,0).begin(),
                           reference::instance().get_mendelian_probs(i,0).end());
        } else if(fGen[index+j] == 'H'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_mendelian_allele(i,1).begin(),
                            reference::instance().get_mendelian_allele(i,1).end());
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_mendelian_probs(i,1).begin(),
                           reference::instance().get_mendelian_probs(i,1).end());
        } else if(fGen[index+j] == 'R'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_mendelian_allele(i,2).begin(),
                            reference::instance().get_mendelian_allele(i,2).end());
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_mendelian_probs(i,2).begin(),
                           reference::instance().get_mendelian_probs(i,2).end());
        } else if(fGen[index+j] == 'S'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_mendelian_allele(i,3).begin(),
                            reference::instance().get_mendelian_allele(i,3).end());
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_mendelian_probs(i,3).begin(),
                           reference::instance().get_mendelian_probs(i,3).end());
        }
        
      } // end allele loop 
    } // end females
    
    
    // MALES
    if( mScore[i] ){
      // homing
      // loop over alleles, all loci diploid
      for(size_t j=0; j<2; ++j){
        // fill allele and probs
        if(mGen[index+j] == 'W'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_homing_allele(i,0).begin(),
                            reference::instance().get_homing_allele(i,0).end());
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_homing_probs(i,0).begin(),
                           reference::instance().get_homing_probs(i,0).end());
        } else if(mGen[index+j] == 'H'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_homing_allele(i,1).begin(),
                            reference::instance().get_homing_allele(i,1).end());
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_homing_probs(i,1).begin(),
                           reference::instance().get_homing_probs(i,1).end());
        } else if(mGen[index+j] == 'R'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_homing_allele(i,2).begin(),
                            reference::instance().get_homing_allele(i,2).end());
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_homing_probs(i,2).begin(),
                           reference::instance().get_homing_probs(i,2).end());
        } else if(mGen[index+j] == 'S'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_homing_allele(i,3).begin(),
                            reference::instance().get_homing_allele(i,3).end());
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_homing_probs(i,3).begin(),
                           reference::instance().get_homing_probs(i,3).end());
        }
        
      } // end allele loop
    } else {
      // no homing
      // loop over alleles, all loci diploid
      for(size_t j=0; j<2; ++j){
        // fill allele and probs
        if(mGen[index+j] == 'W'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_mendelian_allele(i,0).begin(),
                            reference::instance().get_mendelian_allele(i,0).end());
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_mendelian_probs(i,0).begin(),
                           reference::instance().get_mendelian_probs(i,0).end());
        } else if(mGen[index+j] == 'H'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_mendelian_allele(i,1).begin(),
                            reference::instance().get_mendelian_allele(i,1).end());
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_mendelian_probs(i,1).begin(),
                           reference::instance().get_mendelian_probs(i,1).end());
        } else if(mGen[index+j] == 'R'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_mendelian_allele(i,2).begin(),
                            reference::instance().get_mendelian_allele(i,2).end());
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_mendelian_probs(i,2).begin(),
                           reference::instance().get_mendelian_probs(i,2).end());
        } else if(mGen[index+j] == 'S'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_mendelian_allele(i,3).begin(),
                            reference::instance().get_mendelian_allele(i,3).end());
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_mendelian_probs(i,3).begin(),
                           reference::instance().get_mendelian_probs(i,3).end());
        }
        
      } // end allele loop
    } // end males
    
  } // end loci loop
  
  /*****************************************************************************/
  // End Next-Gen alleles
  /*****************************************************************************/
  
  /*****************************************************************************/
  // All Combinations of Male/Female for each Loci
  /*****************************************************************************/
  std::string holdAllele;
  
  // Aggregate duplicates before storing
  std::unordered_map<std::string, double>   duplicates;
  std::unordered_map<std::string, double>::iterator value;
  
  
  // loop over alleles
  for(size_t i=0; i<numAlleles; ++i){
    
    // all combinations of female/male alleles, and probs
    for(size_t fem=0; fem<fAllele[i].size(); ++fem){
      for(size_t mal=0; mal<mAllele[i].size(); ++mal){
        // combine and sort allele
        holdAllele = fAllele[i][fem] + mAllele[i][mal];
        std::sort(holdAllele.begin(), holdAllele.end() );
        
        // Here we aggregate non-unique values
        // check if it is in map
        value = duplicates.find(holdAllele);
        if(value == duplicates.end()){
          // not in map, add it
          duplicates.insert(std::make_pair(holdAllele, fProbs[i][fem]*mProbs[i][mal]));
        } else {
          // is in map, combine with existing value
          value->second += fProbs[i][fem]*mProbs[i][mal];
        }
        
      } // end male loop
    } // end female loop
    
    // clear fAllele and fProbs for re-use
    fAllele[i].clear();
    fProbs[i].clear();
    
    // Take unique values in map, return to matrix
    for(auto elem : duplicates){
      fAllele[i].push_back(elem.first);
      fProbs[i].push_back(elem.second);
    }
    
    // clear map
    duplicates.clear();
  } // end 
  
  /*****************************************************************************/
  // End All Combinations of MAle/Female for Each Loci
  /*****************************************************************************/
  
  /*****************************************************************************/
  // Cartesian Product of All Loci
  /*****************************************************************************/
  
  sVec finalGenotypes;
  dVec finalProbs;
  
  sVec holdGens;
  dVec holdProbs;
  
  // fill first index
  for(index=0; index<fAllele[0].size(); ++index){
    finalGenotypes.push_back(fAllele[0][index]);
    finalProbs.push_back(fProbs[0][index]);
  }
  
  // loop over rest of them
  for(index=1; index<numAlleles; ++index){
    // loop over current elements in return, and elements in next vector
    for(size_t current=0; current<finalGenotypes.size(); ++current){
      for(size_t newOne=0; newOne<fAllele[index].size(); ++newOne){
        // combine strings and probabilities
        holdGens.push_back(finalGenotypes[current] + fAllele[index][newOne]);
        holdProbs.push_back(finalProbs[current] * fProbs[index][newOne]);
      } // end loop over new things
    } // end loop over current things
    
    // set returns
    finalGenotypes = holdGens;
    finalProbs = holdProbs; 
    
    // clear holders
    holdGens.clear();
    holdProbs.clear();
  } // end loop over loci
  
  /*****************************************************************************/
  // End Cartesian Product of All Loci
  /*****************************************************************************/
  
} // end function







/**************************************
 * oneLocus Generating Function
 **************************************/




