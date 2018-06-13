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
 * Constructor & Destructor
******************************************************************************/
Daisy::Daisy(const int& patchID_,
             const Rcpp::ListOf<Rcpp::List>& aTypes,
             const Rcpp::ListOf<Rcpp::List>& maleReleases_,
             const Rcpp::ListOf<Rcpp::List>& femaleReleases_,
             const Rcpp::ListOf<Rcpp::List>& larvaeReleases_) : Patch::Patch()
{

  patchID = patchID_;
  
  /****************
  * SET POPULATIONS
  ****************/
  // holder objects, these are for distribution function
  int minAge;
  dVec ageDist;
  ageDist.reserve(parameters::instance().get_stage_sum(1));
  
  // eggs
  eggs.reserve(2*parameters::instance().get_larva_eq(patchID));
  minAge = 0;
  ageDist.assign(parameters::instance().get_stage_time(0)+1,1);
  CreateMosquitoes2Allele(parameters::instance().get_larva_eq(patchID),
                         minAge, ageDist, aTypes, eggs);
  
  // larva
  larva.reserve(2*parameters::instance().get_larva_eq(patchID));
  minAge = parameters::instance().get_stage_time(0)+1;
  
  ageDist.clear();
  for(int power = minAge; power <= parameters::instance().get_stage_sum(1); ++power){
   ageDist.push_back(std::pow(1.0-parameters::instance().get_mu(1), power));
  }
  
  CreateMosquitoes2Allele(parameters::instance().get_larva_eq(patchID),
                         minAge, ageDist, aTypes, larva);
  
  // pupa
  pupa.reserve(2*parameters::instance().get_adult_pop_eq(patchID));
  
  // adults
  adult_male.reserve(2*parameters::instance().get_adult_pop_eq(patchID));
  adult_female.reserve(2*parameters::instance().get_adult_pop_eq(patchID));
  unmated_female.reserve(2*parameters::instance().get_adult_pop_eq(patchID));
  
  minAge = parameters::instance().get_stage_sum(2)+1;
  ageDist.assign(parameters::instance().get_stage_sum(3) - minAge +1,1);
  CreateMosquitoes2Allele(parameters::instance().get_adult_pop_eq(patchID)/2,
                         minAge, ageDist, aTypes, adult_male);
  CreateMosquitoes2Allele(parameters::instance().get_adult_pop_eq(patchID)/2,
                         minAge, ageDist, aTypes, unmated_female);
  
  
  /****************
  * RELEASES
  ****************/
  
  // male releases
  if(maleReleases_.size()>0){
   size_t mR = maleReleases_.size();
   releaseM.reserve(mR);
   for(size_t i=0; i<mR; i++){
     releaseM[i] = release_event(Rcpp::as<sVec>(maleReleases_[i]["genVec"]),
                                 Rcpp::as<iVec>(maleReleases_[i]["ageVec"]),
                                 Rcpp::as<int>(maleReleases_[i]["tRelease"])
     );
   }
   std::sort(releaseM.begin(), releaseM.end(), [](release_event a, release_event b){
     return a.release_time > b.release_time;
   });
  }
  
  // female releases
  if(femaleReleases_.size()>0){
   size_t mR = femaleReleases_.size();
   releaseF.reserve(mR);
   for(size_t i=0; i<mR; i++){
     releaseF[i] = release_event(Rcpp::as<sVec>(femaleReleases_[i]["genVec"]),
                                 Rcpp::as<iVec>(femaleReleases_[i]["ageVec"]),
                                 Rcpp::as<int>(femaleReleases_[i]["tRelease"])
     );
   }
   std::sort(releaseF.begin(), releaseF.end(), [](release_event a, release_event b){
     return a.release_time > b.release_time;
   });
  }
  
  // larva releases
  if(larvaeReleases_.size()>0){
   size_t mR = larvaeReleases_.size();
   releaseL.reserve(mR);
   for(size_t i=0; i<mR; i++){
     releaseL[i] = release_event(Rcpp::as<sVec>(larvaeReleases_[i]["genVec"]),
                                 Rcpp::as<iVec>(larvaeReleases_[i]["ageVec"]),
                                 Rcpp::as<int>(larvaeReleases_[i]["tRelease"])
     );
   }
   std::sort(releaseL.begin(), releaseL.end(), [](release_event a, release_event b){
     return a.release_time > b.release_time;
   });
  }
  
  // Things to hold for reset
  releaseM0 = releaseM;
  releaseF0 = releaseF;
  releaseL0 = releaseL;
  
};

Daisy::~Daisy(){};

/******************************************************************************
 * Default (compiler-generated) move semantics
******************************************************************************/
Daisy::Daisy(Daisy&& d) = default;
Daisy& Daisy::operator=(Daisy&& d) = default;

/******************************************************************************
 * Reset
******************************************************************************/
void Daisy::reset_Patch(const Rcpp::ListOf<Rcpp::List>& aTypes){
  
  /****************
   * RESET POPULATIONS
   ****************/
  // holder objects, these are for distribution function
  int minAge;
  dVec ageDist;
  ageDist.reserve(parameters::instance().get_stage_sum(1));
  
  // eggs
  eggs.clear();
  minAge = 0;
  ageDist.assign(parameters::instance().get_stage_time(0)+1,1);
  CreateMosquitoes2Allele(parameters::instance().get_larva_eq(patchID),
                          minAge, ageDist, aTypes, eggs);
  
  // larva
  larva.clear();
  minAge = parameters::instance().get_stage_time(0)+1;
  
  ageDist.clear();
  for(int power = minAge; power <= parameters::instance().get_stage_sum(1); ++power){
    ageDist.push_back(std::pow(1.0-parameters::instance().get_mu(1), power));
  }
  
  CreateMosquitoes2Allele(parameters::instance().get_larva_eq(patchID),
                          minAge, ageDist, aTypes, larva);
  
  // pupa
  pupa.clear();
  
  // adults
  adult_male.clear();
  adult_female.clear();
  unmated_female.clear();
  
  minAge = parameters::instance().get_stage_sum(2)+1;
  ageDist.assign(parameters::instance().get_stage_sum(3) - minAge +1,1);
  CreateMosquitoes2Allele(parameters::instance().get_adult_pop_eq(patchID)/2,
                          minAge, ageDist, aTypes, adult_male);
  CreateMosquitoes2Allele(parameters::instance().get_adult_pop_eq(patchID)/2,
                          minAge, ageDist, aTypes, unmated_female);
  
  
  /****************
   * RESET RELEASES
   ****************/
  releaseM = releaseM0;
  releaseF = releaseF0;
  releaseL = releaseL0;
}


/******************************************************************************
 * Lay eggs
******************************************************************************/
void Daisy::oneDay_layEggs(){
  
  
  Rcpp::Rcout<<"You chose the daisy function!"<<std::endl;
  
  
  
}

/******************************************************************************
 * FUNCTIONS
 ******************************************************************************/

/**************************************
 * SETUP
 **************************************/
  // This is also used in the multiLocusClass setup and reset functions
void CreateMosquitoes2Allele(const int& numMos, const int& minAge, const dVec& ageDist,
                             const Rcpp::ListOf<Rcpp::List>& aTypes, popVec returnPop){
  
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

/**************************************
 * GENERATING FUNCTION
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



