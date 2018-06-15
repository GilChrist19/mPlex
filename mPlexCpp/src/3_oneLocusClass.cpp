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
oneLocus::oneLocus(const int& patchID_,
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
  CreateMosquitoes2Loci(parameters::instance().get_larva_eq(patchID),
                         minAge, ageDist, aTypes, eggs);
  
  // larva
  larva.reserve(2*parameters::instance().get_larva_eq(patchID));
  minAge = parameters::instance().get_stage_time(0)+1;
  
  ageDist.clear();
  for(int power = minAge; power <= parameters::instance().get_stage_sum(1); ++power){
   ageDist.push_back(std::pow(1.0-parameters::instance().get_mu(1), power));
  }
  
  CreateMosquitoes2Loci(parameters::instance().get_larva_eq(patchID),
                         minAge, ageDist, aTypes, larva);
  
  // pupa
  pupa.reserve(2*parameters::instance().get_adult_pop_eq(patchID));
  
  // adults
  adult_male.reserve(2*parameters::instance().get_adult_pop_eq(patchID));
  adult_female.reserve(2*parameters::instance().get_adult_pop_eq(patchID));
  unmated_female.reserve(2*parameters::instance().get_adult_pop_eq(patchID));
  
  minAge = parameters::instance().get_stage_sum(2)+1;
  ageDist.assign(parameters::instance().get_stage_sum(3) - minAge +1,1);
  CreateMosquitoes2Loci(parameters::instance().get_adult_pop_eq(patchID)/2,
                         minAge, ageDist, aTypes, adult_male);
  CreateMosquitoes2Loci(parameters::instance().get_adult_pop_eq(patchID)/2,
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


oneLocus::~oneLocus(){};

/******************************************************************************
 * Default (compiler-generated) move semantics
 ******************************************************************************/
oneLocus::oneLocus(oneLocus&& d) = default;
oneLocus& oneLocus::operator=(oneLocus&& d) = default;

/******************************************************************************
 * Reset
 ******************************************************************************/
void oneLocus::reset_Patch(const Rcpp::ListOf<Rcpp::List>& aTypes){
  
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
  CreateMosquitoes2Loci(parameters::instance().get_larva_eq(patchID),
                          minAge, ageDist, aTypes, eggs);
  
  // larva
  larva.clear();
  minAge = parameters::instance().get_stage_time(0)+1;
  
  ageDist.clear();
  for(int power = minAge; power <= parameters::instance().get_stage_sum(1); ++power){
    ageDist.push_back(std::pow(1.0-parameters::instance().get_mu(1), power));
  }
  
  CreateMosquitoes2Loci(parameters::instance().get_larva_eq(patchID),
                          minAge, ageDist, aTypes, larva);
  
  // pupa
  pupa.clear();
  
  // adults
  adult_male.clear();
  adult_female.clear();
  unmated_female.clear();
  
  minAge = parameters::instance().get_stage_sum(2)+1;
  ageDist.assign(parameters::instance().get_stage_sum(3) - minAge +1,1);
  CreateMosquitoes2Loci(parameters::instance().get_adult_pop_eq(patchID)/2,
                          minAge, ageDist, aTypes, adult_male);
  CreateMosquitoes2Loci(parameters::instance().get_adult_pop_eq(patchID)/2,
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
void oneLocus::oneDay_layEggs(){
  
  
  Rcpp::Rcout<<"You chose the oneLocus function!"<<std::endl;
  
  iVec newEggs;
  
  for(auto female : adult_female){
    
    // calculate genotypes and probs of offspring
    MultiplexOffspring_oLocus(female.get_genotype(), female.get_mate());
    
    // pull eggs over offspring probs
    newEggs = prng::instance().get_rmultinom(parameters::instance().get_beta()
                                               * reference::instance().get_s(female.get_genotype()),
                                               holdProbs1);
    
    // create new eggs
    for(size_t eggIndex=0; eggIndex<newEggs.size(); ++eggIndex){
      for(size_t it=0; it<eggIndex; ++it){
        eggs.emplace_back(Mosquito(0, holdGens1[eggIndex]));
      } // end loop over number of eggs per genotype
    } // end loop over newEggs vector
  } // end loop over females
  
}

/******************************************************************************
 * FUNCTIONS
 ******************************************************************************/

/**************************************
 * SETUP
 **************************************/
void CreateMosquitoes2Loci(const int& numMos, const int& minAge, const dVec& ageDist,
                           const Rcpp::ListOf<Rcpp::List>& aTypes, popVec returnPop){
  
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

/**************************************
 * GENERATING FUNCTION
 **************************************/
void oneLocus::MultiplexOffspring_oLocus(const std::string& fGen, const std::string& mGen){

  
  
  
  // get number of alleles, this is reused every time and resets here
  numAlleles = fGen.size()/2;
  
  /*****************************************************************************/
  // Score Each Allele
  /*****************************************************************************/
  // these get reused, and set here
  fScore = false;
  mScore = false;
  index = 0;
  
  //loop over loci, separate alleles and score
  while(!fScore || !mScore){
    // female score
    if(fGen[index] == 'H') {fScore = true;}
    // male score
    if(mGen[index] == 'H') {mScore = true;}
    // increment index
    ++index;
  }
  /*****************************************************************************/
  // End Score Each Allele
  /*****************************************************************************/
  
  /*****************************************************************************/
  // Determine Next-Gen alleles
  /*****************************************************************************/
  // reset things that are reused
  index=0;
  
  holdGens1.clear();
  holdProbs1.clear();
  
  fAllele.clear();
  fProbs.clear();
  mAllele.clear();
  mProbs.clear();
  
  
  // FEMALES
  if(fScore){
    // TRUE - there is homing
    // loop over alleles, only two
    for(size_t allele=0; allele<2; ++allele, index+=numAlleles){
      // loop over each target in allele
      for(size_t i=0; i<numAlleles; ++i){
        // fill allele and probs
        if(fGen[i+index] == 'W'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_homing_allele(i,0).begin(),
                                 reference::instance().get_homing_allele(i,0).end());
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_homing_probs(i,0).begin(),
                           reference::instance().get_homing_probs(i,0).end());
        } else if(fGen[i+index] == 'H'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_homing_allele(i,1).begin(),
                            reference::instance().get_homing_allele(i,1).end());
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_homing_probs(i,1).begin(),
                           reference::instance().get_homing_probs(i,1).end());
        } else if(fGen[i+index] == 'R'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_homing_allele(i,2).begin(),
                            reference::instance().get_homing_allele(i,2).end());
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_homing_probs(i,2).begin(),
                           reference::instance().get_homing_probs(i,2).end());
        } else if(fGen[i+index] == 'S'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_homing_allele(i,3).begin(),
                            reference::instance().get_homing_allele(i,3).end());
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_homing_probs(i,3).begin(),
                           reference::instance().get_homing_probs(i,3).end());
        }
        
        /********************************
        // All Combinations of Target Sites
        *********************************/
        // make combinations
        for(size_t first=0; first < holdGens1.size(); ++first){
          for(size_t second=0; second < holdGens2.size(); ++second){
            holdGens3.push_back(holdGens1[first] + holdGens2[second]);
            holdProbs3.push_back(holdProbs1[first] * holdProbs2[second]);
          }
        }
         
        // set combinations in return matrix
        //  Is there a better way to set the first one?
        if(holdGens1.size() == 0){
          holdGens1 = holdGens2;
          holdProbs1 = holdProbs2;
        } else {
          holdGens1 = holdGens3;
          holdProbs1 = holdProbs3;
        }
        
        // clear hold values
        holdGens2.clear();
        holdProbs2.clear();
        holdGens3.clear();
        holdProbs3.clear();
        /********************************
        // End All Combinations of Target Sites
        *********************************/
      } // end loop over target sites at each allele
      
      /********************************
      // "Unlist"
      *********************************/
      // put finished set of alleles from each locus into a final vector, they 
      //  can't mix with each other
      //  There could be duplicates here, do I care?
      fAllele.insert(fAllele.end(), holdGens1.begin(), holdGens1.end());
      fProbs.insert(fProbs.end(), holdProbs1.begin(), holdProbs1.end());
      
      // clear holder
      holdGens1.clear();
      holdProbs1.clear();
      
      /********************************
      // End "Unlist"
      *********************************/
    } // end allele loop
    
  } else {
    // FALSE - there is not homing
    // loop over alleles, only two
    for(size_t allele=0; allele<2; ++allele, index+=numAlleles){
      // loop over each target in allele
      for(size_t i=0; i<numAlleles; ++i){
        // fill allele and probs
        if(fGen[i+index] == 'W'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_mendelian_allele(i,0).begin(),
                                 reference::instance().get_mendelian_allele(i,0).end());
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_mendelian_probs(i,0).begin(),
                                reference::instance().get_mendelian_probs(i,0).end());
        } else if(fGen[i+index] == 'H'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_mendelian_allele(i,1).begin(),
                                 reference::instance().get_mendelian_allele(i,1).end());
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_mendelian_probs(i,1).begin(),
                                reference::instance().get_mendelian_probs(i,1).end());
        } else if(fGen[i+index] == 'R'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_mendelian_allele(i,2).begin(),
                                 reference::instance().get_mendelian_allele(i,2).end());
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_mendelian_probs(i,2).begin(),
                                reference::instance().get_mendelian_probs(i,2).end());
        } else if(fGen[i+index] == 'S'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_mendelian_allele(i,3).begin(),
                                 reference::instance().get_mendelian_allele(i,3).end());
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_mendelian_probs(i,3).begin(),
                                reference::instance().get_mendelian_probs(i,3).end());
        }
        
        /********************************
        // All Combinations of Target Sites
        *********************************/
        // make combinations
        for(size_t first=0; first < holdGens1.size(); ++first){
          for(size_t second=0; second < holdGens2.size(); ++second){
            holdGens3.push_back(holdGens1[first] + holdGens2[second]);
            holdProbs3.push_back(holdProbs1[first] * holdProbs2[second]);
          }
        }
        
        // set combinations in return matrix
        //  Is there a better way to set the first one?
        if(holdGens1.size() == 0){
          holdGens1 = holdGens2;
          holdProbs1 = holdProbs2;
        } else {
          holdGens1 = holdGens3;
          holdProbs1 = holdProbs3;
        }
        
        // clear hold values
        holdGens2.clear();
        holdProbs2.clear();
        holdGens3.clear();
        holdProbs3.clear();
        /********************************
        // End All Combinations of Target Sites
        *********************************/
      } // end loop over target sites at each allele
      
      /********************************
      // "Unlist"
      *********************************/
      // put finished set of alleles from each locus into a final vector, they 
      //  can't mix with each other
      //  There could be duplicates here, do I care?
      fAllele.insert(fAllele.end(), holdGens1.begin(), holdGens1.end());
      fProbs.insert(fProbs.end(), holdProbs1.begin(), holdProbs1.end());
      
      // clear holder
      holdGens1.clear();
      holdProbs1.clear();
      
      /********************************
      // End "Unlist"
      *********************************/
    } // end allele loop
    
  } // end female if statement
  
  // reset index
  index=0;
  
  // MALES
  if(mScore){
    // TRUE - there is homing
    // loop over alleles, only two
    for(size_t allele=0; allele<2; ++allele, index+=numAlleles){
      // loop over each target in allele
      for(size_t i=0; i<numAlleles; ++i){
        // fill allele and probs
        if(mGen[i+index] == 'W'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_homing_allele(i,0).begin(),
                                 reference::instance().get_homing_allele(i,0).end());
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_homing_probs(i,0).begin(),
                                reference::instance().get_homing_probs(i,0).end());
        } else if(mGen[i+index] == 'H'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_homing_allele(i,1).begin(),
                                 reference::instance().get_homing_allele(i,1).end());
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_homing_probs(i,1).begin(),
                                reference::instance().get_homing_probs(i,1).end());
        } else if(mGen[i+index] == 'R'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_homing_allele(i,2).begin(),
                                 reference::instance().get_homing_allele(i,2).end());
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_homing_probs(i,2).begin(),
                                reference::instance().get_homing_probs(i,2).end());
        } else if(mGen[i+index] == 'S'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_homing_allele(i,3).begin(),
                                 reference::instance().get_homing_allele(i,3).end());
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_homing_probs(i,3).begin(),
                                reference::instance().get_homing_probs(i,3).end());
        }
        
        /********************************
        // All Combinations of Target Sites
        *********************************/
        // make combinations
        for(size_t first=0; first < holdGens1.size(); ++first){
          for(size_t second=0; second < holdGens2.size(); ++second){
            holdGens3.push_back(holdGens1[first] + holdGens2[second]);
            holdProbs3.push_back(holdProbs1[first] * holdProbs2[second]);
          }
        }
        
        // set combinations in return matrix
        //  Is there a better way to set the first one?
        if(holdGens1.size() == 0){
          holdGens1 = holdGens2;
          holdProbs1 = holdProbs2;
        } else {
          holdGens1 = holdGens3;
          holdProbs1 = holdProbs3;
        }
        
        // clear hold values
        holdGens2.clear();
        holdProbs2.clear();
        holdGens3.clear();
        holdProbs3.clear();
        /********************************
        // End All Combinations of Target Sites
        *********************************/
      } // end loop over target sites at each allele
      
      /********************************
      // "Unlist"
      *********************************/
      // put finished set of alleles from each locus into a final vector, they 
      //  can't mix with each other
      //  There could be duplicates here, do I care?
      mAllele.insert(mAllele.end(), holdGens1.begin(), holdGens1.end());
      mProbs.insert(mProbs.end(), holdProbs1.begin(), holdProbs1.end());
      
      // clear holder
      holdGens1.clear();
      holdProbs1.clear();
      
      /********************************
      // End "Unlist"
      *********************************/
    } // end allele loop
    
  } else {
    
    // FALSE - there is no homing
    // loop over alleles, only two
    for(size_t allele=0; allele<2; ++allele, index+=numAlleles){
      // loop over each target in allele
      for(size_t i=0; i<numAlleles; ++i){
        // fill allele and probs
        if(mGen[i+index] == 'W'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_mendelian_allele(i,0).begin(),
                                 reference::instance().get_mendelian_allele(i,0).end());
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_mendelian_probs(i,0).begin(),
                                reference::instance().get_mendelian_probs(i,0).end());
        } else if(mGen[i+index] == 'H'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_mendelian_allele(i,1).begin(),
                                 reference::instance().get_mendelian_allele(i,1).end());
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_mendelian_probs(i,1).begin(),
                                reference::instance().get_mendelian_probs(i,1).end());
        } else if(mGen[i+index] == 'R'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_mendelian_allele(i,2).begin(),
                                 reference::instance().get_mendelian_allele(i,2).end());
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_mendelian_probs(i,2).begin(),
                                reference::instance().get_mendelian_probs(i,2).end());
        } else if(mGen[i+index] == 'S'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_mendelian_allele(i,3).begin(),
                                 reference::instance().get_mendelian_allele(i,3).end());
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_mendelian_probs(i,3).begin(),
                                reference::instance().get_mendelian_probs(i,3).end());
        }
        
        /********************************
        // All Combinations of Target Sites
        *********************************/
        // make combinations
        for(size_t first=0; first < holdGens1.size(); ++first){
          for(size_t second=0; second < holdGens2.size(); ++second){
            holdGens3.push_back(holdGens1[first] + holdGens2[second]);
            holdProbs3.push_back(holdProbs1[first] * holdProbs2[second]);
          }
        }
        
        // set combinations in return matrix
        //  Is there a better way to set the first one?
        if(holdGens1.size() == 0){
          holdGens1 = holdGens2;
          holdProbs1 = holdProbs2;
        } else {
          holdGens1 = holdGens3;
          holdProbs1 = holdProbs3;
        }
        
        // clear hold values
        holdGens2.clear();
        holdProbs2.clear();
        holdGens3.clear();
        holdProbs3.clear();
        /********************************
        // End All Combinations of Target Sites
        *********************************/
      } // end loop over target sites at each allele
      
      /********************************
      // "Unlist"
      *********************************/
      // put finished set of alleles from each locus into a final vector, they 
      //  can't mix with each other
      //  There could be duplicates here, do I care?
      mAllele.insert(mAllele.end(), holdGens1.begin(), holdGens1.end());
      mProbs.insert(mProbs.end(), holdProbs1.begin(), holdProbs1.end());
      
      // clear holder
      holdGens1.clear();
      holdProbs1.clear();
      
      /********************************
      // End "Unlist"
      *********************************/
    } // end allele loop
    
  } // end male if statement
  /*****************************************************************************/
  // End Next-Gen alleles
  /*****************************************************************************/
   
  /*****************************************************************************/
  // Cartesian Product of male/female loci
  /*****************************************************************************/
  // These get reused, but two don't need cleared
  // holdAllele
  duplicates.clear();
  // std::unordered_map<std::string, double>::iterator value;
  
  
  // all combinations of female/male alleles and probs
  for(size_t fem=0; fem<fAllele.size(); ++fem){
    for(size_t mal=0; mal<mAllele.size(); ++mal){
      
      // combine alleles
      holdAllele = fAllele[fem] + mAllele[mal];
      
      // Here we aggregate non-unique values
      // check if it is in map
      value = duplicates.find(holdAllele);
      if(value == duplicates.end()){
        // not in map, add it
        duplicates.insert(std::make_pair(holdAllele, fProbs[fem] * mProbs[mal]));
      } else {
        // is in map, combine with existing value
        value->second += fProbs[fem]*mProbs[mal];
      }
      
    } // end male loop
  } //  end female loop
  
  // clear for reuse
  holdGens1.clear();
  holdProbs1.clear();
  
  // Take unique values in map, return to matrix
  for(auto elem : duplicates){
    holdGens1.push_back(elem.first);
    holdProbs1.push_back(elem.second);
  }
  
  /*****************************************************************************/
  // End Cartesian Product of male/female loci
  /*****************************************************************************/
  
} // end function








