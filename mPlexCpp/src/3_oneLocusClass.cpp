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
void MultiplexOffspring_oLocus(const std::string& fGen, const std::string& mGen){
  
  
  
  
  
  
  
  // get number of alleles
  int numAlleles = fGen.size()/2;
  
  /*****************************************************************************/
  // Score Each Allele
  /*****************************************************************************/
  bool fScore = false, mScore = false;
  int index;
  
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
  sMat fAllele(2), mAllele(2);
  dMat fProbs(2), mProbs(2);
  index=0;
  
  // FEMALES
  if(fScore){
    // TRUE - there is homing
    // loop over alleles, only two
    for(size_t allele=0; allele<2; ++allele, index+=numAlleles){
      // loop over each target in allele
      for(size_t i=0; i<numAlleles; ++i){
        // fill allele and probs
        if(fGen[i+index] == 'W'){
          fAllele[allele].insert(fAllele[allele].end(), reference::instance().get_homing_allele(i,0).begin(),
                            reference::instance().get_homing_allele(i,0).end());
          fProbs[allele].insert(fProbs[allele].end(), reference::instance().get_homing_probs(i,0).begin(),
                           reference::instance().get_homing_probs(i,0).end());
        } else if(fGen[i+index] == 'H'){
          fAllele[allele].insert(fAllele[allele].end(), reference::instance().get_homing_allele(i,1).begin(),
                            reference::instance().get_homing_allele(i,1).end());
          fProbs[allele].insert(fProbs[allele].end(), reference::instance().get_homing_probs(i,1).begin(),
                           reference::instance().get_homing_probs(i,1).end());
        } else if(fGen[i+index] == 'R'){
          fAllele[allele].insert(fAllele[allele].end(), reference::instance().get_homing_allele(i,2).begin(),
                            reference::instance().get_homing_allele(i,2).end());
          fProbs[allele].insert(fProbs[allele].end(), reference::instance().get_homing_probs(i,2).begin(),
                           reference::instance().get_homing_probs(i,2).end());
        } else if(fGen[i+index] == 'S'){
          fAllele[allele].insert(fAllele[allele].end(), reference::instance().get_homing_allele(i,3).begin(),
                            reference::instance().get_homing_allele(i,3).end());
          fProbs[allele].insert(fProbs[allele].end(), reference::instance().get_homing_probs(i,3).begin(),
                           reference::instance().get_homing_probs(i,3).end());
        }
      } // end loop over target sites at each allele
    } // end allele loop
    
  } else {
    // FALSE - there is not homing
    // loop over alleles, only two
    for(size_t allele=0; allele<2; ++allele, index+=numAlleles){
      // loop over each target in allele
      for(size_t i=0; i<numAlleles; ++i){
        // fill allele and probs
        if(fGen[i+index] == 'W'){
          fAllele[allele].insert(fAllele[allele].end(), reference::instance().get_mendelian_allele(i,0).begin(),
                                 reference::instance().get_mendelian_allele(i,0).end());
          fProbs[allele].insert(fProbs[allele].end(), reference::instance().get_mendelian_probs(i,0).begin(),
                                reference::instance().get_mendelian_probs(i,0).end());
        } else if(fGen[i+index] == 'H'){
          fAllele[allele].insert(fAllele[allele].end(), reference::instance().get_mendelian_allele(i,1).begin(),
                                 reference::instance().get_mendelian_allele(i,1).end());
          fProbs[allele].insert(fProbs[allele].end(), reference::instance().get_mendelian_probs(i,1).begin(),
                                reference::instance().get_mendelian_probs(i,1).end());
        } else if(fGen[i+index] == 'R'){
          fAllele[allele].insert(fAllele[allele].end(), reference::instance().get_mendelian_allele(i,2).begin(),
                                 reference::instance().get_mendelian_allele(i,2).end());
          fProbs[allele].insert(fProbs[allele].end(), reference::instance().get_mendelian_probs(i,2).begin(),
                                reference::instance().get_mendelian_probs(i,2).end());
        } else if(fGen[i+index] == 'S'){
          fAllele[allele].insert(fAllele[allele].end(), reference::instance().get_mendelian_allele(i,3).begin(),
                                 reference::instance().get_mendelian_allele(i,3).end());
          fProbs[allele].insert(fProbs[allele].end(), reference::instance().get_mendelian_probs(i,3).begin(),
                                reference::instance().get_mendelian_probs(i,3).end());
        }
      } // end loop over target sites at each allele
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
          mAllele[allele].insert(mAllele[allele].end(), reference::instance().get_homing_allele(i,0).begin(),
                                 reference::instance().get_homing_allele(i,0).end());
          mProbs[allele].insert(mProbs[allele].end(), reference::instance().get_homing_probs(i,0).begin(),
                                reference::instance().get_homing_probs(i,0).end());
        } else if(mGen[i+index] == 'H'){
          mAllele[allele].insert(mAllele[allele].end(), reference::instance().get_homing_allele(i,1).begin(),
                                 reference::instance().get_homing_allele(i,1).end());
          mProbs[allele].insert(mProbs[allele].end(), reference::instance().get_homing_probs(i,1).begin(),
                                reference::instance().get_homing_probs(i,1).end());
        } else if(mGen[i+index] == 'R'){
          mAllele[allele].insert(mAllele[allele].end(), reference::instance().get_homing_allele(i,2).begin(),
                                 reference::instance().get_homing_allele(i,2).end());
          mProbs[allele].insert(mProbs[allele].end(), reference::instance().get_homing_probs(i,2).begin(),
                                reference::instance().get_homing_probs(i,2).end());
        } else if(mGen[i+index] == 'S'){
          mAllele[allele].insert(mAllele[allele].end(), reference::instance().get_homing_allele(i,3).begin(),
                                 reference::instance().get_homing_allele(i,3).end());
          mProbs[allele].insert(mProbs[allele].end(), reference::instance().get_homing_probs(i,3).begin(),
                                reference::instance().get_homing_probs(i,3).end());
        }
      } // end loop over target sites at each allele
    } // end allele loop
    
  } else {
    
    // FALSE - there is no homing
    // loop over alleles, only two
    for(size_t allele=0; allele<2; ++allele, index+=numAlleles){
      // loop over each target in allele
      for(size_t i=0; i<numAlleles; ++i){
        // fill allele and probs
        if(mGen[i+index] == 'W'){
          mAllele[allele].insert(mAllele[allele].end(), reference::instance().get_mendelian_allele(i,0).begin(),
                                 reference::instance().get_mendelian_allele(i,0).end());
          mProbs[allele].insert(mProbs[allele].end(), reference::instance().get_mendelian_probs(i,0).begin(),
                                reference::instance().get_mendelian_probs(i,0).end());
        } else if(mGen[i+index] == 'H'){
          mAllele[allele].insert(mAllele[allele].end(), reference::instance().get_mendelian_allele(i,1).begin(),
                                 reference::instance().get_mendelian_allele(i,1).end());
          mProbs[allele].insert(mProbs[allele].end(), reference::instance().get_mendelian_probs(i,1).begin(),
                                reference::instance().get_mendelian_probs(i,1).end());
        } else if(mGen[i+index] == 'R'){
          mAllele[allele].insert(mAllele[allele].end(), reference::instance().get_mendelian_allele(i,2).begin(),
                                 reference::instance().get_mendelian_allele(i,2).end());
          mProbs[allele].insert(mProbs[allele].end(), reference::instance().get_mendelian_probs(i,2).begin(),
                                reference::instance().get_mendelian_probs(i,2).end());
        } else if(mGen[i+index] == 'S'){
          mAllele[allele].insert(mAllele[allele].end(), reference::instance().get_mendelian_allele(i,3).begin(),
                                 reference::instance().get_mendelian_allele(i,3).end());
          mProbs[allele].insert(mProbs[allele].end(), reference::instance().get_mendelian_probs(i,3).begin(),
                                reference::instance().get_mendelian_probs(i,3).end());
        }
      } // end loop over target sites at each allele
    } // end allele loop
    
  } // end male if statement
  /*****************************************************************************/
  // End Next-Gen alleles
  /*****************************************************************************/
  
  /*****************************************************************************/
  // All Combinations of Target Sites for Each Allele
  /*****************************************************************************/
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
} // end function



