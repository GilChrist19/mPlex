///////////////////////////////////////////////////////////////////////////////
//                           ____  __          ______          
//                ____ ___  / __ \/ /__  _  __/ ____/___  ____ 
//               / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
//              / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
//             /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/ 
//                                              /_/   /_/      
///////////////////////////////////////////////////////////////////////////////

#include "3_DriveDefinitions.hpp"



/******************************************************************************
 * Constructor & Destructor
 ******************************************************************************/
multiLocus::multiLocus(const int& patchID_,
                       const Rcpp::ListOf<Rcpp::List>& aTypes,
                       const Rcpp::List& maleReleases_,
                       const Rcpp::List& femaleReleases_,
                       const Rcpp::List& larvaeReleases_) : Patch::Patch()
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
       releaseM.emplace_back(release_event(Rcpp::as<Rcpp::List>(maleReleases_[i])["genVec"],
                                           Rcpp::as<Rcpp::List>(maleReleases_[i])["ageVec"],
                                           Rcpp::as<Rcpp::List>(maleReleases_[i])["tRelease"]
       ));
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
       releaseF.emplace_back(release_event(Rcpp::as<Rcpp::List>(femaleReleases_[i])["genVec"],
                                           Rcpp::as<Rcpp::List>(femaleReleases_[i])["ageVec"],
                                           Rcpp::as<Rcpp::List>(femaleReleases_[i])["tRelease"]
       ));
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
       releaseL.emplace_back(release_event(Rcpp::as<Rcpp::List>(larvaeReleases_[i])["genVec"],
                                           Rcpp::as<Rcpp::List>(larvaeReleases_[i])["ageVec"],
                                           Rcpp::as<Rcpp::List>(larvaeReleases_[i])["tRelease"]
       ));
     }
     std::sort(releaseL.begin(), releaseL.end(), [](release_event a, release_event b){
       return a.release_time > b.release_time;
     });
   }
   
   
   // Things to hold for reset
   releaseM0 = releaseM;
   releaseF0 = releaseF;
   releaseL0 = releaseL;
   
   // migration setup
   maleMigration.resize(parameters::instance().get_n_patch());
   femaleMigration.resize(parameters::instance().get_n_patch());
   
   // Reproduction setup
   numAlleles = aTypes.size();
   fProbs.resize(numAlleles);
   mProbs.resize(numAlleles);
   fAllele.resize(numAlleles);
   mAllele.resize(numAlleles);

};

multiLocus::~multiLocus(){};

/******************************************************************************
 * Default (compiler-generated) move semantics
 ******************************************************************************/
multiLocus::multiLocus(multiLocus&& d) = default;
multiLocus& multiLocus::operator=(multiLocus&& d) = default;

/******************************************************************************
 * Reset
 ******************************************************************************/
void multiLocus::reset_Patch(const Rcpp::ListOf<Rcpp::List>& aTypes){
  
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
void multiLocus::oneDay_layEggs(){
  
  for(auto female : adult_female){
    
    // calculate genotypes and probs of offspring
    MultiplexOffspring_mLoci(female.get_genotype(), female.get_mate());
    
    // pull eggs over offspring probs
    newEggs = prng::instance().get_rmultinom(parameters::instance().get_beta()
                                               * reference::instance().get_s(female.get_genotype()),
                                               finalProbs);
    
    // create new eggs
    for(size_t eggIndex=0; eggIndex<newEggs.size(); ++eggIndex){
      for(size_t it=0; it<newEggs[eggIndex]; ++it){
        eggs.emplace_back(Mosquito(0, finalGenotypes[eggIndex]));
      } // end loop over number of eggs per genotype
    } // end loop over newEggs vector
  } // end loop over females
  
}

/******************************************************************************
 * FUNCTIONS
 ******************************************************************************/

/**************************************
 * GENERATING FUNCTION
 **************************************/
void multiLocus::MultiplexOffspring_mLoci(const std::string& fGen, const std::string& mGen){



  //// list of objects I create every time ////
  // vector<bool> fScore
  // vector<bool> mScore






  // get number of alleles
  // don't need to reset because it's redefined every time
  // numAlleles = fGen.size()/2;

  /*****************************************************************************/
  // Score Each Allele
  /*****************************************************************************/
  std::vector<bool> fScore(numAlleles, false), mScore(numAlleles, false);

  //loop over loci, separate alleles and score
  index=0;
  for(size_t i=0; i<fGen.size(); i+=2, ++index){
    // female score
    if( (fGen[i] == 'H') || (fGen[i+1] == 'H')) {fScore[index] = true;}
    // male score
    if( (mGen[i] == 'H') || (mGen[i+1] == 'H')) {mScore[index] = true;}

  } // end scoring loop
  /*****************************************************************************/
  // End Score Each Allele
  /*****************************************************************************/

  /*****************************************************************************/
  // Determine Next-Gen alleles
  /*****************************************************************************/
  // these get reused, so clear them first
  for(index=0; index < numAlleles; ++index){
    fProbs[index].clear();
    mProbs[index].clear();
    fAllele[index].clear();
    mAllele[index].clear();
  }


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
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_homing_allele_begin(i,0),
                            reference::instance().get_homing_allele_end(i,0));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_homing_probs_begin(i,0),
                           reference::instance().get_homing_probs_end(i,0));
        } else if(fGen[index+j] == 'H'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_homing_allele_begin(i,1),
                            reference::instance().get_homing_allele_end(i,1));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_homing_probs_begin(i,1),
                           reference::instance().get_homing_probs_end(i,1));
        } else if(fGen[index+j] == 'R'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_homing_allele_begin(i,2),
                            reference::instance().get_homing_allele_end(i,2));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_homing_probs_begin(i,2),
                           reference::instance().get_homing_probs_end(i,2));
        } else if(fGen[index+j] == 'S'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_homing_allele_begin(i,3),
                            reference::instance().get_homing_allele_end(i,3));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_homing_probs_begin(i,3),
                           reference::instance().get_homing_probs_end(i,3));
        }

      } // end allele loop
    } else {
      // no homing
      // loop over alleles, all loci diploid
      for(size_t j=0; j<2; ++j){
        // fill allele and probs
        if(fGen[index+j] == 'W'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_mendelian_allele_begin(i,0),
                            reference::instance().get_mendelian_allele_end(i,0));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_mendelian_probs_begin(i,0),
                           reference::instance().get_mendelian_probs_end(i,0));
        } else if(fGen[index+j] == 'H'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_mendelian_allele_begin(i,1),
                            reference::instance().get_mendelian_allele_end(i,1));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_mendelian_probs_begin(i,1),
                           reference::instance().get_mendelian_probs_end(i,1));
        } else if(fGen[index+j] == 'R'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_mendelian_allele_begin(i,2),
                            reference::instance().get_mendelian_allele_end(i,2));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_mendelian_probs_begin(i,2),
                           reference::instance().get_mendelian_probs_end(i,2));
        } else if(fGen[index+j] == 'S'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_mendelian_allele_begin(i,3),
                            reference::instance().get_mendelian_allele_end(i,3));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_mendelian_probs_begin(i,3),
                           reference::instance().get_mendelian_probs_end(i,3));
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
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_homing_allele_begin(i,0),
                            reference::instance().get_homing_allele_end(i,0));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_homing_probs_begin(i,0),
                           reference::instance().get_homing_probs_end(i,0));
        } else if(mGen[index+j] == 'H'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_homing_allele_begin(i,1),
                            reference::instance().get_homing_allele_end(i,1));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_homing_probs_begin(i,1),
                           reference::instance().get_homing_probs_end(i,1));
        } else if(mGen[index+j] == 'R'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_homing_allele_begin(i,2),
                            reference::instance().get_homing_allele_end(i,2));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_homing_probs_begin(i,2),
                           reference::instance().get_homing_probs_end(i,2));
        } else if(mGen[index+j] == 'S'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_homing_allele_begin(i,3),
                            reference::instance().get_homing_allele_end(i,3));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_homing_probs_begin(i,3),
                           reference::instance().get_homing_probs_end(i,3));
        }

      } // end allele loop
    } else {
      // no homing
      // loop over alleles, all loci diploid
      for(size_t j=0; j<2; ++j){
        // fill allele and probs
        if(mGen[index+j] == 'W'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_mendelian_allele_begin(i,0),
                            reference::instance().get_mendelian_allele_end(i,0));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_mendelian_probs_begin(i,0),
                           reference::instance().get_mendelian_probs_end(i,0));
        } else if(mGen[index+j] == 'H'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_mendelian_allele_begin(i,1),
                            reference::instance().get_mendelian_allele_end(i,1));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_mendelian_probs_begin(i,1),
                           reference::instance().get_mendelian_probs_end(i,1));
        } else if(mGen[index+j] == 'R'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_mendelian_allele_begin(i,2),
                            reference::instance().get_mendelian_allele_end(i,2));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_mendelian_probs_begin(i,2),
                           reference::instance().get_mendelian_probs_end(i,2));
        } else if(mGen[index+j] == 'S'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_mendelian_allele_begin(i,3),
                            reference::instance().get_mendelian_allele_end(i,3));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_mendelian_probs_begin(i,3),
                           reference::instance().get_mendelian_probs_end(i,3));
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

  // Aggregate duplicates before storing
  // these get reused, so save them
  duplicates.clear();
  // value
  // holdAllele is defined in class, always gets reset


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
  // these get reused, clear them
  finalGenotypes.clear();
  finalProbs.clear();

  holdGens.clear();
  holdProbs.clear();

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










