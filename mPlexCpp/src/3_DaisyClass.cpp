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



