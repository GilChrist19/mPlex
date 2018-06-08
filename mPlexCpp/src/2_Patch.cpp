//      __  _____________  ____  ______
//     /  |/  / ____/ __ \/ __ \/ ____/___  ____
//    / /|_/ / / __/ / / / / / / /   / __ \/ __ \
//   / /  / / /_/ / /_/ / /_/ / /___/ /_/ / /_/ /
//  /_/  /_/\____/_____/_____/\____/ .___/ .___/
//                                /_/   /_/

#include "1_Mosquito.hpp"
#include "2_Patch.hpp"
#include "mPlex-PRNG.hpp"
#include "mPlex-Parameters.hpp"
#include "mPlex-Reference.hpp"

#include <math.h>















/******************************************************************************
* constructor & destructor
******************************************************************************/
Patch::Patch(const int& patchID_,
             const int& simTime_,
             const popVec& eggs_t0_,
             const popVec& larva_t0_,
             const popVec& pupa_t0_,
             const popVec& adult_male_t0_,
             const popVec& unmated_female_t0_,
             const Rcpp::ListOf<Rcpp::List> maleReleases_,
             const Rcpp::ListOf<Rcpp::List> femaleReleases_,
             const Rcpp::ListOf<Rcpp::List> larvaeReleases_)
{



  // set populations
  eggs_t0 = eggs_t0_;
  larva_t0 = larva_t0_;
  pupa_t0 = pupa_t0_;
  adult_male_t0 = adult_male_t0_;
  unmated_female_t0 = unmated_female_t0_;

  eggs.reserve(2*eggs_t0_.size());
  eggs = eggs_t0_;
  
  larva.reserve(2*larva_t0_.size());
  larva = larva_t0_;
  
  pupa.reserve(pupa_t0_.size());
  pupa = pupa_t0_;
  
  adult_male.reserve(adult_male_t0_.size());
  adult_male = adult_male_t0_;
  
  adult_female.reserve(adult_male_t0_.size());
  
  unmated_female = unmated_female_t0_;
  
  
  
  

  
  
  
  
  
  
  
  
  
  
  // modified mosquito releases

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

Patch::~Patch(){};

/******************************************************************************
* default (compiler-generated) move semantics
******************************************************************************/
Patch::Patch(Patch&& rhs) = default;
Patch& Patch::operator=(Patch&& rhs) = default;






/******************************************************************************
* Population Dynamics for One Day
******************************************************************************/
/* population dynamics */
void Patch::oneDay_popDynamics(Patch* P){

  // Death
  oneDay_eggDeath(P);
  oneDay_larvaDeath(P);
  oneDay_pupaDeath(P);
  oneDay_adultDeath(P);
  
  // Age
  oneDay_eggAge(P);
  oneDay_larvaeAge(P);
  oneDay_pupaAge(P);
  oneDay_adultAge(P);
  
  // Mature
  oneDay_pupaMaturation(P);
  oneDay_larvaMaturation(P);
  oneDay_eggMaturation(P);
  
  // Mate
  oneDay_mating(P);
  
  // New Eggs
  oneDay_layEggs(P);
  
  // Releases
  oneDay_Releases(P);

};

/******************************************************************************
* Functions
******************************************************************************/
/**************************************
 * Death
***************************************/
void Patch::oneDay_eggDeath(Patch* P){

  // Loop over all eggs in the vector
  for(auto it = P->eggs.rbegin(); it != P->eggs.rend(); ++it){
    // If it is your time, swap and remove
    if(prng::instance().get_rbinom(1, parameters::instance().get_mu(0))){
      std::swap(*it, P->eggs.back());
      P->eggs.pop_back();
    }
  } // end loop
  
}

void Patch::oneDay_larvaDeath(Patch* P){
  
  int stage = parameters::instance().get_stage_time(1);
  double alpha = parameters::instance().get_alpha(P->patchID);
  alpha = std::pow(alpha/(alpha + P->larva.size()), 1.0/stage);
  alpha = 1.0-alpha*(1.0-parameters::instance().get_mu(1));
  
  // Loop over all larva in the vector
  for(auto it = P->larva.rbegin(); it != P->larva.rend(); ++it){
    // If it is your time, swap and remove
    if(prng::instance().get_rbinom(1, alpha) ){
      std::swap(*it, P->larva.back());
      P->larva.pop_back();
    }
  } // end loop

}

void Patch::oneDay_pupaDeath(Patch* P){
  
  // Loop over all pupa in the vector
  for(auto it = P->pupa.rbegin(); it != P->pupa.rend(); ++it){
    // If it is your time, swap and remove
    if(prng::instance().get_rbinom(1, parameters::instance().get_mu(2))){
      std::swap(*it, P->pupa.back());
      P->pupa.pop_back();
    }
  } // end loop
  
}

void Patch::oneDay_adultDeath(Patch* P){
  
  double probs;
  double minusMu = 1.0 - parameters::instance().get_mu(3);
  
  // Loop over all adult males in the vector
  for(auto it = P->adult_male.rbegin(); it != P->adult_male.rend(); ++it){
    // get genotype specific chance of death
    probs = minusMu * reference::instance().get_omega(it->get_genotype());
    // If it is your time, swap and remove
    if(prng::instance().get_rbinom(1, 1.0-probs )){
      std::swap(*it, P->adult_male.back());
      P->adult_male.pop_back();
    }
  } // end loop
  
  // Loop over all adult females in the vector
  for(auto it = P->adult_female.rbegin(); it != P->adult_female.rend(); ++it){
    // get genotype specific chance of death
    probs = minusMu * reference::instance().get_omega(it->get_genotype());
    // If it is your time, swap and remove
    if(prng::instance().get_rbinom(1, 1.0-probs )){
      std::swap(*it, P->adult_female.back());
      P->adult_female.pop_back();
    }
  } // end loop
  
}

/**************************************
 * Age
***************************************/
void Patch::oneDay_eggAge(Patch* P){
  // loop over eggs, age them all one day
  for(auto& it : P->eggs){
    it.age_one_day();
  }
}

void Patch::oneDay_larvaeAge(Patch* P){
  // loop over larva, age them all one day
  for(auto& it : P->larva){
    it.age_one_day();
  }
}

void Patch::oneDay_pupaAge(Patch* P){
  // loop over pupa, age them all one day
  for(auto& it : P->pupa){
    it.age_one_day();
  }
}

void Patch::oneDay_adultAge(Patch* P){
  // loop over males
  for(auto& it : P->adult_male){
    it.age_one_day();
  }
  // loop over females
  for(auto& it : P-> adult_female){
    it.age_one_day();
  }
}

/**************************************
 * Mature
***************************************/
void Patch::oneDay_pupaMaturation(Patch* P){
  
  int age = parameters::instance().get_stage_sum(2);
  
  // Loop over all pupa in the vector
  for(auto it = P->pupa.rbegin(); it != P->pupa.rend(); ++it){
    
    // If it is your time, we move in
    if(it->get_age() > age){
      
      // gender mill, mwuahahahaha
      if(prng::instance().get_rbinom(1, reference::instance().get_phi(it->get_genotype()) )){
        // female!
        // check if it actually makes it
        if(prng::instance().get_rbinom(1, reference::instance().get_xiF(it->get_genotype()) )){
          // you become an unmated adult female
          P->unmated_female.push_back(*it);
        }
        // you actually died
        std::swap(*it, P->pupa.back());
        P->pupa.pop_back();
      } else {
        // male!
        // see if you actually make it
        if(prng::instance().get_rbinom(1, reference::instance().get_xiM(it->get_genotype()) )){
          // you become and adult male
          P->adult_male.push_back(*it);
        }
        // you died
        std::swap(*it, P->pupa.back());
        P->pupa.pop_back();
        
      } // end gender mill
      
    } // end age statement
    
  } // end loop
  
}

void Patch::oneDay_larvaMaturation(Patch* P){
  
  int age = parameters::instance().get_stage_sum(1);
  
  // Loop over all larva in the vector
  for(auto it = P->larva.rbegin(); it != P->larva.rend(); ++it){
    // If it is your time, swap and remove
    if(it->get_age() > age){
      // put maturing egg into larva
      P->pupa.push_back(*it);
      // swap position and remove egg
      std::swap(*it, P->larva.back());
      P->larva.pop_back();
    }
  } // end loop
  
}

void Patch::oneDay_eggMaturation(Patch* P){
  
  int age = parameters::instance().get_stage_time(0);
  
  // Loop over all eggs in the vector
  for(auto it = P->eggs.rbegin(); it != P->eggs.rend(); ++it){
    // If it is your time, swap and remove
    if(it->get_age() > age){
      // put maturing egg into larva
      P->larva.push_back(*it);
      // swap position and remove egg
      std::swap(*it, P->eggs.back());
      P->eggs.pop_back();
    }
  } // end loop

}

/**************************************
 * Mate
***************************************/
void Patch::oneDay_mating(Patch* P){
  
  if((P->adult_male.size() != 0) && (P->unmated_female.size() !=0) ){
    
    // genotype and probs vectors
    std::vector<std::string> genotypes(P->adult_male.size());
    std::vector<double> probs(P->adult_male.size());
    std::string mate;
    
    // loop over males, get genotypes and mating fitness
    for(size_t it=0; it<P->adult_male.size(); ++it){
      genotypes[it] = P->adult_male[it].get_genotype();
      probs[it] = reference::instance().get_eta(genotypes[it]);
    }
    
    // loop over unmated females
    for(auto it= P->unmated_female.rbegin(); it != P->unmated_female.rend(); ++it){
      //get and set mate
      mate = genotypes[prng::instance().get_oneSample(probs)];
      it->set_mate(mate);
      // move to real adult females and remove
      P->adult_female.push_back(*it);
      P->unmated_female.pop_back();
    } // end loop over unmated females 

  } // only do if there are unmated females and males to mate with
  
}

/**************************************
 * Releases
***************************************/
void Patch::oneDay_Releases(Patch* P){
  
  /****************
   * MALE
  ****************/
  // if there are releases left, and the time is appropriate
  if( (!P->releaseM.empty()) && (P->releaseM.back().release_time <= parameters::instance().get_t_now()) ){
    
    for(size_t it = 0; it < P->releaseM.back().pop_ages.size(); ++it){
      
      P->adult_male.push_back(Mosquito(P->releaseM.back().pop_ages[it],
                                       P->releaseM.back().pop_names[it])
                              );
    } // end loop over released mosquitoes
    
    // remove release
    P->releaseM.pop_back();
  } // end male release
  
  /****************
   * FEMALE
  ****************/
  if( (!P->releaseF.empty()) && (P->releaseF.back().release_time <= parameters::instance().get_t_now()) ){
    
    
    for(size_t it = 0; it < P->releaseF.back().pop_ages.size(); ++it){
      
      P->adult_female.push_back(Mosquito(P->releaseF.back().pop_ages[it],
                                         P->releaseF.back().pop_names[it])
                                );
    } // end loop over releases
    
    // remove release
    P->releaseF.pop_back();
  } // end female release
  
  /****************
   * LARVA
  ****************/
  if( (!P->releaseL.empty()) && (P->releaseL.back().release_time <= parameters::instance().get_t_now()) ){
    
    
    for(size_t it = 0; it < P->releaseL.back().pop_ages.size(); ++it){
      
      P->adult_female.push_back(Mosquito(P->releaseL.back().pop_ages[it],
                                         P->releaseL.back().pop_names[it])
      );
    } // end loop over releases
    
    // remove release
    P->releaseL.pop_back();
  } // end larva release
  
}

/**************************************
 * Migration
***************************************/
void Patch::oneDay_migrationOut(Patch* P){
  
  // clear old migration
  P->maleMigration.clear();
  P->femaleMigration.clear();
  
  // used a lot
  int patch;
  dVec Probs = prng::instance().get_rdirichlet(parameters::instance().get_male_migration(P->patchID)); 
  
  /****************
   * MALE
  ****************/
  // loop over all adult males
  for(auto it = P->adult_male.rbegin(); it != P->adult_male.rend(); ++it){
    //get which patch he goes to
    patch = prng::instance().get_oneSample(Probs);
    
    // if not this patch
    if(patch != P->patchID){ 
      // add to new patch it goes to
      P->maleMigration[patch].push_back(*it);
      // remove
      std::swap(*it, P->adult_male.back());
      P->adult_male.pop_back();
    } // end if this patch
    
  } // end loop over adult males
  
  
  /****************
   * FEMALE
  ****************/
  // get slightly more variance in probability
  Probs = prng::instance().get_rdirichlet(parameters::instance().get_female_migration(P->patchID));
  
  // loop over all females
  for(auto it = P->adult_female.rbegin(); it != P->adult_female.rend(); ++it){
    // get which patch she goes to
    patch = prng::instance().get_oneSample(Probs);
    
    // if not this patch
    if(patch != P->patchID){
      // add to new place
      P->femaleMigration[patch].push_back(*it);
      // remove
      std::swap(*it, P->adult_female.back());
      P->adult_female.pop_back();
    } // end if
    
  } // end loop over females
  
}

void Patch::oneDay_migrationIn(const popVec& male, const popVec& female){
  // append new males
  adult_male.insert(adult_male.end(), male.begin(), male.end());
  // append females
  adult_female.insert(adult_female.end(), female.begin(), female.end() );
}

/**************************************
 * Reset
***************************************/
void Patch::reset_Patch(){
  
  // reset aquatic phases
  eggs = eggs_t0;
  larva = larva_t0;
  pupa = pupa_t0;
  
  // reset adults
  adult_male = adult_male_t0;
  unmated_female = unmated_female_t0;
  adult_female.clear();
  
  // reset releases
  releaseM = releaseM0;
  releaseF = releaseF0;
  releaseL = releaseL0;
}

/**************************************
 * Print Functions
***************************************/
void Patch::init_output(){
  
  
  
  
  
}


void Patch::oneDay_writeOutput(){
  
  
  
  
  
  
  
}

















