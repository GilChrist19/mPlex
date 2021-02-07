///////////////////////////////////////////////////////////////////////////////
//                           ____  __          ______          
//                ____ ___  / __ \/ /__  _  __/ ____/___  ____ 
//               / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
//              / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
//             /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/ 
//                                              /_/   /_/      
///////////////////////////////////////////////////////////////////////////////

#include "2_Patch.hpp"


/******************************************************************************
* constructor & destructor
******************************************************************************/
Patch::Patch(const int& patchID_,
             const Rcpp::List& maleReleases_,
             const Rcpp::List& femaleReleases_,
             const Rcpp::List& eggReleases_,
             const Rcpp::List& matedFemaleReleases_)
{
  
  // patch ID
  patchID = patchID_;
  
  /****************
   * RELEASES
   ****************/
  
  size_t mR;
  
  // male releases
  if(maleReleases_.size()>0){
    mR = maleReleases_.size();
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
    mR = femaleReleases_.size();
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
  if(eggReleases_.size()>0){
    mR = eggReleases_.size();
    releaseE.reserve(mR);
    for(size_t i=0; i<mR; i++){
      releaseE.emplace_back(release_event(Rcpp::as<Rcpp::List>(eggReleases_[i])["genVec"],
                                          Rcpp::as<Rcpp::List>(eggReleases_[i])["ageVec"],
                                          Rcpp::as<Rcpp::List>(eggReleases_[i])["tRelease"]
      ));
    }
    std::sort(releaseE.begin(), releaseE.end(), [](release_event a, release_event b){
      return a.release_time > b.release_time;
    });
  }
  
  // mated female releases
  if(matedFemaleReleases_.size()>0){
    mR = matedFemaleReleases_.size();
    releaseMF.reserve(mR);
    for(size_t i=0; i<mR; i++){
      releaseMF.emplace_back(release_mate(Rcpp::as<Rcpp::List>(matedFemaleReleases_[i])["genVec"],
                                          Rcpp::as<Rcpp::List>(matedFemaleReleases_[i])["mateVec"],
                                          Rcpp::as<Rcpp::List>(matedFemaleReleases_[i])["ageVec"],
                                          Rcpp::as<Rcpp::List>(matedFemaleReleases_[i])["tRelease"]
      ));
    }
    std::sort(releaseMF.begin(), releaseMF.end(), [](release_mate a, release_mate b){
      return a.release_time > b.release_time;
    });
  }
  
  // Things to hold for reset
  releaseM0 = releaseM;
  releaseF0 = releaseF;
  releaseE0 = releaseE;
  releaseMF0 = releaseMF;
  
  // set migration size objects
  maleMigration.resize(parameters::instance().get_n_patch());
  femaleMigration.resize(parameters::instance().get_n_patch());
  
  // set mating objects
  genNames.reserve(2*parameters::instance().get_adult_pop_eq(patchID));
  genProbs.reserve(2*parameters::instance().get_adult_pop_eq(patchID));

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
void Patch::oneDay_popDynamics(prng& myPRNG){

  // Death and Age
  oneDay_eggDeathAge(myPRNG);
  oneDay_larvaDeathAge(myPRNG);
  oneDay_pupaDeathAge(myPRNG);
  oneDay_adultDeathAge(myPRNG);
  
  // Mature
  oneDay_pupaMaturation(myPRNG);
  oneDay_larvaMaturation();
  oneDay_eggMaturation();
  
  // Mate
  oneDay_mating(myPRNG);
  
  // New Eggs
  oneDay_layEggs(myPRNG);
  
  // Releases
  oneDay_Releases();
  
  // Migration out
  oneDay_migrationOut(myPRNG);

};

/******************************************************************************
* Functions
******************************************************************************/
// cute swap/pop statement came from here
// https://gamedev.stackexchange.com/questions/33888/what-is-the-most-efficient-container-to-store-dynamic-game-objects-in
// 20210207 - JB
//  Updating swap operation.
//  Reality is, we are removing one object from the list. 
//  Currently, we swap the desired object with the end object, then remove it.
//  But, we could just overwrite the desired object, so we now have 2 copies of 
//   the end object, and then remove the end object, thereby releiving the extra copy.
//  This is less operations than a swap, and ends with the same effect.
/**************************************
 * Death
***************************************/
void Patch::oneDay_eggDeathAge(prng& myPRNG){

  // set bernoulli
  myPRNG.set_cBern(parameters::instance().get_mu(0));
  
  // Loop over all eggs in the vector
  for(auto it = eggs.rbegin(); it != eggs.rend(); ++it){
    // If it is your time, swap and remove
    if(myPRNG.get_cBern() ){
      // std::swap(*it, eggs.back());
      *it = std::move(eggs.back());
      eggs.pop_back();
    } else {
      it->age_one_day();
    }
  } // end loop
  
}

void Patch::oneDay_larvaDeathAge(prng& myPRNG){
  
  holdInt = parameters::instance().get_stage_time(1);
  holdDbl = parameters::instance().get_alpha(patchID);
  holdDbl = std::pow(holdDbl/(holdDbl + larva.size()), 1.0/holdInt );
  holdDbl = 1.0-holdDbl*(1.0-parameters::instance().get_mu(1));
  
  // set bernoulli
  myPRNG.set_cBern(holdDbl);
  

  // Loop over all larva in the vector
  for(auto it = larva.rbegin(); it != larva.rend(); ++it){
    // If it is your time, swap and remove
    if(myPRNG.get_cBern() ){
      // std::swap(*it, larva.back());
      *it = std::move(larva.back());
      larva.pop_back();
    } else {
      it->age_one_day();
    }
  } // end loop
  
}

void Patch::oneDay_pupaDeathAge(prng& myPRNG){
  
  // set bernoulli
  myPRNG.set_cBern(parameters::instance().get_mu(2) );

  // Loop over all pupa in the vector
  for(auto it = pupa.rbegin(); it != pupa.rend(); ++it){
    // If it is your time, swap and remove
    if(myPRNG.get_cBern() ){
      // std::swap(*it, pupa.back());
      *it = std::move(pupa.back());
      pupa.pop_back();
    } else {
      it->age_one_day();
    }
  } // end loop
  
}

void Patch::oneDay_adultDeathAge(prng& myPRNG){
  
  double probs;
  holdDbl = 1.0 - parameters::instance().get_mu(3);
  
  // Loop over all adult males in the vector
  for(auto it = adult_male.rbegin(); it != adult_male.rend(); ++it){
    // get genotype specific chance of death
    probs = holdDbl * reference::instance().get_omega(it->get_genotype());
    // If it is your time, swap and remove
    if(myPRNG.get_rBern(1.0-probs) ){
      // std::swap(*it, adult_male.back());
      *it = std::move(adult_male.back());
      adult_male.pop_back();
    } else {
      it->age_one_day();
    }
  } // end loop
  
  // Loop over all adult females in the vector
  for(auto it = adult_female.rbegin(); it != adult_female.rend(); ++it){
    // get genotype specific chance of death
    probs = holdDbl * reference::instance().get_omega(it->get_genotype());
    // If it is your time, swap and remove
    if(myPRNG.get_rBern(1.0-probs) ){
      // std::swap(*it, adult_female.back());
      *it = std::move(adult_female.back());
      adult_female.pop_back();
    } else {
      it->age_one_day();
    }
  } // end loop
  
}

/**************************************
 * Mature
***************************************/
void Patch::oneDay_pupaMaturation(prng& myPRNG){
  
  holdInt = parameters::instance().get_stage_sum(2);
  
  // set bernoulli
  myPRNG.set_cBern(parameters::instance().get_mu(3));
  
  // Loop over all pupa in the vector
  for(auto it = pupa.rbegin(); it != pupa.rend(); ++it){
    
    // If it is your time, we move in
    if(it->get_age() >= holdInt){
      
      // one extra death probability, in the math
      if(myPRNG.get_cBern() ){
        // death!
        // std::swap(*it, pupa.back());
        *it = std::move(pupa.back());
        pupa.pop_back();
        
        // gender mill, mwuahahahaha
      } else if(myPRNG.get_rBern(reference::instance().get_phi(it->get_genotype())) ){
        // female!
        // check if it actually makes it
        if(myPRNG.get_rBern(reference::instance().get_xiF(it->get_genotype())) ){
          // you become an unmated adult female
          unmated_female.push_back(*it);
        }
        // you actually died
        // std::swap(*it, pupa.back());
        *it = std::move(pupa.back());
        pupa.pop_back();
      } else {
        // male!
        // see if you actually make it
        if(myPRNG.get_rBern(reference::instance().get_xiM(it->get_genotype())) ){
          // you become and adult male
          adult_male.push_back(*it);
        }
        // you died
        // std::swap(*it, pupa.back());
        *it = std::move(pupa.back());
        pupa.pop_back();
        
      } // end gender mill
      
    } // end age statement
    
  } // end loop
  
}

void Patch::oneDay_larvaMaturation(){
  
  holdInt = parameters::instance().get_stage_sum(1);
  
  // Loop over all larva in the vector
  for(auto it = larva.rbegin(); it != larva.rend(); ++it){
    // If it is your time, swap and remove
    if(it->get_age() >= holdInt){
      // put maturing egg into larva
      pupa.push_back(*it);
      // swap position and remove egg
      // std::swap(*it, larva.back());
      *it = std::move(larva.back());
      larva.pop_back();
    }
  } // end loop
  
}

void Patch::oneDay_eggMaturation(){
  
  holdInt = parameters::instance().get_stage_sum(0);
  
  // Loop over all eggs in the vector
  for(auto it = eggs.rbegin(); it != eggs.rend(); ++it){
    // If it is your time, swap and remove
    if(it->get_age() >= holdInt){
      // put maturing egg into larva
      larva.push_back(*it);
      // swap position and remove egg
      // std::swap(*it, eggs.back());
      *it = std::move(eggs.back());
      eggs.pop_back();
    }
  } // end loop

}

/**************************************
 * Mate
***************************************/
void Patch::oneDay_mating(prng& myPRNG){
  
  if(!adult_male.empty() && !unmated_female.empty() ){
    
    // clear holder objects for reuse
    genNames.clear();
    genProbs.clear();

    for(auto& male : adult_male){
      genNames.push_back(male.get_genotype());
      genProbs.push_back(reference::instance().get_eta(male.get_genotype()));
    }

    // set sample
    myPRNG.set_oneSample(genProbs);

    // loop over unmated females
    for(auto it = unmated_female.rbegin(); it != unmated_female.rend(); ++it){
      //get and set mate
      mateName = genNames[myPRNG.get_oneSample()];
      it->set_mate(mateName);
      // move to real adult females and remove
      adult_female.push_back(*it);
      unmated_female.pop_back();
    } // end loop over unmated females

    // only do if there are unmated females and males to mate with
  } else if(adult_male.empty() && !unmated_female.empty() ){
    
    double probs;
    holdDbl = 1.0 - parameters::instance().get_mu(3);
    
    // Loop over all unmated females in the vector
    for(auto it = unmated_female.rbegin(); it != unmated_female.rend(); ++it){
      // get genotype specific chance of death
      probs = holdDbl * reference::instance().get_omega(it->get_genotype());
      // If it is your time, swap and remove
      if(myPRNG.get_rBern(1.0-probs) ){
        // std::swap(*it, unmated_female.back());
        *it = std::move(unmated_female.back());
        unmated_female.pop_back();
      } else {
        it->age_one_day();
      }
    } // end loop
    
  } // age/kill them for being alive, but don't move to full adults, still unmated
  
}

/**************************************
 * Releases
***************************************/
void Patch::oneDay_Releases(){
  
  /****************
   * MALE
  ****************/
  // if there are releases left, and the time is appropriate
  if( (!releaseM.empty()) && (releaseM.back().release_time <= parameters::instance().get_t_now()) ){
    
    for(size_t it = 0; it < releaseM.back().pop_ages.size(); ++it){
      
      adult_male.push_back(Mosquito(releaseM.back().pop_ages[it],
                                    releaseM.back().pop_names[it]) );
    } // end loop over released mosquitoes
    
    // remove release
    releaseM.pop_back();
  } // end male release
  
  /****************
   * FEMALE
  ****************/
  if( (!releaseF.empty()) && (releaseF.back().release_time <= parameters::instance().get_t_now()) ){
    
    for(size_t it = 0; it < releaseF.back().pop_ages.size(); ++it){
      
      unmated_female.push_back(Mosquito(releaseF.back().pop_ages[it],
                                        releaseF.back().pop_names[it]) );
    } // end loop over releases
    
    // remove release
    releaseF.pop_back();
  } // end female release
  
  /****************
   * Egg release
  ****************/
  if( (!releaseE.empty()) && (releaseE.back().release_time <= parameters::instance().get_t_now()) ){
    
    for(size_t it = 0; it < releaseE.back().pop_ages.size(); ++it){
      
      eggs.push_back(Mosquito(releaseE.back().pop_ages[it],
                              releaseE.back().pop_names[it]) );
    } // end loop over releases
    
    // remove release
    releaseE.pop_back();
  } // end larva release
  
  /****************
   * Mated Female release
  ****************/
  if( (!releaseMF.empty()) && (releaseMF.back().release_time <= parameters::instance().get_t_now()) ){
    
    for(size_t it = 0; it < releaseMF.back().pop_ages.size(); ++it){
      
      adult_female.push_back(Mosquito(releaseMF.back().pop_ages[it],
                                      releaseMF.back().pop_names[it],
                                      releaseMF.back().mate_names[it]) );
    } // end loop over releases
    
    // remove release
    releaseMF.pop_back();
  } // end larva release
  
}

/**************************************
 * Migration
***************************************/
void Patch::oneDay_migrationOut(prng& myPRNG){
  
  // clear old migration
  for(size_t i=0; i < parameters::instance().get_n_patch(); ++i){
    maleMigration[i].clear();
    femaleMigration[i].clear();
  }

  
  /****************
   * MALE
  ****************/
  myPRNG.set_oneSample(parameters::instance().get_male_migration(patchID));
  
  // loop over all adult males
  for(auto it = adult_male.rbegin(); it != adult_male.rend(); ++it){
    //get which patch he goes to
    holdInt = myPRNG.get_oneSample();
    
    // if not this patch
    if(holdInt != patchID){ 
      // add to new patch it goes to
      maleMigration[holdInt].push_back(*it);
      // remove
      // std::swap(*it, adult_male.back());
      *it = std::move(adult_male.back());
      adult_male.pop_back();
    } // end if this patch
    
  } // end loop over adult males
  
  
  /****************
   * FEMALE
  ****************/
  // reset oneSample for female probs
  myPRNG.set_oneSample(parameters::instance().get_female_migration(patchID));
  
  // loop over all females
  for(auto it = adult_female.rbegin(); it != adult_female.rend(); ++it){
    // get which patch she goes to
    holdInt = myPRNG.get_oneSample();
    
    // if not this patch
    if(holdInt != patchID){
      // add to new place
      femaleMigration[holdInt].push_back(*it);
      // remove
      // std::swap(*it, adult_female.back());
      *it = std::move(adult_female.back());
      adult_female.pop_back();
    } // end if
    
  } // end loop over females
  
  
  /****************
   * BATCH
  ****************/
  if(myPRNG.get_rBern(parameters::instance().get_batchProbs(patchID)) ){
    
    // which patch will they migrate to
    // set oneSample, then get it
    myPRNG.set_oneSample(parameters::instance().get_batchLocation(patchID));
    holdInt = myPRNG.get_oneSample();
    
    // set binomial for inside loop
    myPRNG.set_cBern(parameters::instance().get_batchMale(patchID));

    // do male migration
    for(auto mos = adult_male.rbegin(); mos != adult_male.rend(); ++mos){
      // see if it moves
      if(myPRNG.get_cBern() ){
        // put him in new patch migration set
        maleMigration[holdInt].push_back(*mos);
        // remove him from here
        // std::swap(*mos, adult_male.back());
        *mos = std::move(adult_male.back());
        adult_male.pop_back();
      }
    } // end males
    
    // reset binom for females
    myPRNG.set_cBern(parameters::instance().get_batchFemale(patchID));

    // female batch migration
    for(auto mos = adult_female.rbegin(); mos != adult_female.rend(); ++mos){
      // see if she moves
      if(myPRNG.get_cBern() ){
        // put her in the new patch
        femaleMigration[holdInt].push_back(*mos);
        // remove from current patch
        // std::swap(*mos, adult_female.back());
        *mos = std::move(adult_female.back());
        adult_female.pop_back();
      }
    } // end female batch

  } // end batch migration
  
}

void Patch::oneDay_migrationIn(const popVec& male, const popVec& female){
  // append new males
  adult_male.insert(adult_male.end(), male.begin(), male.end());
  // append females
  adult_female.insert(adult_female.end(), female.begin(), female.end() );
}

/**************************************
 * Print Functions
***************************************/
void Patch::init_output(std::ofstream& ADM_log, std::ofstream& ADF_log){
  
  // write headers
  ADM_log << "Time,Age,Genotype\n";
  ADF_log << "Time,Age,Genotype,Mate\n";
  
  // write males
  for(auto& mos : adult_male){
    ADM_log << parameters::instance().get_t_now() << "," << mos.print_male();
  }
  
  // write unmated females
  for(auto& mos : unmated_female){
    ADF_log << parameters::instance().get_t_now() << "," << mos.print_female();
  }
  
}


void Patch::oneDay_writeOutput(std::ofstream& ADM_log, std::ofstream& ADF_log){
  
  // write males
  if(!adult_male.empty()){
    for(auto& mos : adult_male){
      ADM_log << parameters::instance().get_t_now() << "," << mos.print_male();
    }
  } else {
    ADM_log << parameters::instance().get_t_now() << ",,\n";
  }

  
  // write adult females
  if(!adult_female.empty()){
    for(auto& mos : adult_female){
      ADF_log << parameters::instance().get_t_now() << "," << mos.print_female();
    }
  } else {
    ADF_log << parameters::instance().get_t_now() << ",,,\n";
  }

}

















