///////////////////////////////////////////////////////////////////////////////
//     ________ __ __  _______ 
//    / ____/ //_//  |/  / __ \
//   / /   / ,<  / /|_/ / /_/ /
//  / /___/ /| |/ /  / / _, _/ 
//  \____/_/ |_/_/  /_/_/ |_|  
//   
///////////////////////////////////////////////////////////////////////////////

#include "2_Patch.hpp"


/******************************************************************************
* constructor & destructor
******************************************************************************/
Patch::Patch(const int& patchID_,
             const Rcpp::List& maleReleases_,
             const Rcpp::List& femaleReleases_,
             const Rcpp::List& eggReleases_)
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
  
  
  // Things to hold for reset
  releaseM0 = releaseM;
  releaseF0 = releaseF;
  releaseE0 = releaseE;
  
  
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
void Patch::oneDay_popDynamics(prng& myPRNG, bigBrother& myBB,
                               std::vector<std::ofstream *>& myLogFiles){

  // Death and Age
  oneDay_eggDeathAge(myPRNG);
  oneDay_larvaDeathAge(myPRNG);
  oneDay_pupaDeathAge(myPRNG);
  oneDay_adultDeathAge(myPRNG);
  
  // Mature
  oneDay_pupaMaturation(myPRNG, *myLogFiles[0]);
  oneDay_larvaMaturation();
  oneDay_eggMaturation();
  
  // Mate
  oneDay_mating(myPRNG, *myLogFiles[1]);
  
  // New Eggs
  oneDay_layEggs(myPRNG, myBB);
  
  // Releases
  oneDay_Releases(*myLogFiles[0]);
  
  // Flush Buffer for Safety
  (*myLogFiles[0]).flush();
  (*myLogFiles[1]).flush();
  
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
//   the end object, and then remove the end object, thereby relieving the extra copy.
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
  int maxAge(parameters::instance().get_male_max_age()); // instantiate to male max age
  holdDbl = 1.0 - parameters::instance().get_mu(3);
  
  // Loop over all adult males in the vector
  for(auto it = adult_male.rbegin(); it != adult_male.rend(); ++it){
    // get myID specific chance of death
    probs = holdDbl * reference::instance().get_omega(it->get_myID());
    // If it is your time, swap and remove
    if(myPRNG.get_rBern(1.0-probs) || (it->get_age() >= maxAge)){
      // std::swap(*it, adult_male.back());
      *it = std::move(adult_male.back());
      adult_male.pop_back();
    } else {
      it->age_one_day();
    }
  } // end loop
  
  // grab female max age
  maxAge = parameters::instance().get_female_max_age();
  
  // Loop over all adult females in the vector
  for(auto it = adult_female.rbegin(); it != adult_female.rend(); ++it){
    // get myID specific chance of death
    probs = holdDbl * reference::instance().get_omega(it->get_myID());
    // If it is your time, swap and remove
    if(myPRNG.get_rBern(1.0-probs) || (it->get_age() >= maxAge)){
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
void Patch::oneDay_pupaMaturation(prng& myPRNG, std::ofstream& mFile){
  
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
      } else if(myPRNG.get_rBern(reference::instance().get_phi(it->get_myID())) ){
        // female!
        // check if it actually makes it
        if(myPRNG.get_rBern(reference::instance().get_xiF(it->get_myID())) ){
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
        if(myPRNG.get_rBern(reference::instance().get_xiM(it->get_myID())) ){
          // you become and adult male
          adult_male.push_back(*it);
          // log mosquito
          mFile << patchID << (*it).print_male_short();
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
void Patch::oneDay_mating(prng& myPRNG, std::ofstream& fFile){
  
  if(!adult_male.empty() && !unmated_female.empty() ){
    
    // clear holder objects for reuse
    genNames.clear();
    genProbs.clear();

    for(auto& male : adult_male){
      genNames.push_back(male.get_myID());
      genProbs.push_back(reference::instance().get_eta(male.get_myID()));
    }

    // set sample
    myPRNG.set_oneSample(genProbs);

    // loop over unmated females
    for(auto it = unmated_female.rbegin(); it != unmated_female.rend(); ++it){
      //get and set mate
      mateName = genNames[myPRNG.get_oneSample()];
      it->set_mate(mateName);
      // move to real adult females
      adult_female.push_back(*it);
      // log mosquito
      fFile << patchID << (*it).print_female_short();
      // remove
      unmated_female.pop_back();
    } // end loop over unmated females

    // only do if there are unmated females and males to mate with
  } else if(adult_male.empty() && !unmated_female.empty() ){
    
    double probs;
    int maxAge(parameters::instance().get_female_max_age()); // instantiate to female max age
    holdDbl = 1.0 - parameters::instance().get_mu(3);
    
    // Loop over all unmated females in the vector
    for(auto it = unmated_female.rbegin(); it != unmated_female.rend(); ++it){
      // get myID specific chance of death
      probs = holdDbl * reference::instance().get_omega(it->get_myID());
      // If it is your time, swap and remove
      if(myPRNG.get_rBern(1.0-probs) || (it->get_age() >= maxAge)){
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
void Patch::oneDay_Releases(std::ofstream& mFile){
  
  /****************
   * MALE
  ****************/
  // if there are releases left, and the time is appropriate
  if( (!releaseM.empty()) && (releaseM.back().release_time <= parameters::instance().get_t_now()) ){
    
    for(size_t it = 0; it < releaseM.back().pop_ages.size(); ++it){
      // make new mosquito
      adult_male.push_back(Mosquito(releaseM.back().pop_ages[it],
                                    releaseM.back().pop_names[it],
                                    "0", "0")
                           );
      
      // log mosquito
      mFile << patchID << adult_male.back().print_male_short();
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
                                      releaseF.back().pop_names[it],
                                      "0", "0")
                             );
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
                              releaseE.back().pop_names[it],
                              "0", "0")
                     );
    } // end loop over releases
    
    // remove release
    releaseE.pop_back();
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
  // set oneSample
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
void Patch::init_output(std::vector<std::ofstream *>& logFiles){
  
  // logFile has 6 ofstreams in it
  // Eggs, Larva, Pupa, Males, Females, poplog_counts
  
  // write headers for Eggs, Larva, Pupa, and Males
  for(size_t i : {0,1,2,3}){
    *logFiles[i] << "Time,Age,myID,momID,dadID\n";
  }
  
  // write header for Females
  *logFiles[4] << "Time,Age,myID,momID,dadID,Mate\n";
  
  // write header for poplog_counts
  *logFiles[5] << "Time,Eggs,Larvae,Pupae,Males,Females\n";
  
  // no population to write because we sample at specific times/days
}




// accessory function for writing mosquitoes
void write_stage(popVec& pop, const int& patchID, const int& stage,
                 std::ofstream& logFile, std::string (Mosquito::*print)(),
                 prng& myPRNG, const int& cTime){

  
  if( parameters::instance().get_sampDay(stage, patchID) ){
    // set bernoulli
    myPRNG.set_cBern(parameters::instance().get_sampCov(stage, patchID));
    // loop over everyone in current stage
    for(auto it = pop.rbegin(); it != pop.rend(); ++it){
      // if we sample it, it then must die
      //  otherwise, nothing happens
      if(myPRNG.get_cBern()){
        // write output
        logFile << cTime << ((*it).*print)();
        
        // kill
        // std::swap(*it, pop.back());
        *it = std::move(pop.back());
        pop.pop_back();
      } // end if statement
    } // end loop over critters
    
    // Flush Buffer for Safety
    logFile.flush();
    
  } // end if statement
} // end function







void Patch::oneDay_writeOutput(std::vector<std::ofstream *>& logFiles, prng& myPRNG){
  // grab current time
  int curTime = parameters::instance().get_t_now();
  
  // output counts
  //  do this first b/c we destructively sample the population below, so this 
  //  gives the true total
  *logFiles[5] << curTime << "," << eggs.size() << "," << larva.size() 
               << "," << pupa.size() << "," << adult_male.size() << "," 
               << adult_female.size() << std::endl;
  
  // check eggs
  write_stage(eggs, patchID, 0, *logFiles[0], &Mosquito::print_aquatic, myPRNG, curTime);
  
  // check larva
  write_stage(larva, patchID, 1, *logFiles[1], &Mosquito::print_aquatic, myPRNG, curTime);
  
  // check pupa
  write_stage(pupa, patchID, 2, *logFiles[2], &Mosquito::print_aquatic, myPRNG, curTime);
  
  // check adult males
  write_stage(adult_male, patchID, 3, *logFiles[3], &Mosquito::print_male, myPRNG, curTime);
  
  // check females
  write_stage(adult_female, patchID, 4, *logFiles[4], &Mosquito::print_female, myPRNG, curTime);

}

















