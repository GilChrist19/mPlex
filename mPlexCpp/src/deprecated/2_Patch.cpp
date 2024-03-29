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
// specific to each child


/******************************************************************************
* default (compiler-generated) move semantics
******************************************************************************/
Patch::Patch(Patch&& rhs) = default;
Patch& Patch::operator=(Patch&& rhs) = default;


/******************************************************************************
* Population Dynamics for One Day
******************************************************************************/
/* population dynamics */
void Patch::oneDay_popDynamics(){

  // Death
  oneDay_eggDeath();
  oneDay_larvaDeath();
  oneDay_pupaDeath();
  oneDay_adultDeath();
  
  // Age
  oneDay_eggAge();
  oneDay_larvaeAge();
  oneDay_pupaAge();
  oneDay_adultAge();
  
  // Mature
  oneDay_pupaMaturation();
  oneDay_larvaMaturation();
  oneDay_eggMaturation();
  
  // Mate
  oneDay_mating();
  
  // New Eggs
  oneDay_layEggs();
  
  // Releases
  oneDay_Releases();
  
  // Migration out
  oneDay_migrationOut();

};

/******************************************************************************
* Functions
******************************************************************************/
/**************************************
 * Death
***************************************/
void Patch::oneDay_eggDeath(){

  // set bernoulli
  prng::instance().set_cBern(parameters::instance().get_mu(0));
  
  // Loop over all eggs in the vector
  for(auto it = eggs.rbegin(); it != eggs.rend(); ++it){
    // If it is your time, swap and remove
    if(prng::instance().get_cBern() ){
      std::swap(*it, eggs.back());
      eggs.pop_back();
    }
  } // end loop
  
}

void Patch::oneDay_larvaDeath(){
  
  int stage = parameters::instance().get_stage_time(1);
  double alpha = parameters::instance().get_alpha(patchID);
  alpha = std::pow(alpha/(alpha + larva.size()), 1.0/stage);
  alpha = 1.0-alpha*(1.0-parameters::instance().get_mu(1));
  
  // set bernoulli
  prng::instance().set_cBern(alpha);

  // Loop over all larva in the vector
  for(auto it = larva.rbegin(); it != larva.rend(); ++it){
    // If it is your time, swap and remove
    if(prng::instance().get_cBern() ){
      std::swap(*it, larva.back());
      larva.pop_back();
    }
  } // end loop

}

void Patch::oneDay_pupaDeath(){
  
  // set bernoulli
  prng::instance().set_cBern(parameters::instance().get_mu(2) );

  // Loop over all pupa in the vector
  for(auto it = pupa.rbegin(); it != pupa.rend(); ++it){
    // If it is your time, swap and remove
    if(prng::instance().get_cBern() ){
      std::swap(*it, pupa.back());
      pupa.pop_back();
    }
  } // end loop
  
}

void Patch::oneDay_adultDeath(){
  
  double probs;
  double minusMu = 1.0 - parameters::instance().get_mu(3);
  
  // Loop over all adult males in the vector
  for(auto it = adult_male.rbegin(); it != adult_male.rend(); ++it){
    // get genotype specific chance of death
    probs = minusMu * reference::instance().get_omega(it->get_genotype());
    // If it is your time, swap and remove
    if(prng::instance().get_rBern(1.0-probs) ){
      std::swap(*it, adult_male.back());
      adult_male.pop_back();
    }
  } // end loop
  
  // Loop over all adult females in the vector
  for(auto it = adult_female.rbegin(); it != adult_female.rend(); ++it){
    // get genotype specific chance of death
    probs = minusMu * reference::instance().get_omega(it->get_genotype());
    // If it is your time, swap and remove
    if(prng::instance().get_rBern(1.0-probs) ){
      std::swap(*it, adult_female.back());
      adult_female.pop_back();
    }
  } // end loop
  
}

/**************************************
 * Age
***************************************/
void Patch::oneDay_eggAge(){
  // loop over eggs, age them all one day
  for(auto& it : eggs){
    it.age_one_day();
  }
}

void Patch::oneDay_larvaeAge(){
  // loop over larva, age them all one day
  for(auto& it : larva){
    it.age_one_day();
  }
}

void Patch::oneDay_pupaAge(){
  // loop over pupa, age them all one day
  for(auto& it : pupa){
    it.age_one_day();
  }
}

void Patch::oneDay_adultAge(){
  // loop over males
  for(auto& it : adult_male){
    it.age_one_day();
  }
  // loop over females
  for(auto& it :  adult_female){
    it.age_one_day();
  }
}

/**************************************
 * Mature
***************************************/
void Patch::oneDay_pupaMaturation(){
  
  int age = parameters::instance().get_stage_sum(2);
  
  // Loop over all pupa in the vector
  for(auto it = pupa.rbegin(); it != pupa.rend(); ++it){
    
    // If it is your time, we move in
    if(it->get_age() > age){
      
      // gender mill, mwuahahahaha
      if(prng::instance().get_rBern(reference::instance().get_phi(it->get_genotype())) ){
        // female!
        // check if it actually makes it
        if(prng::instance().get_rBern(reference::instance().get_xiF(it->get_genotype())) ){
          // you become an unmated adult female
          unmated_female.push_back(*it);
        }
        // you actually died
        std::swap(*it, pupa.back());
        pupa.pop_back();
      } else {
        // male!
        // see if you actually make it
        if(prng::instance().get_rBern(reference::instance().get_xiM(it->get_genotype())) ){
          // you become and adult male
          adult_male.push_back(*it);
        }
        // you died
        std::swap(*it, pupa.back());
        pupa.pop_back();
        
      } // end gender mill
      
    } // end age statement
    
  } // end loop
  
}

void Patch::oneDay_larvaMaturation(){
  
  int age = parameters::instance().get_stage_sum(1);
  
  // Loop over all larva in the vector
  for(auto it = larva.rbegin(); it != larva.rend(); ++it){
    // If it is your time, swap and remove
    if(it->get_age() > age){
      // put maturing egg into larva
      pupa.push_back(*it);
      // swap position and remove egg
      std::swap(*it, larva.back());
      larva.pop_back();
    }
  } // end loop
  
}

void Patch::oneDay_eggMaturation(){
  
  int age = parameters::instance().get_stage_time(0);
  
  // Loop over all eggs in the vector
  for(auto it = eggs.rbegin(); it != eggs.rend(); ++it){
    // If it is your time, swap and remove
    if(it->get_age() > age){
      // put maturing egg into larva
      larva.push_back(*it);
      // swap position and remove egg
      std::swap(*it, eggs.back());
      eggs.pop_back();
    }
  } // end loop

}

/**************************************
 * Mate
***************************************/
void Patch::oneDay_mating(){
  
  if((adult_male.size() != 0) && (unmated_female.size() !=0) ){
    
    // clear holder objects for reuse
    genNames.clear();
    genProbs.clear();

    for(auto& male : adult_male){
      genNames.push_back(male.get_genotype());
      genProbs.push_back(reference::instance().get_eta(male.get_genotype()));
    }

    // set sample
    prng::instance().set_oneSample(genProbs);

    // loop over unmated females
    for(auto it= unmated_female.rbegin(); it != unmated_female.rend(); ++it){
      //get and set mate
      mateName = genNames[prng::instance().get_oneSample()];
      it->set_mate(mateName);
      // move to real adult females and remove
      adult_female.push_back(*it);
      unmated_female.pop_back();
    } // end loop over unmated females

  } // only do if there are unmated females and males to mate with
  
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
                                       releaseM.back().pop_names[it])
                              );
    } // end loop over released mosquitoes
    
    // remove release
    releaseM.pop_back();
  } // end male release
  
  /****************
   * FEMALE
  ****************/
  if( (!releaseF.empty()) && (releaseF.back().release_time <= parameters::instance().get_t_now()) ){
    
    
    for(size_t it = 0; it < releaseF.back().pop_ages.size(); ++it){
      
      adult_female.push_back(Mosquito(releaseF.back().pop_ages[it],
                                         releaseF.back().pop_names[it])
                                );
    } // end loop over releases
    
    // remove release
    releaseF.pop_back();
  } // end female release
  
  /****************
   * LARVA
  ****************/
  if( (!releaseE.empty()) && (releaseE.back().release_time <= parameters::instance().get_t_now()) ){
    
    
    for(size_t it = 0; it < releaseE.back().pop_ages.size(); ++it){
      
      eggs.push_back(Mosquito(releaseE.back().pop_ages[it],
                                         releaseE.back().pop_names[it]));
    } // end loop over releases
    
    // remove release
    releaseE.pop_back();
  } // end larva release
  
}

/**************************************
 * Migration
***************************************/
void Patch::oneDay_migrationOut(){
  
  // clear old migration
  for(size_t i=0; i < parameters::instance().get_n_patch(); ++i){
    maleMigration[i].clear();
    femaleMigration[i].clear();
  }

  
  // used a lot
  int patch;
  dVec Probs = prng::instance().get_rdirichlet(parameters::instance().get_male_migration(patchID)); 
  
  /****************
   * MALE
  ****************/
  // set oneSample
  prng::instance().set_oneSample(Probs);
  
  // loop over all adult males
  for(auto it = adult_male.rbegin(); it != adult_male.rend(); ++it){
    //get which patch he goes to
    patch = prng::instance().get_oneSample();
    
    // if not this patch
    if(patch != patchID){ 
      // add to new patch it goes to
      maleMigration[patch].push_back(*it);
      // remove
      std::swap(*it, adult_male.back());
      adult_male.pop_back();
    } // end if this patch
    
  } // end loop over adult males
  
  
  /****************
   * FEMALE
  ****************/
  // get slightly more variance in probability
  Probs = prng::instance().get_rdirichlet(parameters::instance().get_female_migration(patchID));
  
  // reset oneSample for female probs
  prng::instance().set_oneSample(Probs);
  
  // loop over all females
  for(auto it = adult_female.rbegin(); it != adult_female.rend(); ++it){
    // get which patch she goes to
    patch = prng::instance().get_oneSample();
    
    // if not this patch
    if(patch != patchID){
      // add to new place
      femaleMigration[patch].push_back(*it);
      // remove
      std::swap(*it, adult_female.back());
      adult_female.pop_back();
    } // end if
    
  } // end loop over females
  
  
  /****************
   * BATCH
  ****************/
  if(prng::instance().get_rBern(parameters::instance().get_batchProbs(patchID)) ){
    
    // which patch will they migrate to
    // set oneSample, then get it
    prng::instance().set_oneSample(parameters::instance().get_batchLocation(patchID));
    patch = prng::instance().get_oneSample();
    
    // set binomial for inside loop
    prng::instance().set_cBern(parameters::instance().get_batchMale(patchID));

    // do male migration
    for(auto mos = adult_male.rbegin(); mos != adult_male.rend(); ++mos){
      // see if it moves
      if(prng::instance().get_cBern() ){
        // put him in new patch migration set
        maleMigration[patch].push_back(*mos);
        // remove him from here
        std::swap(*mos, adult_male.back());
        adult_male.pop_back();
      }
    } // end males
    
    // reset binom for females
    prng::instance().set_cBern(parameters::instance().get_batchFemale(patchID));

    // female batch migration
    for(auto mos = adult_female.rbegin(); mos != adult_female.rend(); ++mos){
      // see if she moves
      if(prng::instance().get_cBern() ){
        // put her in the new patch
        femaleMigration[patch].push_back(*mos);
        // remove from current patch
        std::swap(*mos, adult_female.back());
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
  
  // write males
  for(auto& mos : adult_male){
    ADM_log << parameters::instance().get_t_now() <<  "," << patchID << "," << mos.print_male();
  }
  
  // write unmated females
  for(auto& mos : unmated_female){
    ADF_log << parameters::instance().get_t_now() <<  "," << patchID << "," << mos.print_female();
  }
  
}


void Patch::oneDay_writeOutput(std::ofstream& ADM_log, std::ofstream& ADF_log){
  
  // write males
  if(adult_male.size() > 0){
    for(auto& mos : adult_male){
      ADM_log << parameters::instance().get_t_now() <<  "," << patchID << "," << mos.print_male();
    }
  } else {
    ADM_log << parameters::instance().get_t_now() <<  "," << patchID << ",,\n";
  }

  
  // write adult females
  if(adult_female.size() > 0){
    for(auto& mos : adult_female){
      ADF_log << parameters::instance().get_t_now() <<  "," << patchID << "," << mos.print_female();
    }
  } else {
    ADF_log << parameters::instance().get_t_now() <<  "," << patchID << ",,,\n";
  }

  
}

















