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
    std::sort(releaseM.begin(), releaseM.end(), [](release_int a, release_int b){
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
    std::sort(releaseF.begin(), releaseF.end(), [](release_int a, release_int b){
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
    std::sort(releaseL.begin(), releaseL.end(), [](release_int a, release_int b){
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

























/**************************************
* Mating
***************************************/
void Patch::oneGen_Mating(Patch* E){

  // if there are males, do it, else, return zeros
  if(arma::sum(E->Population.slice(E->tNow).col(0)) > 0){

    // two things I need
    size_t genotypesN = cube_base::instance().get_genotypesN();

    // probabilities, need conversion to doubles, then normalize
    Column_dbl probs = arma::conv_to<Column_dbl>::from(E->Population.slice(E->tNow).col(0)) % cube_base::instance().get_eta();
    probs = probs/ arma::sum(probs);


    // loop over females
    for(size_t i=0; i<genotypesN; ++i){
      // if there are females, do it, else zeros
      if(E->Population.slice(E->tNow).at(i,1) > 0){
        R::rmultinom(E->Population.slice(E->tNow).at(i,0), probs.begin(), genotypesN, E->holdVec.begin());
        E->matedFemales.row(i) = arma::conv_to<arma::Row<int> >::from(E->holdVec);
      } else {
        E->matedFemales.row(i).zeros();
      }
    } // end for loop
  } else {
    E->matedFemales.zeros();
  }

}

/**************************************
* Lay Eggs
***************************************/
void Patch::oneGen_layEggs(Patch* E){
  // Use several times
  size_t genotypesN = cube_base::instance().get_genotypesN();

  // setup hold matrix
  double offspring;

  // beta and S combination
  Column_dbl betaS = E->beta * cube_base::instance().get_s();

  // loop over all offspring genotypes
  for(size_t slice=0; slice<genotypesN; ++slice){
    /* hadamard product of AF1dly with cube_offspring */
    E->AF1_hadamard = arma::conv_to<Matrix_dbl>::from(E->matedFemales) % cube_base::instance().get_cube_z(slice) % cube_base::instance().get_tau_z(slice);
    /* multiply element-wise by beta * s */
    E->AF1_hadamard.each_col() %= betaS;
    offspring = arma::accu(E->AF1_hadamard);
    /* sum over genotype and put into eggs matrix */
    if(offspring>0){
      E->Population.slice(E->tNow).at(slice,2) = R::rpois(offspring);
    }
  } // end for loop
}

/**************************************
* Male/Female
***************************************/
void Patch::oneGen_sexDetermination(Patch* E){

  // get sex ratios
  Column_dbl phi = cube_base::instance().get_phi();
  size_t genotypesN= cube_base::instance().get_genotypesN();

  // pull binomial by sex ratio, save as male/female
  for(size_t i=0; i<genotypesN; ++i){
    if(E->Population.slice(E->tNow).at(i,2)>0){
      E->Population.slice(E->tNow).at(i,3) = R::rbinom(E->Population.slice(E->tNow).at(i,2),phi[i]);
      E->Population.slice(E->tNow).at(i,4) = E->Population.slice(E->tNow).at(i,2) - E->Population.slice(E->tNow).at(i,3);
    }
  }// end loop
}

/**************************************
* Sample for Tomorrow
***************************************/
void Patch::oneGen_adultSampling(Patch* E){

  //These hypergeometric things are a problem!!
  /* The hypergeometric takes integer bins as input.
   * however, multiplying the integers by doubls should give doubles, which I
   * am casting back as ints. This is creating a truncation that at low densities
   * could cause problems. Is th is ok? It matches the deterministic because
   * there is a cutoff in that one, where anything less than 1 is turned to 0.
   */

  // use several times
  int genotypesN = cube_base::instance().get_genotypesN();
  int drawSize;

  // males for tomorrow
  drawSize = R::rbinom(E->AdPopEQ, 0.5);
  E->holdVec = arma::conv_to<Column_int>::from(arma::conv_to<Column_dbl>::from(E->Population.slice(E->tNow).col(3))
                                                 % cube_base::instance().get_xiM());
  E->Population.slice(E->tNow+1).col(0) = rmvHyperGeo(E->holdVec, drawSize);

  // females for tomorrow
  drawSize = R::rbinom(E->AdPopEQ, 0.5);
  E->holdVec = arma::conv_to<Column_int>::from(arma::conv_to<Column_dbl>::from(E->Population.slice(E->tNow).col(4))
                                                 % cube_base::instance().get_xiF());
  E->Population.slice(E->tNow+1).col(1) = rmvHyperGeo(E->holdVec, drawSize);

}

/**************************************
 * Release Functions
 ***************************************/
void Patch::oneGen_maleReleases(Patch* E){

  //These hypergeometric things are a problem!!
  /* The hypergeometric takes integer bins as input.
   * however, multiplying the integers by doubls should give doubles, which I
   * am casting back as ints. This is creating a truncation that at low densities
   * could cause problems. Is th is ok? It matches the deterministic because
   * there is a cutoff in that one, where anything less than 1 is turned to 0.
   */

  // releases queued up
  if(!E->releaseM.empty()){
    if(E->releaseM.back().release_time <= E->tNow){

      // adult mosquitoes replaced by release
      int drawSize = R::rbinom(2*arma::sum(E->releaseM.back().release_pop), 0.5);
      E->holdVec = arma::conv_to<Column_int>::from(arma::conv_to<Column_dbl>::from(E->Population.slice(E->tNow).col(0))
                                                     % cube_base::instance().get_xiM());
      E->holdVec = rmvHyperGeo(E->holdVec, drawSize);


      // add new and remove replaces
      E->Population.slice(E->tNow).col(0) += E->releaseM.back().release_pop;
      E->Population.slice(E->tNow).col(0) -= E->holdVec;

      // remove release from list
      E->releaseM.pop_back();
    }
  }
}

void Patch::oneGen_femaleReleases(Patch* E){

  // releases queued up
  if(!E->releaseF.empty()){
    if(E->releaseF.back().release_time <= E->tNow){

      // adult mosquitoes replaced by release
      int drawSize = R::rbinom(2*arma::sum(E->releaseF.back().release_pop), 0.5);
      E->holdVec = arma::conv_to<Column_int>::from(arma::conv_to<Column_dbl>::from(E->Population.slice(E->tNow).col(1))
                                                     % cube_base::instance().get_xiF());
      E->holdVec = rmvHyperGeo(E->holdVec, drawSize);

      // add new and remove old
      E->Population.slice(E->tNow).col(1) += E->releaseF.back().release_pop;
      E->Population.slice(E->tNow).col(1) -= E->holdVec;

      // remove release from list
      E->releaseF.pop_back();
    }
  }
}

/**************************************
 * Reset Function
***************************************/
void Patch::reset_Experiment(){

  // reset time
  this->tNow = 0;

  // rezero population cube, set initial adults
  this->Population.zeros();
  this->Population.slice(0).col(0) = ADMt0;
  this->Population.slice(0).col(1) = AF1t0;

  // reset releases
  this->releaseF = releaseF0;
  this->releaseM = releaseM0;
}














