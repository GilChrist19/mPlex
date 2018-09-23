///////////////////////////////////////////////////////////////////////////////
//                           ____  __          ______          
//                ____ ___  / __ \/ /__  _  __/ ____/___  ____ 
//               / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
//              / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
//             /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/ 
//                                              /_/   /_/      
///////////////////////////////////////////////////////////////////////////////

#include "3_DriveDefinitions.hpp"
#include "4_BigBrother.hpp"


/******************************************************************************
 * Constructor & Destructor
******************************************************************************/
Family::Family(const int& patchID_,
              const Rcpp::List& maleReleases_,
              const Rcpp::List& femaleReleases_,
              const Rcpp::List& eggReleases_) : Patch::Patch(patchID_,
                                                              maleReleases_,
                                                              femaleReleases_,
                                                              eggReleases_)
{
  
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
  ageDist.assign(parameters::instance().get_stage_time(0),1);
  CreateMosquitoesFamily(parameters::instance().get_larva_eq(patchID),
                         minAge, ageDist, eggs);

  
  // larva
  larva.reserve(2*parameters::instance().get_larva_eq(patchID));
  minAge = parameters::instance().get_stage_time(0)+1;

  ageDist.clear();
  int counter(0);
  for(int power = minAge; power <= parameters::instance().get_stage_sum(1); ++power, counter+=2){
   ageDist.push_back(std::pow(1.0-parameters::instance().get_mu(1), counter));
  }
  
  CreateMosquitoesFamily(parameters::instance().get_larva_eq(patchID),
                         minAge, ageDist, larva);

  // pupa
  pupa.reserve(2*parameters::instance().get_adult_pop_eq(patchID));

  // adults
  adult_male.reserve(2*parameters::instance().get_adult_pop_eq(patchID));
  adult_female.reserve(2*parameters::instance().get_adult_pop_eq(patchID));
  unmated_female.reserve(2*parameters::instance().get_adult_pop_eq(patchID));
  
  minAge = parameters::instance().get_stage_sum(2)+1;
  ageDist.assign(parameters::instance().get_stage_sum(3) - minAge +1,1);
  
  CreateMosquitoesFamily(parameters::instance().get_adult_pop_eq(patchID)/2,
                         minAge, ageDist, adult_male);
  CreateMosquitoesFamily(parameters::instance().get_adult_pop_eq(patchID)/2,
                         minAge, ageDist, unmated_female);

};

Family::~Family(){};


/******************************************************************************
 * Default (compiler-generated) move semantics
******************************************************************************/
Family::Family(Family&& d) = default;
Family& Family::operator=(Family&& d) = default;

/******************************************************************************
 * Reset
******************************************************************************/
void Family::reset_Patch(const Rcpp::ListOf<Rcpp::List>& aTypes){
  
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
  ageDist.assign(parameters::instance().get_stage_time(0),1);
  CreateMosquitoesFamily(parameters::instance().get_larva_eq(patchID),
                          minAge, ageDist, eggs);
  
  // larva
  larva.clear();
  minAge = parameters::instance().get_stage_time(0)+1;
  
  ageDist.clear();
  int counter(0);
  for(int power = minAge; power <= parameters::instance().get_stage_sum(1); ++power, counter+=2){
    ageDist.push_back(std::pow(1.0-parameters::instance().get_mu(1), counter));
  }
  
  CreateMosquitoesFamily(parameters::instance().get_larva_eq(patchID),
                          minAge, ageDist, larva);
  
  // pupa
  pupa.clear();
  
  // adults
  adult_male.clear();
  adult_female.clear();
  unmated_female.clear();
  
  minAge = parameters::instance().get_stage_sum(2)+1;
  ageDist.assign(parameters::instance().get_stage_sum(3) - minAge +1,1);
  CreateMosquitoesFamily(parameters::instance().get_adult_pop_eq(patchID)/2,
                          minAge, ageDist, adult_male);
  CreateMosquitoesFamily(parameters::instance().get_adult_pop_eq(patchID)/2,
                          minAge, ageDist, unmated_female);
  
  
  /****************
   * RESET RELEASES
   ****************/
  releaseM = releaseM0;
  releaseF = releaseF0;
  releaseE = releaseE0;
};

/******************************************************************************
 * Class Functions
******************************************************************************/
void Family::oneDay_layEggs(){
  
  // loop over mated, adult females
  for(auto female : adult_female){
    
    // get number of new offspring
    index = prng::instance().get_rpois(parameters::instance().get_beta()
                                        * reference::instance().get_s(female.get_genotype()));
    
    // create new eggs
    for(size_t it=0; it<index; ++it){
      eggs.emplace_back(Mosquito(0, BigBrother::instance().get_ID()));
      eggs.back().set_parents(female.get_genotype(), female.get_mate());
    } // end making eggs
  } // end loop over females
  
}

void Family::init_output(std::ofstream& ADM_log, std::ofstream& ADF_log){
  
  // write headers
  if(patchID == 0){
    ADM_log << "Time,Patch,Age,myID,momID,dadID\n";
    ADF_log << "Time,Patch,Age,myID,momID,dadID,Mate\n";
  }
  
  // write males
  for(auto& mos : adult_male){
    ADM_log << parameters::instance().get_t_now() <<  "," << patchID << "," << mos.print_maleFam();
  }
  
  // write unmated females
  for(auto& mos : unmated_female){
    ADF_log << parameters::instance().get_t_now() <<  "," << patchID << "," << mos.print_femaleFam();
  }
  
}

void  Family::oneDay_writeOutput(std::ofstream& ADM_log, std::ofstream& ADF_log){
  
  // write males
  if(adult_male.size() > 0){
    for(auto& mos : adult_male){
      ADM_log << parameters::instance().get_t_now() <<  "," << patchID << "," << mos.print_maleFam();
    }
  } else {
    ADM_log << parameters::instance().get_t_now() <<  "," << patchID << ",,,,\n";
  }

  
  // write adult females
  if(adult_female.size() > 0){
    for(auto& mos : adult_female){
      ADF_log << parameters::instance().get_t_now() <<  "," << patchID << "," << mos.print_femaleFam();
    }
  } else {
    ADF_log << parameters::instance().get_t_now() <<  "," << patchID << ",,,,,\n";
  }
  
}

/******************************************************************************
 * FUNCTIONS
 ******************************************************************************/

/**************************************
 * SETUP
 **************************************/
void CreateMosquitoesFamily(const int& numMos, const int& minAge,
                            const dVec& ageDist, popVec& returnPop){
  
  // holders
  int age;
  
  // set age dist
  prng::instance().set_oneSample(ageDist);
  
  // loop over number of mosquioes to create
  for(size_t count=0; count < numMos; ++count){
    
    // get age
    age = minAge + prng::instance().get_oneSample();
    
    // add new mosquito
    // don't get mom/pop
    returnPop.emplace_back(Mosquito(age, BigBrother::instance().get_ID()));
    
  } // end loop over mosquitoes
  
}















