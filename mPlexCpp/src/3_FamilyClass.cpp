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
               const Rcpp::List& eggReleases_,
               const Rcpp::List& matedFemaleReleases_) : Patch::Patch(patchID_,
                                                                      maleReleases_,
                                                                      femaleReleases_,
                                                                      eggReleases_,
                                                                      matedFemaleReleases_)
{
  
  /****************
  * SET POPULATIONS
  ****************/
  fillPopulation(patchID, eggs, larva, pupa,
                 adult_male, adult_female, unmated_female,
                 CreateMosquitoesFamily);

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
void Family::reset_Patch(){
  
  /****************
   * RESET POPULATIONS
   ****************/
  eggs.clear();
  larva.clear();
  pupa.clear();
  adult_male.clear();
  adult_female.clear();
  unmated_female.clear();
  
  /****************
   * SET POPULATIONS
   ****************/
  fillPopulation(patchID, eggs, larva, pupa,
                 adult_male, adult_female, unmated_female,
                 CreateMosquitoesFamily);
  
  /****************
   * RESET RELEASES
   ****************/
  releaseM = releaseM0;
  releaseF = releaseF0;
  releaseE = releaseE0;
  releaseMF = releaseMF0;
};

/******************************************************************************
 * Class Functions
******************************************************************************/
void Family::oneDay_layEggs(prng& myPRNG){
  
  // loop over mated, adult females
  for(auto female : adult_female){
    
    // get number of new offspring
    index = myPRNG.get_rpois(parameters::instance().get_beta()
                                        * reference::instance().get_s(female.get_genotype()));
    
    // create new eggs
    for(size_t it=0; it<index; ++it){
      // protect increment and return of new id
      #pragma omp critical (newIDLock)
      newID = BigBrother::instance().get_ID();
      
      eggs.emplace_back(Mosquito(0, newID));
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
  if(!adult_male.empty()){
    for(auto& mos : adult_male){
      ADM_log << parameters::instance().get_t_now() <<  "," << patchID << "," << mos.print_maleFam();
    }
  } else {
    ADM_log << parameters::instance().get_t_now() <<  "," << patchID << ",,,,\n";
  }

  
  // write adult females
  if(!adult_female.empty()){
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
void CreateMosquitoesFamily(const int& Leq, const int& minAge,
                            const dVec& ageDist, popVec& returnPop){
  
  // loop over each age that mosquitoes can be
  for(size_t age = 0; age < ageDist.size(); age++){
    // loop over number of mosquitoes to create
    for(size_t num=0; num < round(Leq * ageDist[age]); num++){
      
      // add new mosquito to population
      returnPop.emplace_back(Mosquito(minAge + age, BigBrother::instance().get_ID()));
      
    } // end loop over number of mosquitoes
  } // end loop over age distribution
  
}















