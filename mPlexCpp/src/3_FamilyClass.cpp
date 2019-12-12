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
  int minAge, maxAge;
  dVec distHold;
  
  // solve aquatic distribution
  dVec aDist = popDist(parameters::instance().get_mu(0),
                       parameters::instance().get_alpha(patchID),
                       parameters::instance().get_larva_eq(patchID),
                       {parameters::instance().get_stage_time(0),
                        parameters::instance().get_stage_time(1),
                        parameters::instance().get_stage_time(2)});
  
  // eggs
  minAge = 0;
  maxAge = parameters::instance().get_stage_time(0) - 1;
  
  distHold.resize(maxAge+1);
  std::copy(aDist.begin(), aDist.begin() + maxAge+1, distHold.begin());
  
  eggs.reserve(round(1.1*parameters::instance().get_larva_eq(patchID) * accumulate(distHold.begin(), distHold.end(), 0.0)));
  
  CreateMosquitoesFamily(parameters::instance().get_larva_eq(patchID),
                          minAge, distHold, eggs);
  
  
  // larva
  minAge = parameters::instance().get_stage_time(0);
  maxAge = parameters::instance().get_stage_sum(1)-1;
  
  distHold.resize(maxAge+1 - minAge);
  std::copy(aDist.begin() + minAge, aDist.begin() + maxAge + 1, distHold.begin());
  
  larva.reserve(round(1.1*parameters::instance().get_larva_eq(patchID) * accumulate(distHold.begin(), distHold.end(), 0.0)));
  
  CreateMosquitoesFamily(parameters::instance().get_larva_eq(patchID),
                          minAge, distHold, larva);
  
  
  // pupa
  minAge = parameters::instance().get_stage_sum(1);
  maxAge = parameters::instance().get_stage_sum(2) - 1;
  
  distHold.resize(maxAge+1 - minAge);
  std::copy(aDist.begin() + minAge, aDist.begin() + maxAge + 1, distHold.begin());
  
  pupa.reserve(round(1.1*parameters::instance().get_larva_eq(patchID) * accumulate(distHold.begin(), distHold.end(), 0.0)));
  
  CreateMosquitoesFamily(parameters::instance().get_larva_eq(patchID),
                          minAge, distHold, pupa);
  
  /***********************************/
  // Solve adult distribution
  /***********************************/
  // basically stolen from popDist() and adapted for adults
  // setup matrix to solve
  arma::Mat<double> markovMat(parameters::instance().get_stage_time(3) * 2,
                              parameters::instance().get_stage_time(3) * 2,
                              arma::fill::eye);
  
  // create and fill offDiagonal vector
  arma::Col<double> offDiag(parameters::instance().get_stage_time(3) * 2-1);
  offDiag.fill(parameters::instance().get_mu(3) - 1.0);
  
  // put off diagonal in matrix
  markovMat.diag(1) = offDiag;
  
  // invert matrix
  markovMat = markovMat.i();
  
  // get normalized vector of larval ratios
  arma::Row<double> solVec(markovMat.row(0)/arma::sum(markovMat.row(0)));
  
  // store as standard vector, both for return and because arma::Row doesn't 
  //  have some of the functions I need
  std::vector<double> hold(solVec.begin(), solVec.end());
  
  
  int popSize(round(1.1 * parameters::instance().get_adult_pop_eq(patchID) * accumulate(hold.begin(), hold.end(), 0.0)));
  
  adult_male.reserve(popSize);
  adult_female.reserve(popSize);
  unmated_female.reserve(popSize);
  
  minAge = parameters::instance().get_stage_sum(2);
  
  CreateMosquitoesFamily(parameters::instance().get_adult_pop_eq(patchID)/2,
                          minAge, hold, adult_male);
  CreateMosquitoesFamily(parameters::instance().get_adult_pop_eq(patchID)/2,
                          minAge, hold, unmated_female);

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
  // holder objects, these are for distribution function
  int minAge, maxAge;
  dVec distHold;
  
  // solve aquatic distribution
  dVec aDist = popDist(parameters::instance().get_mu(0),
                       parameters::instance().get_alpha(patchID),
                       parameters::instance().get_larva_eq(patchID),
                       {parameters::instance().get_stage_time(0),
                        parameters::instance().get_stage_time(1),
                        parameters::instance().get_stage_time(2)});
  
  // eggs
  minAge = 0;
  maxAge = parameters::instance().get_stage_time(0) - 1;
  
  distHold.resize(maxAge+1);
  std::copy(aDist.begin(), aDist.begin() + maxAge+1, distHold.begin());
  
  eggs.reserve(round(1.1*parameters::instance().get_larva_eq(patchID) * accumulate(distHold.begin(), distHold.end(), 0.0)));
  
  CreateMosquitoesFamily(parameters::instance().get_larva_eq(patchID),
                          minAge, distHold, eggs);
  
  
  // larva
  minAge = parameters::instance().get_stage_time(0);
  maxAge = parameters::instance().get_stage_sum(1)-1;
  
  distHold.resize(maxAge+1 - minAge);
  std::copy(aDist.begin() + minAge, aDist.begin() + maxAge + 1, distHold.begin());
  
  larva.reserve(round(1.1*parameters::instance().get_larva_eq(patchID) * accumulate(distHold.begin(), distHold.end(), 0.0)));
  
  CreateMosquitoesFamily(parameters::instance().get_larva_eq(patchID),
                          minAge, distHold, larva);
  
  
  // pupa
  minAge = parameters::instance().get_stage_sum(1);
  maxAge = parameters::instance().get_stage_sum(2) - 1;
  
  distHold.resize(maxAge+1 - minAge);
  std::copy(aDist.begin() + minAge, aDist.begin() + maxAge + 1, distHold.begin());
  
  pupa.reserve(round(1.1*parameters::instance().get_larva_eq(patchID) * accumulate(distHold.begin(), distHold.end(), 0.0)));
  
  CreateMosquitoesFamily(parameters::instance().get_larva_eq(patchID),
                          minAge, distHold, pupa);
  
  /***********************************/
  // Solve adult distribution
  /***********************************/
  // basically stolen from popDist() and adapted for adults
  // setup matrix to solve
  arma::Mat<double> markovMat(parameters::instance().get_stage_time(3) * 2,
                              parameters::instance().get_stage_time(3) * 2,
                              arma::fill::eye);
  
  // create and fill offDiagonal vector
  arma::Col<double> offDiag(parameters::instance().get_stage_time(3) * 2-1);
  offDiag.fill(parameters::instance().get_mu(3) - 1.0);
  
  // put off diagonal in matrix
  markovMat.diag(1) = offDiag;
  
  // invert matrix
  markovMat = markovMat.i();
  
  // get normalized vector of larval ratios
  arma::Row<double> solVec(markovMat.row(0)/arma::sum(markovMat.row(0)));
  
  // store as standard vector, both for return and because arma::Row doesn't 
  //  have some of the functions I need
  std::vector<double> hold(solVec.begin(), solVec.end());
  
  
  int popSize(round(1.1 * parameters::instance().get_adult_pop_eq(patchID) * accumulate(hold.begin(), hold.end(), 0.0)));
  
  adult_male.reserve(popSize);
  adult_female.reserve(popSize);
  unmated_female.reserve(popSize);
  
  minAge = parameters::instance().get_stage_sum(2);
  
  CreateMosquitoesFamily(parameters::instance().get_adult_pop_eq(patchID)/2,
                          minAge, hold, adult_male);
  CreateMosquitoesFamily(parameters::instance().get_adult_pop_eq(patchID)/2,
                          minAge, hold, unmated_female);
  
  
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















