///////////////////////////////////////////////////////////////////////////////
//     ________ __ __  _______ 
//    / ____/ //_//  |/  / __ \
//   / /   / ,<  / /|_/ / /_/ /
//  / /___/ /| |/ /  / / _, _/ 
//  \____/_/ |_/_/  /_/_/ |_|  
//   
///////////////////////////////////////////////////////////////////////////////

#include "3_DriveDefinitions.hpp"


/******************************************************************************
 * Constructor & Destructor
******************************************************************************/
Family::Family(const int& patchID_,
              const Rcpp::List& maleReleases_,
              const Rcpp::List& femaleReleases_,
              const Rcpp::List& eggReleases_,
              bigBrother& myBB) : Patch::Patch(patchID_, maleReleases_,
                                               femaleReleases_, eggReleases_)
{
  
  /****************
  * SET POPULATIONS
  ****************/
  fillPopulation(patchID, eggs, larva, pupa,
                 adult_male, adult_female, unmated_female,
                 myBB, CreateMosquitoesFamily);

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
void Family::reset_Patch(bigBrother& myBB){
  
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
                 myBB, CreateMosquitoesFamily);  
  
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
void Family::oneDay_layEggs(prng& myPRNG, bigBrother& myBB){
  
  // loop over mated, adult females
  for(auto female : adult_female){
    
    // get number of new offspring
    //  check if deterministic or if it follows a poisson distribution
    if(parameters::instance().get_beta_const()){
      // constant beta, no distribution
      index = parameters::instance().get_beta() * reference::instance().get_s(female.get_myID());
    } else {
      // poisson distributed beta
      index = myPRNG.get_rpois(parameters::instance().get_beta()
                               * reference::instance().get_s(female.get_myID()));
    }
    
    // create new eggs
    for(size_t it=0; it<index; ++it){
      // // protect increment and return of new id
      // #pragma omp critical (newIDLock)
      // newID = BigBrother::instance().get_ID();
      
      eggs.emplace_back(Mosquito(0, myBB.get_ID(), female.get_myID(), female.get_mate()));
    } // end making eggs
  } // end loop over females
  
}



/******************************************************************************
 * FUNCTIONS
 ******************************************************************************/

/**************************************
 * SETUP
 **************************************/
std::vector<double> markovDist(const double& life, const int& time){
  
  // setup matrix to solve
  arma::Mat<double> markovMat(time, time, arma::fill::eye);
  
  // create and fill offDiagonal vector
  arma::Col<double> offDiag(time-1);
  offDiag.fill((1.0 - life) - 1.0);
  
  // put off diagonal in matrix
  markovMat.diag(1) = offDiag;
  
  // invert matrix
  markovMat = markovMat.i();
  
  // get normalized vector of larval ratios
  arma::Row<double> solVec(markovMat.row(0)/arma::sum(markovMat.row(0)));
  
  // store as standard vector, both for return and because arma::Row doesn't 
  //  have some of the functions I need
  std::vector<double> retVec(solVec.begin(), solVec.end());
  
  return(retVec);
}

std::vector<double> popDist(const double& mu, const double& alpha, const int& larvalEQ, const iVec& timeAq) {
  // This comes from MGDrivE, where Sean solved the initial conditions as a CTMC
  
  /***********************************/
  // setup different death rates
  /***********************************/
  double epLife(1.0 - mu);
  double lLife(std::pow(alpha/(alpha + larvalEQ), 1.0/timeAq[1]) * epLife );
  
  /***********************************/
  // Solve larval distribution
  /***********************************/
  std::vector<double> hold(markovDist(lLife, timeAq[1]));
  
  /***********************************/
  // Egg distribution at the front
  /***********************************/
  for(int i = timeAq[0]; i > 0; i--){
    hold.insert(hold.begin(), hold.front()/epLife);
  }
  
  /***********************************/
  // Pupal distribution at the end
  /***********************************/
  hold.push_back(hold.back() * lLife);
  for(int i = timeAq[2]-1; i > 0; i--){
    hold.push_back(hold.back() * epLife);
  }
  
  // return normalized vector of initial conditions
  return(hold);
}


void fillPopulation(const int& patchID, popVec& eggVec, popVec& larvaVec, popVec& pupaVec,
                    popVec& aMaleVec, popVec& aFemaleVec, popVec& unFemaleVec, bigBrother& myBB,
                    void (*populationFill)(const int&, const int&, const dVec&, bigBrother&, popVec&)
){
  
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
  
  
  /*********/
  // eggs
  /*********/
  minAge = 0;
  maxAge = parameters::instance().get_stage_time(0) - 1;
  
  // set age distribution vector
  distHold.resize(maxAge+1);
  std::copy(aDist.begin(), aDist.begin() + maxAge+1, distHold.begin());
  
  // reserve estimated population size
  eggVec.reserve(round(1.1*parameters::instance().get_larva_eq(patchID) * accumulate(distHold.begin(), distHold.end(), 0.0)));
  
  // fill initial population
  populationFill(parameters::instance().get_larva_eq(patchID),
                 minAge, distHold, myBB, eggVec);
  
  
  /*********/
  // larva
  /*********/
  minAge = parameters::instance().get_stage_time(0);
  maxAge = parameters::instance().get_stage_sum(1)-1;
  
  // set age distribution vector
  distHold.resize(maxAge+1 - minAge);
  std::copy(aDist.begin() + minAge, aDist.begin() + maxAge + 1, distHold.begin());
  
  // reserve estimated population size
  larvaVec.reserve(round(1.1*parameters::instance().get_larva_eq(patchID) * accumulate(distHold.begin(), distHold.end(), 0.0)));
  
  // fill initial population
  populationFill(parameters::instance().get_larva_eq(patchID),
                 minAge, distHold, myBB, larvaVec);
  
  
  /*********/
  // pupa
  /*********/
  minAge = parameters::instance().get_stage_sum(1);
  maxAge = parameters::instance().get_stage_sum(2) - 1;
  
  // set age distribution vector
  distHold.resize(maxAge+1 - minAge);
  std::copy(aDist.begin() + minAge, aDist.begin() + maxAge + 1, distHold.begin());
  
  // reserve estimated population size
  pupaVec.reserve(round(1.1*parameters::instance().get_larva_eq(patchID) * accumulate(distHold.begin(), distHold.end(), 0.0)));
  
  // fill initial population
  populationFill(parameters::instance().get_larva_eq(patchID),
                 minAge, distHold, myBB, pupaVec);
  
  
  /***********************************/
  // Solve adult distribution
  /***********************************/
  minAge = parameters::instance().get_stage_sum(2);
  
  /*********/
  // male
  /*********/
  // set age distribution vector
  maxAge = round(std::min(parameters::instance().get_male_max_age() - minAge,
                          parameters::instance().get_stage_time(3) * 3));
  aDist = markovDist(1.0 - parameters::instance().get_mu(3), maxAge);
  
  // reserve estimated population size
  aMaleVec.reserve(round(0.75*parameters::instance().get_adult_pop_eq(patchID)));
  
  // fill initial population
  populationFill(parameters::instance().get_adult_pop_eq(patchID)/2,
                 minAge, aDist, myBB, aMaleVec);
  
  /*********/
  // female
  /*********/
  // set age distribution vector
  maxAge = round(std::min(parameters::instance().get_female_max_age() - minAge,
                          parameters::instance().get_stage_time(3) * 3));
  aDist = markovDist(1.0 - parameters::instance().get_mu(3), maxAge);
  
  // reserve estimated population size
  int popSize(round(0.75*parameters::instance().get_adult_pop_eq(patchID)));
  aFemaleVec.reserve(popSize);
  unFemaleVec.reserve(popSize);
  
  // fill initial population
  populationFill(parameters::instance().get_adult_pop_eq(patchID)/2,
                 minAge, aDist, myBB, unFemaleVec);
  
}


void CreateMosquitoesFamily(const int& Leq, const int& minAge,
                            const dVec& ageDist, bigBrother& myBB,
                            popVec& returnPop){
  
  // loop over each age that mosquitoes can be
  for(size_t age = 0; age < ageDist.size(); age++){
    // loop over number of mosquitoes to create
    for(size_t num=0; num < round(Leq * ageDist[age]); num++){
      
      // add new mosquito to population
      returnPop.emplace_back(Mosquito(minAge + age, myBB.get_ID(), "0", "0"));
      
    } // end loop over number of mosquitoes
  } // end loop over age distribution
  
}















