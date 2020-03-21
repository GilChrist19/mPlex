///////////////////////////////////////////////////////////////////////////////
//     ________ __ __  _______ 
//    / ____/ //_//  |/  / __ \
//   / /   / ,<  / /|_/ / /_/ /
//  / /___/ /| |/ /  / / _, _/ 
//  \____/_/ |_/_/  /_/_/ |_|  
//   
///////////////////////////////////////////////////////////////////////////////

#include "2_Patch.hpp"
#include "4_bigBrother.hpp"


using dMat = std::vector<dVec>;
using sMat = std::vector<sVec>;

/******************************************************************************
 * Familial Class
******************************************************************************/
#ifndef FAMILY_MPLEX
#define FAMILY_MPLEX

/****************
 * SETUP FUNCTION
 ****************/
// only works for setting up family relations
void CreateMosquitoesFamily(const int& Leq, const int& minAge,
                            const dVec& ageDist, bigBrother& myBB, popVec& returnPop);

// setup proper age distribution in initial population
std::vector<double> markovDist(const double& life, const int& time);
std::vector<double> popDist(const double& mu, const double& alpha, const int& larvalEQ, const iVec& timeAq);

void fillPopulation(const int& patchID, popVec& eggVec, popVec& larvaVec, popVec& pupaVec,
                    popVec& aMaleVec, popVec& aFemaleVec, popVec& unFemaleVec, bigBrother& myBB,
                    void (*populationFill)(const int&, const int&, const dVec&, bigBrother&, popVec&) );

/****************
 * CLASS
 ****************/
class Family : public Patch{
public:
  
  //Constructor and destructor
  Family(const int& patchID_,
        const Rcpp::List& maleReleases_,
        const Rcpp::List& femaleReleases_,
        const Rcpp::List& larvaeReleases_,
        bigBrother& myBB);
  virtual ~Family();
  
  // delete copy constructor and assignment operator
  Family(const Family&) = delete;
  Family& operator=(const Family&&) = delete;
  
  // default move semantics
  Family(Family&&);
  Family& operator=(Family&&);
  
  // functions for this class
  void  oneDay_layEggs(prng& myPRNG, bigBrother& myBB);
  void  reset_Patch(bigBrother& myBB);
  
private:
  
  // used in reproduction
  int index;
  std::string newID;

};

#endif






