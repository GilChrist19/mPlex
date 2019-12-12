///////////////////////////////////////////////////////////////////////////////
//                           ____  __          ______          
//                ____ ___  / __ \/ /__  _  __/ ____/___  ____ 
//               / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
//              / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
//             /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/ 
//                                              /_/   /_/      
///////////////////////////////////////////////////////////////////////////////

#include "2_Patch.hpp"


using dMat = std::vector<dVec>;
using sMat = std::vector<sVec>;


/****************
 * SETUP FUNCTION
 ****************/
std::vector<double> popDist(const double& mu, const double& alpha, const int& larvalEQ, const iVec& timeAq);

/******************************************************************************
 * Daisy Class
******************************************************************************/
#ifndef DAISY_MPLEX
#define DAISY_MPLEX


/****************
 * SETUP FUNCTION
 ****************/
// for setting up Daisy and multiLocus patches
void CreateMosquitoes2Allele(const int& Leq, const int& minAge, const dVec& ageDist,
                             const Rcpp::ListOf<Rcpp::List>& aTypes, popVec& returnPop);

/****************
 * CLASS
 ****************/
class Daisy: public Patch{
public:
  
  
  // constructor
  Daisy(const int& patchID_,
        const Rcpp::List& maleReleases_,
        const Rcpp::List& femaleReleases_,
        const Rcpp::List& larvaeReleases_);
  // destructor
  virtual ~Daisy();
  
  // delete copy constructor and assignment operator
  Daisy(const Daisy&) = delete;
  Daisy& operator=(const Daisy&) = delete;
  
  // default move semantics
  Daisy(Daisy&&);
  Daisy& operator=(Daisy&&);
  
  // functions specific to this class
  void oneDay_layEggs();
  void DaisyOffspring(const std::string& fGen, const std::string& mGen);
  void reset_Patch(); 
  
private:
  // these are values that are used in the mating function
  int numAlleles;
  int index;
  dMat fProbs;
  dMat mProbs;
  sMat fAllele;
  sMat mAllele;
  
  std::string holdAllele;
  std::unordered_map<std::string, double>   duplicates;
  std::unordered_map<std::string, double>::iterator value;
  
  sVec finalGenotypes;
  dVec finalProbs;
  
  sVec holdGens;
  dVec holdProbs;
  
  // used in reproduction
  iVec newEggs;

};

#endif

/******************************************************************************
 * Multiplex multiLocus Class
 ******************************************************************************/
#ifndef MULTILOCUS_MPLEX
#define MULTILOCUS_MPLEX


/****************
 * CLASS
 ****************/
class multiLocus: public Patch{
public:
  
  // constructor
  multiLocus(const int& patchID_,
              const Rcpp::List& maleReleases_,
              const Rcpp::List& femaleReleases_,
              const Rcpp::List& larvaeReleases_);
  // destructor
  virtual ~multiLocus();
  
  // delete copy constructor and assignment operator
  multiLocus(const multiLocus&) = delete;
  multiLocus& operator=(const multiLocus&) = delete;
  
  // default move semantics
  multiLocus(multiLocus&&);
  multiLocus& operator=(multiLocus&&);
  
  // functions specific to this class
  void oneDay_layEggs();
  void MultiplexOffspring_mLoci(const std::string& fGen, const std::string& mGen);
  void reset_Patch();
  
private:
  // these are values that are used in the mating function
  int numAlleles;
  int index;
  dMat fProbs;
  dMat mProbs;
  sMat fAllele;
  sMat mAllele;
  
  std::string holdAllele;
  std::unordered_map<std::string, double>   duplicates;
  std::unordered_map<std::string, double>::iterator value;
    
  sVec finalGenotypes;
  dVec finalProbs;
  
  sVec holdGens;
  dVec holdProbs;
  
  // used in reproduction
  iVec newEggs;
  
};

#endif

/******************************************************************************
 * Multiplex oneLocus Class
******************************************************************************/
#ifndef ONELOCUS_MPLEX
#define ONELOCUS_MPLEX


/****************
 * SETUP FUNCTION
 ****************/
// for setting up oneLocus patch
void CreateMosquitoes2Loci(const int& Leq, const int& minAge, const dVec& ageDist,
                           const Rcpp::ListOf<Rcpp::List>& aTypes, popVec& returnPop);

/****************
 * CLASS
 ****************/
class oneLocus: public Patch{
public:
  
  // constructor
  oneLocus(const int& patchID_,
             const Rcpp::List& maleReleases_,
             const Rcpp::List& femaleReleases_,
             const Rcpp::List& larvaeReleases_);
  // destructor
  virtual ~oneLocus();
  
  // delete copy constructor and assignment operator
  oneLocus(const oneLocus&) = delete;
  oneLocus& operator=(const oneLocus&) = delete;
  
  // default move semantics
  oneLocus(oneLocus&&);
  oneLocus& operator=(oneLocus&&);
  
  // functions specific to this class
  void oneDay_layEggs();
  void MultiplexOffspring_oLocus(const std::string& fGen, const std::string& mGen);
  void reset_Patch();
  
private:
  // these are values that are used in the mating function
  int numLoci;
  int index;
  bool fScore;
  bool mScore;
  
  sVec holdGens1;
  sVec holdGens2;
  sVec holdGens3;
  sVec fAllele;
  sVec mAllele;
  
  dVec holdProbs1;
  dVec holdProbs2;
  dVec holdProbs3;
  dVec fProbs;
  dVec mProbs;
  
  std::string holdAllele;
  std::unordered_map<std::string, double>   duplicates;
  std::unordered_map<std::string, double>::iterator value;
  
  // used in reproduction
  iVec newEggs;
  
};

#endif

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
                            const dVec& ageDist, popVec& returnPop);

/****************
 * CLASS
 ****************/
class Family : public Patch{
public:
  
  //Constructor and destructor
  Family(const int& patchID_,
        const Rcpp::List& maleReleases_,
        const Rcpp::List& femaleReleases_,
        const Rcpp::List& larvaeReleases_);
  virtual ~Family();
  
  // delete copy constructor and assignment operator
  Family(const Family&) = delete;
  Family& operator=(const Family&&) = delete;
  
  // default move semantics
  Family(Family&&);
  Family& operator=(Family&&);
  
  // functions for this class
  void  oneDay_layEggs();
  void  reset_Patch();
  void  init_output(std::ofstream& ADM_log, std::ofstream& ADF_log);
  void  oneDay_writeOutput(std::ofstream& ADM_log, std::ofstream& ADF_log);
  
private:
  
  // used in reproduction
  int index;

};

#endif






