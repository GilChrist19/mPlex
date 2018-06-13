//      __  _____________  ____  ______
//     /  |/  / ____/ __ \/ __ \/ ____/___  ____
//    / /|_/ / / __/ / / / / / / /   / __ \/ __ \
//   / /  / / /_/ / /_/ / /_/ / /___/ /_/ / /_/ /
//  /_/  /_/\____/_____/_____/\____/ .___/ .___/
//                                /_/   /_/


#include "2_Patch.hpp"


using dMat = std::vector<dVec>;
using sMat = std::vector<sVec>;




/******************************************************************************
 * Daisy Class
******************************************************************************/
#ifndef DAISY_MPLEX
#define DAISY_MPLEX

class Daisy: public Patch{
public:
  
  // constructor
  Daisy(const int& patchID_,
        const Rcpp::ListOf<Rcpp::List>& aTypes,
        const Rcpp::ListOf<Rcpp::List>& maleReleases_,
        const Rcpp::ListOf<Rcpp::List>& femaleReleases_,
        const Rcpp::ListOf<Rcpp::List>& larvaeReleases_);
  // destructor
  ~Daisy();
  
  // delete copy constructor and assignment operator
  Daisy(const Daisy&) = delete;
  Daisy& operator=(const Daisy&) = delete;
  
  // default move semantics
  Daisy(Daisy&&);
  Daisy& operator=(Daisy&&);
  
  // single function this class exists for
  void oneDay_layEggs();
  void reset_Patch(const Rcpp::ListOf<Rcpp::List>& aTypes); 

};

#endif

/******************************************************************************
 * Multiplex multiLocus Class
 ******************************************************************************/
#ifndef MULTILOCUS_MPLEX
#define MULTILOCUS_MPLEX

class multiLocus: public Patch{
public:
  
  // constructor
  multiLocus(const int& patchID_,
              const Rcpp::ListOf<Rcpp::List>& aTypes,
              const Rcpp::ListOf<Rcpp::List>& maleReleases_,
              const Rcpp::ListOf<Rcpp::List>& femaleReleases_,
              const Rcpp::ListOf<Rcpp::List>& larvaeReleases_);
  // destructor
  virtual ~multiLocus();
  
  // delete copy constructor and assignment operator
  multiLocus(const multiLocus&) = delete;
  multiLocus& operator=(const multiLocus&) = delete;
  
  // default move semantics
  multiLocus(multiLocus&&);
  multiLocus& operator=(multiLocus&&);
  
  // single function this class exists for
  void oneDay_layEggs();
  void reset_Patch(const Rcpp::ListOf<Rcpp::List>& aTypes);
  
};

#endif

/******************************************************************************
 * Multiplex oneLocus Class
******************************************************************************/
#ifndef ONELOCUS_MPLEX
#define ONELOCUS_MPLEX

class oneLocus: public Patch{
public:
  
  // constructor
  oneLocus(const int& patchID_,
             const Rcpp::ListOf<Rcpp::List>& aTypes,
             const Rcpp::ListOf<Rcpp::List>& maleReleases_,
             const Rcpp::ListOf<Rcpp::List>& femaleReleases_,
             const Rcpp::ListOf<Rcpp::List>& larvaeReleases_);
  // destructor
  virtual ~oneLocus();
  
  // delete copy constructor and assignment operator
  oneLocus(const oneLocus&) = delete;
  oneLocus& operator=(const oneLocus&) = delete;
  
  // default move semantics
  oneLocus(oneLocus&&);
  oneLocus& operator=(oneLocus&&);
  
  // single function this class exists for
  void oneDay_layEggs();
  void reset_Patch(const Rcpp::ListOf<Rcpp::List>& aTypes);
  
};

#endif

/******************************************************************************
 * Mosquito Creation Functions
******************************************************************************/
// These are for initializing patches. It's ugly here, but I won't lose it.
#ifndef MOSQUITO_INITIALIZATION_MPLEX
#define MOSQUITO_INITIALIZATION_MPLEX

// for setting up Daisy and multiLocus patches
void CreateMosquitoes2Allele(int numMos, int minAge, dVec ageDist,
                             Rcpp::ListOf<Rcpp::List> aTypes, popVec returnPop);

// for setting up oneLocus patch
void CreateMosquitoes2Loci(int numMos, int minAge, dVec ageDist,
                           Rcpp::ListOf<Rcpp::List> aTypes, popVec returnPop);

// Daisy generating function
void DaisyOffspring(const std::string& fGen, const std::string& mGen);

// multiLocus Generating function
void MultiplexOffspring_mLoci(const std::string& fGen, const std::string& mGen);

// oneLocus Generating Function



#endif

