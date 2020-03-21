///////////////////////////////////////////////////////////////////////////////
//     ________ __ __  _______ 
//    / ____/ //_//  |/  / __ \
//   / /   / ,<  / /|_/ / /_/ /
//  / /___/ /| |/ /  / / _, _/ 
//  \____/_/ |_/_/  /_/_/ |_|  
//      
///////////////////////////////////////////////////////////////////////////////
#ifndef PATCH_MPLEX
#define PATCH_MPLEX

#include <RcppArmadillo.h>

#include <vector>
#include <fstream>
#include <math.h>

#include "1_Mosquito.hpp"
#include "4_Reference.hpp"
#include "4_PRNG.hpp"
#include "4_Parameters.hpp"
#include "4_bigBrother.hpp"


/**************************************
 * forward declarations & alias
**************************************/
using popVec = std::vector<Mosquito>;
using dVec = std::vector<double>;
using iVec = std::vector<int>; 
using sVec = std::vector<std::string>;


/**************************************
 * single release event struct
**************************************/
using release_event = struct release_event {

  // three things needed, genotype, age, and release time
  sVec  pop_names;
  iVec  pop_ages;
  int   release_time;

  // constructor and destructor
  release_event(const sVec& pop_names_, const iVec& pop_ages_, const int& release_time_) : 
    pop_names(pop_names_), pop_ages(pop_ages_), release_time(release_time_) {};
  ~release_event() {};

};

/**************************************
 * patch class declaration
**************************************/
class Patch {

public:

  /* constructor & destructor */
  Patch(const int& patchID_,
        const Rcpp::List& maleReleases_,
        const Rcpp::List& femaleReleases_,
        const Rcpp::List& eggReleases_);
  virtual ~Patch();

  /* delete all copy semantics: ensures we get legible compile-time errors if we do something stupid */
  Patch(const Patch&) = delete;
  Patch& operator=(const Patch&) = delete;

  /* default move semantics */
  Patch(Patch&&);
  Patch& operator=(Patch&&);

  
  // getters
  int       get_patchID(){return patchID;};
  
  popVec    get_eggs(){return eggs;};
  popVec    get_larva(){return larva;};
  popVec    get_pupa(){return pupa;};
  popVec    get_adult_male(){return adult_male;};
  popVec    get_adult_female(){return adult_female;};
  popVec    get_unmated_female(){return unmated_female;};
  
  popVec    get_maleMigration(const size_t& patch){return maleMigration[patch];};
  popVec    get_femaleMigration(const size_t& patch){return femaleMigration[patch];};
  
  
  // Population dynamics
  void      oneDay_popDynamics(prng& myPRNG, bigBrother& myBB);
  
  // Death
  void      oneDay_eggDeathAge(prng& myPRNG);
  void      oneDay_larvaDeathAge(prng& myPRNG);
  void      oneDay_pupaDeathAge(prng& myPRNG);
  void      oneDay_adultDeathAge(prng& myPRNG);
  
  // Maturation
  void      oneDay_pupaMaturation(prng& myPRNG);
  void      oneDay_larvaMaturation();
  void      oneDay_eggMaturation();
  
  // Mating
  void      oneDay_mating(prng& myPRNG);
  
  // New Eggs
  virtual void  oneDay_layEggs(prng& myPRNG, bigBrother& myBB) = 0;
  
  // Releases 
  void      oneDay_Releases();

  // migration
  void      oneDay_migrationOut(prng& myPRNG);
  void      oneDay_migrationIn(const popVec& male, const popVec& female);
  
  // extras
  virtual void  reset_Patch(bigBrother& myBB) = 0;
  
  void  init_output(std::vector<std::ofstream *>& logFiles);
  void  oneDay_writeOutput(std::vector<std::ofstream *>& logFiles, prng& myPRNG);

  
protected:

  // patchy things
  int patchID;

  // store current populations
  popVec eggs;
  popVec larva;
  popVec pupa;
  popVec adult_male;
  popVec adult_female;
  popVec unmated_female;
  
  // migration things
  std::vector<popVec> maleMigration;
  std::vector<popVec> femaleMigration;
  dVec                probsMigration;
  
  // mating things
  std::vector<std::string> genNames;
  std::vector<double> genProbs;
  std::string mateName;
  
  // holder things
  double holdDbl;
  int    holdInt;
  
  
  // releases
  std::vector<release_event> releaseM0; //male initial releases
  std::vector<release_event> releaseF0; //female initial releases
  std::vector<release_event> releaseE0; //larval initial releases
  
  std::vector<release_event> releaseM; //male releases
  std::vector<release_event> releaseF; //female releases
  std::vector<release_event> releaseE; //larval releases
  
};

#endif
