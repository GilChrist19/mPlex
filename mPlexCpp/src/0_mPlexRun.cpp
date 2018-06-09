

#include <RcppArmadillo.h>

#include <progress.hpp>
#include <progress_bar.hpp>

#include "1_Mosquito.hpp"
#include "2_Patch.hpp"
#include "3_DriveDefinitions.hpp"
#include "mPlex-Parameters.hpp"
#include "mPlex-Reference.hpp"
#include "mPlex-PRNG.hpp"




/******************************************************************************
 * Initialization function for mosquitoes
******************************************************************************/
popVec CreateMosquitoes(int numMos, int minAge, dVec ageDist, Rcpp::List aTypes){
  
  //I DON'T KNOW THAT THIS WORKS FOR ALL TYPES OF MULTIPLEXING, N
  // NEED TO CHECK LOGIC AGAIN!!!!!!!
  
  popVec returnPop(numMos);
  
  
  Rcpp::CharacterVector(Rcpp::as<Rcpp::List>(aTypes[0]).length() );
  std::vector<std::string> holder(2);
  
  
  int age;
  
  // loop over number of mosquitoes to create
  for(int count=0; count <= numMos; ++count){
    
    // generate each genotype
    for
    
    
    
    
    age = minAge + prng::instance().get_oneSample(ageDist);
    
    returnPop[count] = Mosquito(age, );
  } // end loop over mosquitoes
  
  
  // return population
  return( returnPop );
  
}


/******************************************************************************
 * Run the simulation
******************************************************************************/
//' Run Multi-Plex Experiment
//' 
//' Run mPlex experiments, choosing the type of inheritance desired.
//' 
//' @examples
//' nope
//' none
//' sorry
//' 
// [[Rcpp::export]]
void run_mPlex_repetitions(){
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}
