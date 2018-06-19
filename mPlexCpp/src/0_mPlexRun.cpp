



#include <progress.hpp>
#include <progress_bar.hpp>

#include <Rcpp.h>

#include "1_Mosquito.hpp"
#include "2_Patch.hpp"
#include "3_DriveDefinitions.hpp"
#include "4_Parameters.hpp"
#include "4_PRNG.hpp"
#include "4_Reference.hpp"


/******************************************************************************
 * Run the simulation
******************************************************************************/
//' Run Multi-Plex Experiment
//' 
//' Run mPlex experiments, choosing the type of inheritance desired.
//' 
//' @examples
//' \dontrun{
//' nope
//' none
//' sorry
//' }
//' 
//' 
// [[Rcpp::export]]
void run_mPlex_Cpp(const uint_least32_t& seed,
                   const Rcpp::List& networkParameters_,
                   const Rcpp::List& reproductionReference_,
                   const Rcpp::NumericMatrix& migrationMale_,
                   const Rcpp::NumericMatrix& migrationFemale_,
                   const Rcpp::List& migrationBatch_,
                   const std::string& reproductionType_,
                   const bool& verbose){
  
  
  
  
  
  // BEGIN INITIALIZE PRNG
  if(verbose) {Rcpp::Rcout << "Initializing prng ... ";};
  prng::instance().set_seed(seed);
  if(verbose){Rcpp::Rcout <<  " ... done initializing prng!\n";};
  // END INITIALIZE PRNG
  
  
  
  
  // BEGIN SET REFERENCE
  // check reproduction type
  if((reproductionType_ != "DaisyDrive") && (reproductionType_ != "mPlex_oLocus") && (reproductionType_ != "mPlex_mLoci")){
    Rcpp::stop("\nreproductionType must match one of these choices:\n DaisyDrive\n mPlex_oLocus\n mPlex_mLoci\n");
  }
    
  if(verbose) {Rcpp::Rcout << "Setting reference ... ";};
  // set genotype specific parameters than all drives have
  reference::instance().set_reference(reproductionReference_["eta"], reproductionReference_["phi"],
                                      reproductionReference_["omega"], reproductionReference_["xiF"],
                                      reproductionReference_["xiM"], reproductionReference_["s"]);
  
  Rcpp::Rcout << "\nBasics Done\n";
  
  // set mendelian allele reference
  reference::instance().set_mendelian(Rcpp::as<Rcpp::List> (reproductionReference_["mendelian"]),
                                      Rcpp::as<Rcpp::List> (reproductionReference_["mendelianAlleles"]));
  
  Rcpp::Rcout << "\nMendelian Done\n";
  
  // set homing allele reference
  reference::instance().set_homing(Rcpp::as<Rcpp::List> (reproductionReference_["homing"]),
                                    Rcpp::as<Rcpp::List> (reproductionReference_["homingAlleles"]));
  
  Rcpp::Rcout << "\nHoming Done\n";
  

  // set cutting reference if daisy drive
  if(reproductionType_ == "DaisyDrive"){
    reference::instance().set_cutting(Rcpp::as<Rcpp::List> (reproductionReference_["cutting"]),
                                      Rcpp::as<Rcpp::List> (reproductionReference_["cuttingAlleles"]));
  }
  
  if(verbose){Rcpp::Rcout <<  " ... done setting reference!\n";};
  // END SET REFERENCE
  
  
  
  
  
  
  
  // 
  // 
  // // BEGIN SET PARAMETERS
  // if(verbose) {Rcpp::Rcout << "Setting parameters ... ";};
  // 
  // // setup input matrices and then fill them. They are std::vectors, not Rcpp things
  // int dim = migrationMale_.ncol();
  // Rcpp::NumericVector hold(dim);
  // dMat mHold, fHold, bsHold, bmHold;
  // 
  // // loop over number of dims
  // for(size_t i=0; i<dim; ++i){
  //   // male migration
  //   //  multiply move var by matrix now, since this is all we ever do with it.
  //   hold = networkParameters_["moveVar"] * migrationMale_.row(i);
  //   mHold.push_back(Rcpp::as<dVec>(hold));
  //   // female migration
  //   hold = networkParameters_["moveVar"] * migrationFemale_.row(i);
  //   fHold.push_back(Rcpp::as<dVec>(hold));
  //   // batch sex discrepancy
  //   hold = Rcpp::as<Rcpp::NumericMatrix>(migrationBatch_["sexProbs"]).row(i);
  //   bsHold.push_back(Rcpp::as<dVec>(hold));
  //   // batch migration
  //   hold = Rcpp::as<Rcpp::NumericMatrix>(migrationBatch_["moveProbs"]).row(i);
  //   bmHold.push_back(Rcpp::as<dVec>(hold));
  // }
  // 
  // // now set parameters
  // parameters::instance().set_parameters(networkParameters_["nPatch"],networkParameters_["simTime"],networkParameters_["runID"],
  //                      networkParameters_["stageTime"],networkParameters_["beta"],networkParameters_["mu"],
  //                      networkParameters_["alpha"],networkParameters_["Leq"],networkParameters_["AdPopEQ"],
  //                      mHold, fHold, migrationBatch_["batchProbs"], bsHold, bmHold);
  // 
  // if(verbose){Rcpp::Rcout <<  " ... done setting parameters!\n";};
  // // END SET PARAMETERS
  // 
  // 
  

  
  
  
  
  // BEGIN INITIALIZE LOGGER
  if(verbose) {Rcpp::Rcout << "Initializing logging ... ";};
  
  if(verbose){Rcpp::Rcout <<  " ... done initializing logging!\n";};
  // END INITIALIZE LOGGER
  
  
  
  
  
  // BEGIN INITIALIZE NETWORK
  
  
  // END INITIALIZE NETWORK
  
  
  
  
  
  
  
  
  
  
}
