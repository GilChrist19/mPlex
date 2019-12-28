///////////////////////////////////////////////////////////////////////////////
//                           ____  __          ______          
//                ____ ___  / __ \/ /__  _  __/ ____/___  ____ 
//               / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
//              / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
//             /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/ 
//                                              /_/   /_/      
///////////////////////////////////////////////////////////////////////////////

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <progress.hpp>
#include <progress_bar.hpp>

//#include <gperftools/profiler.h>

#include <omp.h> // for parallel loops
#include "1_Mosquito.hpp"
#include "2_Patch.hpp"
#include "3_DriveDefinitions.hpp"
#include "4_Parameters.hpp"
#include "4_PRNG.hpp"
#include "4_Reference.hpp"
#include "4_BigBrother.hpp"


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
void run_mPlex_Cpp(const uint_least32_t& seed_,
                   const uint_least32_t& numThreads_,
                   const Rcpp::List& networkParameters_,
                   const Rcpp::List& reproductionReference_,
                   const Rcpp::List& patchReleases_,
                   const Rcpp::NumericMatrix& migrationMale_,
                   const Rcpp::NumericMatrix& migrationFemale_,
                   const Rcpp::List& migrationBatch_,
                   const std::string& outputDirectory_,
                   const std::string& reproductionType_,
                   const bool& verbose_){
  
  
  #ifdef BASE_PROFILER_H_
    // make sure to change path!!! But keep file name.
    ProfilerStart("/home/jared/Desktop/OUTPUT/profile.log");
  #endif
  
  
  ////////////////////
  // BEGIN INITIALIZE PRNG
  ////////////////////
  if(verbose_) {Rcpp::Rcout << "Initializing " << numThreads_ << " prngs ... ";};

  int myThread;
  omp_set_num_threads(numThreads_);
  
  std::vector<prng> randInst;
  randInst.reserve(numThreads_);
  for(auto i=0; i < numThreads_; i++){
    randInst.emplace_back(prng(i + seed_));
  }
  
  if(verbose_){Rcpp::Rcout <<  " ... done initializing prngs!\n";};
  ////////////////////
  // END INITIALIZE PRNG
  ////////////////////
  
  
  
  ////////////////////
  // BEGIN SET REFERENCE
  ////////////////////
  // check reproduction type
  if((reproductionType_ != "DaisyDrive") && (reproductionType_ != "mPlex_oLocus")
       && (reproductionType_ != "mPlex_mLoci") && (reproductionType_ != "Family")){
    Rcpp::stop("\nreproductionType must match one of these choices:\n DaisyDrive\n mPlex_oLocus\n mPlex_mLoci\n Family\n");
  }
    
  if(verbose_) {Rcpp::Rcout << "Setting reference ... ";};
  // set genotype specific parameters than all drives have
  reference::instance().set_reference(reproductionReference_["eta"], reproductionReference_["phi"],
                                      reproductionReference_["omega"], reproductionReference_["xiF"],
                                      reproductionReference_["xiM"], reproductionReference_["s"]);
  if(verbose_) {Rcpp::Rcout << "\n\tBasics Done\n";};
  
  if(reproductionType_ != "Family"){
    // set allele types for initialization etc
    reference::instance().set_alleleTypes(Rcpp::as<Rcpp::ListOf<Rcpp::List> >(networkParameters_["alleloTypes"]));
    
    // set mendelian allele reference
    reference::instance().set_mendelian(Rcpp::as<Rcpp::List> (reproductionReference_["mendelian"]),
                                        Rcpp::as<Rcpp::List> (reproductionReference_["mendelianAlleles"]));
    if(verbose_) {Rcpp::Rcout << "\tMendelian Done\n";};
    
    // set homing allele reference
    reference::instance().set_homing(Rcpp::as<Rcpp::List> (reproductionReference_["homing"]),
                                      Rcpp::as<Rcpp::List> (reproductionReference_["homingAlleles"]));
    if(verbose_) {Rcpp::Rcout << "\tHoming Done\n";};
  }

  // set cutting reference if daisy drive
  if(reproductionType_ == "DaisyDrive"){
    reference::instance().set_cutting(Rcpp::as<Rcpp::List> (reproductionReference_["cutting"]),
                                      Rcpp::as<Rcpp::List> (reproductionReference_["cuttingAlleles"]));
    if(verbose_) {Rcpp::Rcout << "\tDaisy set\n";};
  }
  
  if(verbose_){Rcpp::Rcout <<  " ... done setting reference!\n";};
  ////////////////////
  // END SET REFERENCE
  ////////////////////
  
  
  
  ////////////////////
  // BEGIN SET PARAMETERS
  ////////////////////
  if(verbose_) {Rcpp::Rcout << "Setting parameters ... ";};

  // setup input matrices and then fill them. They are std::vectors, not Rcpp things
  int dim = migrationMale_.ncol();
  Rcpp::NumericVector hold(dim);
  dMat mHold, fHold, bsHold, bmHold;

  // loop over number of dims
  for(size_t i=0; i<dim; ++i){
    // male migration
    //  multiply move var by matrix now, since this is all we ever do with it.
    hold = networkParameters_["moveVar"] * migrationMale_.row(i);
    mHold.push_back(Rcpp::as<dVec>(hold));
    // female migration
    hold = networkParameters_["moveVar"] * migrationFemale_.row(i);
    fHold.push_back(Rcpp::as<dVec>(hold));
    // batch sex discrepancy
    hold = Rcpp::as<Rcpp::NumericMatrix>(migrationBatch_["sexProbs"]).row(i);
    bsHold.push_back(Rcpp::as<dVec>(hold));
    // batch migration
    hold = Rcpp::as<Rcpp::NumericMatrix>(migrationBatch_["moveProbs"]).row(i);
    bmHold.push_back(Rcpp::as<dVec>(hold));
  }

  // now set parameters
  parameters::instance().set_parameters(networkParameters_["nPatch"],networkParameters_["simTime"],networkParameters_["runID"],
                       networkParameters_["stageTime"],networkParameters_["beta"],networkParameters_["mu"],
                       networkParameters_["alpha"],networkParameters_["Leq"],networkParameters_["AdPopEQ"],
                       mHold, fHold, migrationBatch_["batchProbs"], bsHold, bmHold);

  if(verbose_){Rcpp::Rcout <<  " ... done setting parameters!\n";};
  ////////////////////
  // END SET PARAMETERS
  ////////////////////
  
  
  
  ////////////////////
  // BEGIN INITIALIZE NETWORK
  ////////////////////
  if(verbose_){Rcpp::Rcout << "Initializing network ... \n";};
  
  size_t numPatches = parameters::instance().get_n_patch();

  // vector of patches
  std::vector<std::unique_ptr<Patch> > patches;
  patches.reserve(numPatches);
  
  Rcpp::List patchRelease;
  
  // setup all patches
  for(size_t np=0; np<numPatches; ++np){
    
    patchRelease = patchReleases_[np];
    
    // make patches of correct child type
    if(reproductionType_ == "DaisyDrive"){
      
      if(verbose_ && np==0){Rcpp::Rcout << "\tInitializing Daisy Drive"<<std::endl;};
      
      patches.emplace_back(std::make_unique<Daisy>(np,
                                                   patchRelease["maleReleases"],
                                                   patchRelease["femaleReleases"],
                                                   patchRelease["eggReleases"]));
      
    } else if(reproductionType_ == "mPlex_oLocus"){
      
      if(verbose_ && np==0){Rcpp::Rcout << "\tInitializing oneLocus Drive"<<std::endl;};
      
      patches.emplace_back(std::make_unique<oneLocus>(np,
                                                   patchRelease["maleReleases"],
                                                   patchRelease["femaleReleases"],
                                                   patchRelease["eggReleases"]));
      
    } else if(reproductionType_ == "mPlex_mLoci"){
      
      if(verbose_ && np==0){Rcpp::Rcout << "\tInitializing multiLocus Drive"<<std::endl;};
      
      patches.emplace_back(std::make_unique<multiLocus>(np,
                                                   patchRelease["maleReleases"],
                                                   patchRelease["femaleReleases"],
                                                   patchRelease["eggReleases"]));
      
    } else if(reproductionType_ == "Family"){
      
      if(verbose_ && np==0){Rcpp::Rcout << "\tBigBrother is watching"<<std::endl;};
      
      patches.emplace_back(std::make_unique<Family>(np,
                                                    patchRelease["maleReleases"],
                                                    patchRelease["femaleReleases"],
                                                    patchRelease["eggReleases"]));
      
    }
    
  } // end network loop
  
  if(verbose_){Rcpp::Rcout <<  " ... done initializing network!\n";};
  ////////////////////
  // END INITIALIZE NETWORK
  ////////////////////
  
  

  ////////////////////
  // BEGIN INITIALIZE LOGGER
  ////////////////////
  if(verbose_) {Rcpp::Rcout << "Initializing logging ... ";};
  
  // setup vectors of ofstreams
  std::vector<std::ofstream> M_output(numPatches), F_output(numPatches);
  std::string sHold;

  // setup strings for file names
  std::string maleFile(outputDirectory_+"/M_Run_"
                       + std::string(3 - std::to_string(parameters::instance().get_run_id()).length(), '0')
                       + std::to_string(parameters::instance().get_run_id())
                       + "_Patch_"),
              femaleFile(outputDirectory_+"/F_Run_"
                       + std::string(3 - std::to_string(parameters::instance().get_run_id()).length(), '0')
                       + std::to_string(parameters::instance().get_run_id())
                       + "_Patch_");
  
  // open all streams with file names
  #pragma omp parallel for private(sHold)
  for(size_t np=0; np<numPatches; ++np){
    // denote patch
    sHold = std::string(3 - std::to_string(np).length(), '0')
            + std::to_string(np)
            + ".csv";
    // open streams
    M_output[np].open(maleFile + sHold);
    F_output[np].open(femaleFile + sHold);
  }

  if(verbose_){Rcpp::Rcout <<  " ... done initializing logging!\n";};
  ////////////////////
  // END INITIALIZE LOGGER
  ////////////////////
  

  
  ////////////////////
  // BEGIN INITIALIZE OUTPUT
  ////////////////////
  if(verbose_){Rcpp::Rcout <<  "Initializing output ... ";};
  
  // output is to different files, so parallel
  #pragma omp parallel for
  for(size_t np=0; np<numPatches; ++np){
    patches[np]->init_output(M_output[patches[np]->get_patchID()], F_output[patches[np]->get_patchID()]);
  }
  
  // increment time to begin
  parameters::instance().increment_t_now();
  if(verbose_){Rcpp::Rcout <<  " ... done initializing output!\n";};
  ////////////////////
  // END INITIALIZE OUTPUT
  ////////////////////
  
  

  ////////////////////
  // BEGIN SIMULATION
  ////////////////////
  if(verbose_){Rcpp::Rcout << "begin simulation ... \n";};
  int tMax = parameters::instance().get_sim_time();
  Progress pb(tMax-1,verbose_);
  
  while(parameters::instance().get_t_now() < tMax){
    
    // test for user interrupt
    if(checkInterrupt()) return;
    
    // Independent daily operations
    #pragma omp parallel for private(myThread)
    for(size_t np=0; np<numPatches; ++np){
      // get unique thread for prng
      myThread = omp_get_thread_num();
      // run today's stuff
      patches[np]->oneDay_popDynamics(randInst[myThread]);
    }
    
    // Do in-bound migration in one large loop, no extra structures
    for(size_t inPatch=0; inPatch < numPatches; ++inPatch){
      for(size_t outPatch=0; outPatch < numPatches; ++outPatch){
          patches[inPatch]->oneDay_migrationIn(patches[outPatch]->get_maleMigration(inPatch),
                                               patches[outPatch]->get_femaleMigration(inPatch));
      }
    }
    
    // Log output
    #pragma omp parallel for
    for(size_t np=0; np<numPatches; ++np){
      // run output
      patches[np]->oneDay_writeOutput(M_output[patches[np]->get_patchID()], F_output[patches[np]->get_patchID()]);
    }
    
    // increment time and progress
    parameters::instance().increment_t_now();
    pb.increment();
  }
  
  if(verbose_){Rcpp::Rcout << "... simulation done!\n";};
  ////////////////////
  // END SIMULATION
  ////////////////////
  
  
  // close files
  // all independent, do in parallel
  #pragma omp parallel for
  for(size_t np=0; np<numPatches; ++np){
    M_output[np].close();
    F_output[np].close();
  }
  
  
  // close profiler
  #ifdef BASE_PROFILER_H_
    ProfilerStop();
  #endif
  
}











//' Run Multi-Plex Experiment
//'
//' Run multiple mPlex experiments, choosing the type of inheritance desired.
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
void run_mPlex_Cpp_repetitions(const uint_least32_t& seed_,
                               const uint_least32_t& numReps_,
                               const uint_least32_t& numThreads_,
                               const Rcpp::List& networkParameters_,
                               const Rcpp::List& reproductionReference_,
                               const Rcpp::List& patchReleases_,
                               const Rcpp::NumericMatrix& migrationMale_,
                               const Rcpp::NumericMatrix& migrationFemale_,
                               const Rcpp::List& migrationBatch_,
                               const std::string& outputDirectory_,
                               const std::string& reproductionType_,
                               const bool& verbose_){

  #ifdef BASE_PROFILER_H_
    // make sure to change path!!! But keep file name.
    ProfilerStart("/home/jared/Desktop/OUTPUT/profile.log");
  #endif

    
  ////////////////////
  // BEGIN INITIALIZE PRNG
  ////////////////////
  if(verbose_) {Rcpp::Rcout << "Initializing " << numThreads_ << " prngs ... ";};

  int myThread;
  omp_set_num_threads(numThreads_);
  
  std::vector<prng> randInst;
  randInst.reserve(numThreads_);
  for(auto i=0; i < numThreads_; i++){
    randInst.emplace_back(prng(i + seed_));
  }
  
  if(verbose_){Rcpp::Rcout <<  " ... done initializing prngs!\n";};
  ////////////////////
  // END INITIALIZE PRNG
  ////////////////////



  ////////////////////
  // BEGIN SET REFERENCE
  ////////////////////
  // check reproduction type
  if((reproductionType_ != "DaisyDrive") && (reproductionType_ != "mPlex_oLocus")
       && (reproductionType_ != "mPlex_mLoci") && (reproductionType_ != "Family")){
    Rcpp::stop("\nreproductionType must match one of these choices:\n DaisyDrive\n mPlex_oLocus\n mPlex_mLoci\n Family\n");
  }
    
  if(verbose_) {Rcpp::Rcout << "Setting reference ... ";};
  // set genotype specific parameters than all drives have
  reference::instance().set_reference(reproductionReference_["eta"], reproductionReference_["phi"],
                                      reproductionReference_["omega"], reproductionReference_["xiF"],
                                      reproductionReference_["xiM"], reproductionReference_["s"]);
  if(verbose_) {Rcpp::Rcout << "\n\tBasics Done\n";};
  
  if(reproductionType_ != "Family"){
    // set allele types for initialization etc
    reference::instance().set_alleleTypes(Rcpp::as<Rcpp::ListOf<Rcpp::List> >(networkParameters_["alleloTypes"]));
    
    // set mendelian allele reference
    reference::instance().set_mendelian(Rcpp::as<Rcpp::List> (reproductionReference_["mendelian"]),
                                        Rcpp::as<Rcpp::List> (reproductionReference_["mendelianAlleles"]));
    if(verbose_) {Rcpp::Rcout << "\tMendelian Done\n";};
    
    // set homing allele reference
    reference::instance().set_homing(Rcpp::as<Rcpp::List> (reproductionReference_["homing"]),
                                      Rcpp::as<Rcpp::List> (reproductionReference_["homingAlleles"]));
    if(verbose_) {Rcpp::Rcout << "\tHoming Done\n";};
  }

  // set cutting reference if daisy drive
  if(reproductionType_ == "DaisyDrive"){
    reference::instance().set_cutting(Rcpp::as<Rcpp::List> (reproductionReference_["cutting"]),
                                      Rcpp::as<Rcpp::List> (reproductionReference_["cuttingAlleles"]));
    if(verbose_) {Rcpp::Rcout << "\tDaisy set\n";};
  }
  
  if(verbose_){Rcpp::Rcout <<  " ... done setting reference!\n";};
  ////////////////////
  // END SET REFERENCE
  ////////////////////



  ////////////////////
  // BEGIN SET PARAMETERS
  ////////////////////
  if(verbose_) {Rcpp::Rcout << "Setting parameters ... ";};

  // setup input matrices and then fill them. They are std::vectors, not Rcpp things
  int dim = migrationMale_.ncol();
  Rcpp::NumericVector hold(dim);
  dMat mHold, fHold, bsHold, bmHold;

  // loop over number of dims
  for(size_t i=0; i<dim; ++i){
    // male migration
    //  multiply move var by matrix now, since this is all we ever do with it.
    hold = networkParameters_["moveVar"] * migrationMale_.row(i);
    mHold.push_back(Rcpp::as<dVec>(hold));
    // female migration
    hold = networkParameters_["moveVar"] * migrationFemale_.row(i);
    fHold.push_back(Rcpp::as<dVec>(hold));
    // batch sex discrepancy
    hold = Rcpp::as<Rcpp::NumericMatrix>(migrationBatch_["sexProbs"]).row(i);
    bsHold.push_back(Rcpp::as<dVec>(hold));
    // batch migration
    hold = Rcpp::as<Rcpp::NumericMatrix>(migrationBatch_["moveProbs"]).row(i);
    bmHold.push_back(Rcpp::as<dVec>(hold));
  }

  // now set parameters
  parameters::instance().set_parameters(networkParameters_["nPatch"],networkParameters_["simTime"],networkParameters_["runID"],
                       networkParameters_["stageTime"],networkParameters_["beta"],networkParameters_["mu"],
                       networkParameters_["alpha"],networkParameters_["Leq"],networkParameters_["AdPopEQ"],
                       mHold, fHold, migrationBatch_["batchProbs"], bsHold, bmHold);

  if(verbose_){Rcpp::Rcout <<  " ... done setting parameters!\n";};
  ////////////////////
  // END SET PARAMETERS
  ////////////////////

  

  ////////////////////
  // BEGIN INITIALIZE NETWORK
  ////////////////////
  if(verbose_){Rcpp::Rcout << "Initializing network ... \n";};
  
  size_t numPatches = parameters::instance().get_n_patch();

  std::vector<std::unique_ptr<Patch> > patches;
  patches.reserve(numPatches);
  
  Rcpp::List patchRelease;
  
  // setup all patches
  for(size_t np=0; np<numPatches; ++np){
    
    patchRelease = patchReleases_[np];
    
    // make patches of correct child type
    if(reproductionType_ == "DaisyDrive"){
      
      if(verbose_ && np==0){Rcpp::Rcout << "\tInitializing Daisy Drive"<<std::endl;};
      
      patches.emplace_back(std::make_unique<Daisy>(np,
                                                   patchRelease["maleReleases"],
                                                   patchRelease["femaleReleases"],
                                                   patchRelease["eggReleases"]));
      
    } else if(reproductionType_ == "mPlex_oLocus"){
      
      if(verbose_ && np==0){Rcpp::Rcout << "\tInitializing oneLocus Drive"<<std::endl;};
      
      patches.emplace_back(std::make_unique<oneLocus>(np,
                                                   patchRelease["maleReleases"],
                                                   patchRelease["femaleReleases"],
                                                   patchRelease["eggReleases"]));
      
    } else if(reproductionType_ == "mPlex_mLoci"){
      
      if(verbose_ && np==0){Rcpp::Rcout << "\tInitializing multiLocus Drive"<<std::endl;};
      
      patches.emplace_back(std::make_unique<multiLocus>(np,
                                                   patchRelease["maleReleases"],
                                                   patchRelease["femaleReleases"],
                                                   patchRelease["eggReleases"]));
      
    } else if(reproductionType_ == "Family"){
      
      if(verbose_ && np==0){Rcpp::Rcout << "\tBigBrother is watching"<<std::endl;};
      
      patches.emplace_back(std::make_unique<Family>(np,
                                                    patchRelease["maleReleases"],
                                                    patchRelease["femaleReleases"],
                                                    patchRelease["eggReleases"]));
      
    }
    
  } // end network loop
  
  if(verbose_){Rcpp::Rcout <<  " ... done initializing network!\n";};
  ////////////////////
  // END INITIALIZE NETWORK
  ////////////////////



  
  
  

  // initialize things for inside the loop
  // std::string maleFile, femaleFile;
  // std::ofstream M_output, F_output;
  int tMax = parameters::instance().get_sim_time();
  
  std::vector<std::ofstream> M_output(numPatches), F_output(numPatches);
  std::string maleFile, femaleFile, sHold;



  ////////////////////
  // BEGIN REPETITION WRAP
  ////////////////////
  for(size_t rep=0; rep<numReps_; ++rep){
    
    if(verbose_){Rcpp::Rcout <<  "begin repetition " << rep+parameters::instance().get_run_id()<< "\n";};

    ////////////////////
    // BEGIN INITIALIZE LOGGER
    ////////////////////
    if(verbose_) {Rcpp::Rcout << "Initializing logging ... ";};

    // setup strings for file names
    //  these increment runs
    maleFile = outputDirectory_ + "/M_Run_"
                + std::string(3 - std::to_string(rep+parameters::instance().get_run_id()).length(), '0')
                + std::to_string(rep+parameters::instance().get_run_id())
                + "_Patch_";
    femaleFile = outputDirectory_ + "/F_Run_"
                + std::string(3 - std::to_string(rep+parameters::instance().get_run_id()).length(), '0')
                + std::to_string(rep+parameters::instance().get_run_id())
                + "_Patch_";
    
    // open all streams with file names
    // try in parallel. It's independent, but idk if this is a speedup or just a waste
    #pragma omp parallel for private(sHold)
    for(size_t np=0; np<numPatches; ++np){
      // denote patch
      sHold = std::string(3 - std::to_string(np).length(), '0')
              + std::to_string(np)
              + ".csv";
      // open streams
      M_output[np].open(maleFile + sHold);
      F_output[np].open(femaleFile + sHold);
    }
    
    if(verbose_){Rcpp::Rcout <<  " ... done initializing logging!\n";};
    /////////////////////
    // END INITIALIZE LOGGER
    ////////////////////



    ////////////////////
    // BEGIN INITIALIZE OUTPUT
    ////////////////////
    if(verbose_){Rcpp::Rcout <<  "Initializing output ... ";};
    
    // output is to different files, so parallel
    #pragma omp parallel for
    for(size_t np=0; np<numPatches; ++np){
      patches[np]->init_output(M_output[patches[np]->get_patchID()], F_output[patches[np]->get_patchID()]);
    }
    
    // increment time to begin
    parameters::instance().increment_t_now();
    if(verbose_){Rcpp::Rcout <<  " ... done initializing output!\n";};
    ////////////////////
    // END INITIALIZE OUTPUT
    ////////////////////

    
    ////////////////////
    // BEGIN SIMULATION
    ////////////////////
    Progress pb(tMax-1,verbose_);
    while(parameters::instance().get_t_now() < tMax){

      // test for user interrupt
      if(checkInterrupt()) return;

      // Independent daily operations
      #pragma omp parallel for private(myThread)
      for(size_t np=0; np<numPatches; ++np){
        // get unique thread for prng
        myThread = omp_get_thread_num();
        // run today's stuff
        patches[np]->oneDay_popDynamics(randInst[myThread]);
      }

      // Do in-bound migration in one large loop, no extra structures
      for(size_t inPatch=0; inPatch < numPatches; ++inPatch){
        for(size_t outPatch=0; outPatch < numPatches; ++outPatch){
          patches[inPatch]->oneDay_migrationIn(patches[outPatch]->get_maleMigration(inPatch),
                                               patches[outPatch]->get_femaleMigration(inPatch));
        }
      }
      
      // Log output
      #pragma omp parallel for
      for(size_t np=0; np<numPatches; ++np){
        patches[np]->oneDay_writeOutput(M_output[patches[np]->get_patchID()], F_output[patches[np]->get_patchID()]);
      }

      // increment time and progress
      parameters::instance().increment_t_now();
      pb.increment();
    } // end one sim loop
    ////////////////////
    // END SIMULATION
    ////////////////////
    
    
    // close files
    // all independent, do in parallel
    #pragma omp parallel for
    for(size_t np=0; np<numPatches; ++np){
      M_output[np].close();
      F_output[np].close();
    }
    
    
    ////////////////////
    // RESET SIMULATION
    ////////////////////
    //  This does involve rebuilding the starting mosquitoes
    parameters::instance().reset_t_now();
    BigBrother::instance().reset();
    #pragma omp parallel for
    for(size_t np=0; np<numPatches; ++np){
      patches[np]->reset_Patch();
    }

    if(verbose_){Rcpp::Rcout <<  "end repetition " << rep+parameters::instance().get_run_id() << "\n\n";};
  } // end repetition loop

  if(verbose_){Rcpp::Rcout << "... end repetitions \n";};
  ////////////////////
  // END REPETITION WRAP
  ////////////////////
  
  
  // close profiler
  #ifdef BASE_PROFILER_H_
    ProfilerStop();
  #endif
}




