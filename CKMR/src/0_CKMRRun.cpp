///////////////////////////////////////////////////////////////////////////////
//     ________ __ __  _______ 
//    / ____/ //_//  |/  / __ \
//   / /   / ,<  / /|_/ / /_/ /
//  / /___/ /| |/ /  / / _, _/ 
//  \____/_/ |_/_/  /_/_/ |_|  
//   
///////////////////////////////////////////////////////////////////////////////

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <progress.hpp>
#include <progress_bar.hpp>

//#include <gperftools/profiler.h>

// I think the plugin includes compiler flags
// Source: https://wbnicholson.wordpress.com/2014/07/10/parallelization-in-rcpp-via-openmp/
//  However, do I even need it? I have most of the flags
#ifdef _OPENMP
#include <omp.h> // for parallel loops
#endif
// [[Rcpp::plugins(openmp)]]

#include "1_Mosquito.hpp"
#include "2_Patch.hpp"
#include "3_DriveDefinitions.hpp"
#include "4_Parameters.hpp"
#include "4_PRNG.hpp"
#include "4_Reference.hpp"
#include "4_bigBrother.hpp"


/******************************************************************************
 * Run the simulation
******************************************************************************/
//' Run CKMR Experiment
//'
//' Run multiple CKMR experiments.
//'
//' @examples
//' \dontrun{
//' nope
//' none
//' sorry
//' }
//'
// [[Rcpp::export]]
void run_CKMR(const std::uint64_t& s1_,
               const std::uint64_t& s2_,
               const std::uint64_t& s3_,
               const std::uint64_t& s4_,
               const uint_least32_t& numReps_,
               const uint_least32_t& numThreads_,
               const Rcpp::List& networkParameters_,
               const Rcpp::List& reproductionReference_,
               const Rcpp::List& patchReleases_,
               const Rcpp::NumericMatrix& migrationMale_,
               const Rcpp::NumericMatrix& migrationFemale_,
               const Rcpp::List& migrationBatch_,
               const Rcpp::List& samplingParameters_,
               const std::string& outputDirectory_,
               const bool& verbose_){

  #ifdef BASE_PROFILER_H_
    // make sure to change path!!! But keep file name.
    ProfilerStart("/home/jared/Desktop/OUTPUT/profile.log");
  #endif


  ////////////////////
  // BEGIN INITIALIZE PRNG
  ////////////////////
  if(verbose_) {Rcpp::Rcout << "Initializing " << numThreads_ << " prngs ... ";};

  #ifdef _OPENMP
  int myThread;
  omp_set_num_threads(numThreads_);
  #endif
  
  std::vector<prng> randInst;
  randInst.reserve(numThreads_);
  std::array<std::uint64_t, 4> seed = {s1_,s2_,s3_,s4_};
  for(size_t i=0; i < numThreads_; i++){
    // iterate all seeds in vector
    for(auto& j : seed){j += i;}
    // set seed
    randInst.emplace_back(prng(seed));
  }
  
  if(verbose_){Rcpp::Rcout <<  " ... done initializing prngs!\n";};
  ////////////////////
  // END INITIALIZE PRNG
  ////////////////////



  ////////////////////
  // BEGIN SET REFERENCE
  ////////////////////
  if(verbose_) {Rcpp::Rcout << "Setting reference ... ";};
  // set genotype specific parameters than all drives have
  reference::instance().set_reference(reproductionReference_["eta"], reproductionReference_["phi"],
                                      reproductionReference_["omega"], reproductionReference_["xiF"],
                                      reproductionReference_["xiM"], reproductionReference_["s"]);
  
  if(verbose_){Rcpp::Rcout <<  "\n\tBasics Done\n ... done setting reference!\n";};
  ////////////////////
  // END SET REFERENCE
  ////////////////////


  
  ////////////////////
  // BEGIN SET PARAMETERS
  ////////////////////
  if(verbose_) {Rcpp::Rcout << "Setting parameters ... ";};

  // setup input matrices and then fill them. They are std::vectors, not Rcpp things
  int dim = migrationMale_.ncol();
  Rcpp::NumericVector hold(dim), hold2(2);
  Rcpp::IntegerVector iHold(dim);
  dMat mHold(dim, dVec(dim)), fHold(dim, dVec(dim)), bsHold(dim, dVec(2)), bmHold(dim, dVec(dim)), sampCov(dim, dVec(5));
  iMat sampDay(dim, iVec(5)); // 5 because 5 life stages; E, L, P, Male, Female

  // loop over number of dims
  for(size_t i=0; i<dim; ++i){
    // male migration
    //  multiply move var by matrix now, since this is all we ever do with it.
    hold = networkParameters_["moveVar"] * migrationMale_.row(i);
    mHold[i] = Rcpp::as<dVec>(hold);

    // female migration
    hold = networkParameters_["moveVar"] * migrationFemale_.row(i);
    fHold[i] = Rcpp::as<dVec>(hold);
    
    // batch sex discrepancy
    hold2 = Rcpp::as<Rcpp::NumericMatrix>(migrationBatch_["sexProbs"]).row(i);
    bsHold[i]= Rcpp::as<dVec>(hold2);
    // batch migration
    hold = Rcpp::as<Rcpp::NumericMatrix>(migrationBatch_["moveProbs"]).row(i);
    bmHold[i] = Rcpp::as<dVec>(hold);
    
    // sampling parameters
    //  updated to new shapes so each patch can be different
    hold = Rcpp::as<Rcpp::NumericMatrix>(samplingParameters_["samplingCoverage"]).row(i);
    sampCov[i] = Rcpp::as<dVec>(hold);
    
    iHold = Rcpp::as<Rcpp::IntegerMatrix>(samplingParameters_["samplingDays"]).row(i);
    sampDay[i] = Rcpp::as<iVec>(iHold);
  }

  // now set parameters
  parameters::instance().set_parameters(networkParameters_["nPatch"],networkParameters_["simTime"],networkParameters_["runID"],
                                        networkParameters_["stageTime"],networkParameters_["beta"],networkParameters_["mu"],
                                        networkParameters_["alpha"],networkParameters_["Leq"],networkParameters_["AdPopEQ"],
                                        mHold, fHold, migrationBatch_["batchProbs"], bsHold, bmHold,
                                        sampDay, sampCov);

  if(verbose_){Rcpp::Rcout <<  " ... done setting parameters!\n";};

  ////////////////////
  // END SET PARAMETERS
  ////////////////////



  ////////////////////
  // BEGIN INITIALIZE NETWORK
  ////////////////////
  if(verbose_){Rcpp::Rcout << "Initializing network ... \n\tbigBrother is watching\n";};
  
  // vector of bigBrothers to get ID variables from
  std::vector<bigBrother> bigGov;
  bigGov.reserve(numThreads_);
  for(size_t i=0; i < numThreads_; i++){
    bigGov.emplace_back(bigBrother(i, numThreads_));
  }

  // vector of patches
  size_t numPatches = parameters::instance().get_n_patch();
  std::vector<std::unique_ptr<Patch> > patches;
  patches.reserve(numPatches);
  
  Rcpp::List patchRelease;

  // setup all patches
  for(size_t np=0; np<numPatches; ++np){

    // pull out  relaeses for this patch
    patchRelease = patchReleases_[np];

    // setup patch
    patches.emplace_back(std::make_unique<Family>(np,
                                                  patchRelease["maleReleases"],
                                                  patchRelease["femaleReleases"],
                                                  patchRelease["eggReleases"],
                                                  bigGov[np % numThreads_]));
  } // end network loop

  if(verbose_){Rcpp::Rcout <<  " ... done initializing network!\n";};
  ////////////////////
  // END INITIALIZE NETWORK
  ////////////////////





  // initialize things for inside the loop
  int tMax = parameters::instance().get_sim_time();

  // setup vectors of ofstreams
  std::vector<std::vector<std::ofstream *> > output(numPatches, std::vector<std::ofstream *>(5));

  // setup strings for file names
  std::vector<std::string> sHoldVec(5);
  std::string sHold;

  ////////////////////
  // BEGIN REPETITION WRAP
  ////////////////////
  for(size_t rep=0; rep<numReps_; ++rep){
    
    if(verbose_){Rcpp::Rcout <<  "Begin repetition " << rep+parameters::instance().get_run_id()<< "\n";};

    ////////////////////
    // BEGIN INITIALIZE LOGGER
    ////////////////////
    if(verbose_) {Rcpp::Rcout << "Initializing logging ... ";};

    dim = 0;
    for(const std::string& stage : {"/E_Run_","/L_Run_","/P_Run_","/M_Run_","/F_Run_"}){
      // base string name for each output
      sHoldVec[dim] = outputDirectory_ + stage
      + std::string(3 - std::to_string(rep+parameters::instance().get_run_id()).length(), '0')
      + std::to_string(rep+parameters::instance().get_run_id())
      + "_Patch_";

      // reusing dim from above, increment here
      dim++;
    } // end loop over file base names

    // open all streams with file names
    // try in parallel. It's independent, but idk if this is a speedup or just a waste
    // maybe schedule as static? 
    #pragma omp parallel for default(shared) private(sHold) schedule(auto)
    for(size_t np=0; np<numPatches; ++np){
      // denote patch
      sHold = std::string(3 - std::to_string(np).length(), '0')
      + std::to_string(np)
      + ".csv";

      // open streams
      for(size_t stage=0; stage < 5; ++stage){
        // initialize ofstream and open with file name
        output[np][stage] = new std::ofstream(sHoldVec[stage] + sHold);
      } // end loop over stages

    } // end loop over patches

    if(verbose_){Rcpp::Rcout <<  " ... done initializing logging!\n";};
    /////////////////////
    // END INITIALIZE LOGGER
    ////////////////////
    
    
    
    ////////////////////
    // BEGIN INITIALIZE OUTPUT
    ////////////////////
    if(verbose_){Rcpp::Rcout <<  "Initializing output ... ";};
    
    // output is to different files, so parallel
    #pragma omp parallel for default(shared) schedule(dynamic)
    for(size_t np=0; np<numPatches; ++np){
      patches[np]->init_output( output[patches[np]->get_patchID()] );
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
      #pragma omp parallel for default(shared) private(myThread) schedule(dynamic)
      for(size_t np=0; np<numPatches; ++np){
        #ifdef _OPENMP
        // get unique thread for prng
        myThread = omp_get_thread_num();
        #endif
        
        // run today's stuff
        patches[np]->oneDay_popDynamics(randInst[myThread], bigGov[myThread]);
      }

      // Do in-bound migration in one large loop, no extra structures
      for(size_t inPatch=0; inPatch < numPatches; ++inPatch){
        for(size_t outPatch=0; outPatch < numPatches; ++outPatch){
          patches[inPatch]->oneDay_migrationIn(patches[outPatch]->get_maleMigration(inPatch),
                                               patches[outPatch]->get_femaleMigration(inPatch));
        }
      }

      // Log output
      #pragma omp parallel for default(shared) private(myThread) schedule(dynamic)
      for(size_t np=0; np<numPatches; ++np){
        #ifdef _OPENMP
        // get unique thread for prng
        myThread = omp_get_thread_num();
        #endif
        
        // run output
        patches[np]->oneDay_writeOutput(output[patches[np]->get_patchID()], randInst[myThread]);
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
    #pragma omp parallel for default(shared) schedule(auto)
    for(size_t np=0; np<numPatches; ++np){
      for(size_t stage=0; stage < 5; ++stage){
        output[np][stage]->close();
      } // end stage loop
    } // end patch loop


    ////////////////////
    // RESET SIMULATION
    ////////////////////
    //  This does involve rebuilding the starting mosquitoes
    // skip if it's the last rep
    if((rep+1) == numReps_){
      if(verbose_){Rcpp::Rcout <<  "End repetition " << rep+parameters::instance().get_run_id() << "\n";};
      break;
    }

    parameters::instance().reset_t_now();
    
    #pragma omp parallel for default(shared) schedule(auto)
    for(size_t i=0; i < numThreads_; ++i){
      bigGov[i].reset();
    }
    
    #pragma omp parallel for default(shared) private(myThread) schedule(dynamic)
    for(size_t np=0; np<numPatches; ++np){
      #ifdef _OPENMP
      // get unique thread for prng
      myThread = omp_get_thread_num();
      #endif
      
      patches[np]->reset_Patch(bigGov[myThread]);
    }

    if(verbose_){Rcpp::Rcout <<  "End repetition " << rep+parameters::instance().get_run_id() << "\n\n";};
  } // end repetition loop
  
  ////////////////////
  // END REPETITION WRAP
  ////////////////////
  

  // close profiler
  #ifdef BASE_PROFILER_H_
    ProfilerStop();
  #endif
}

