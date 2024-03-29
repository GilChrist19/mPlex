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
#include "4_BigBrother.hpp"


/******************************************************************************
 * Run the simulation
******************************************************************************/
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
// [[Rcpp::export]]
void run_mPlex(const std::uint64_t& s1_,
               const std::uint64_t& s2_,
               const std::uint64_t& s3_,
               const std::uint64_t& s4_,
               const uint_least32_t& numReps_,
               const uint_least32_t& numThreads_,
               const Rcpp::List& networkParameters_,
               const Rcpp::List& reproductionReference_,
               const Rcpp::ListOf<Rcpp::List>& initAlleles_, 
               const Rcpp::List& patchReleases_,
               const Rcpp::NumericMatrix& migrationMale_,
               const Rcpp::NumericMatrix& migrationFemale_,
               const Rcpp::List& migrationBatch_,
               const std::string& reproductionType_,
               const std::string& outputDirectory_,
               const bool& verbose_){

  #ifdef BASE_PROFILER_H_
    // make sure to change path!!! But keep file name.
    ProfilerStart("/home/gilchrist/Desktop/OUTPUT/profile.log");
  #endif

  // define outside ifdef so that it works when _OPENMP isn't define
  int myThread(0);
  
  // set threads if using openMP
  #ifdef _OPENMP
    omp_set_num_threads(numThreads_);
  #else
    const_cast<uint_least32_t&>(numThreads_) = 1;
  #endif
    
    
    
  ////////////////////
  // BEGIN INITIALIZE PRNG
  ////////////////////
  if(verbose_) {Rcpp::Rcout << "Initializing " << numThreads_ << " prngs ... ";};

  std::vector<prng> randInst;
  randInst.reserve(numThreads_);
  std::array<std::uint64_t, 4> seed = {s1_,s2_,s3_,s4_};
  for(auto i=0; i < numThreads_; i++){
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
  if(verbose_) {Rcpp::Rcout << "\n\tBasics Done\n";};
  
  if(reproductionType_ != "Family"){
    // set allele types for initialization etc
    reference::instance().set_alleleTypes(initAlleles_);
    
    // set mendelian allele reference
    // reference::instance().set_mendelian(Rcpp::as<Rcpp::List> (reproductionReference_["mendelian"]),
    //                                     Rcpp::as<Rcpp::List> (reproductionReference_["mendelianAlleles"]));
    reference::instance().set_mendelian(Rcpp::as<Rcpp::ListOf<Rcpp::List> > (reproductionReference_["mendelian"])["female"],
                                        Rcpp::as<Rcpp::ListOf<Rcpp::List> > (reproductionReference_["mendelianAlleles"])["female"],
                                        Rcpp::as<Rcpp::ListOf<Rcpp::List> > (reproductionReference_["mendelian"])["male"],
                                        Rcpp::as<Rcpp::ListOf<Rcpp::List> > (reproductionReference_["mendelianAlleles"])["male"]);
    if(verbose_) {Rcpp::Rcout << "\tMendelian Done\n";};
    
    // set homing allele reference
    // reference::instance().set_homing(Rcpp::as<Rcpp::List> (reproductionReference_["homing"]),
    //                                   Rcpp::as<Rcpp::List> (reproductionReference_["homingAlleles"]));
    reference::instance().set_homing(Rcpp::as<Rcpp::ListOf<Rcpp::List> > (reproductionReference_["homing"])["female"],
                                     Rcpp::as<Rcpp::ListOf<Rcpp::List> > (reproductionReference_["homingAlleles"])["female"],
                                     Rcpp::as<Rcpp::ListOf<Rcpp::List> > (reproductionReference_["homing"])["male"],
                                     Rcpp::as<Rcpp::ListOf<Rcpp::List> > (reproductionReference_["homingAlleles"])["male"]);
    if(verbose_) {Rcpp::Rcout << "\tHoming Done\n";};
  }

  // set cutting reference if daisy drive
  if(reproductionType_ == "DaisyDrive"){
    // reference::instance().set_cutting(Rcpp::as<Rcpp::List> (reproductionReference_["cutting"]),
    //                                   Rcpp::as<Rcpp::List> (reproductionReference_["cuttingAlleles"]));
    reference::instance().set_cutting(Rcpp::as<Rcpp::ListOf<Rcpp::List> > (reproductionReference_["cutting"])["female"],
                                      Rcpp::as<Rcpp::ListOf<Rcpp::List> > (reproductionReference_["cuttingAlleles"])["female"],
                                      Rcpp::as<Rcpp::ListOf<Rcpp::List> > (reproductionReference_["cutting"])["male"],
                                      Rcpp::as<Rcpp::ListOf<Rcpp::List> > (reproductionReference_["cuttingAlleles"])["male"]);
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
  parameters::instance().set_parameters(networkParameters_["nPatch"],networkParameters_["simTime"],
                       networkParameters_["sampTime"],networkParameters_["runID"],
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
                                                   patchRelease["eggReleases"],
                                                   patchRelease["matedFemaleReleases"]));
      
    } else if(reproductionType_ == "mPlex_oLocus"){
      
      if(verbose_ && np==0){Rcpp::Rcout << "\tInitializing oneLocus Drive"<<std::endl;};
      
      patches.emplace_back(std::make_unique<oneLocus>(np,
                                                   patchRelease["maleReleases"],
                                                   patchRelease["femaleReleases"],
                                                   patchRelease["eggReleases"],
                                                   patchRelease["matedFemaleReleases"]));
      
    } else if(reproductionType_ == "mPlex_mLoci"){
      
      if(verbose_ && np==0){Rcpp::Rcout << "\tInitializing multiLocus Drive"<<std::endl;};
      
      patches.emplace_back(std::make_unique<multiLocus>(np,
                                                   patchRelease["maleReleases"],
                                                   patchRelease["femaleReleases"],
                                                   patchRelease["eggReleases"],
                                                   patchRelease["matedFemaleReleases"]));
      
    } else if(reproductionType_ == "Family"){
      
      if(verbose_ && np==0){Rcpp::Rcout << "\tBigBrother is watching"<<std::endl;};
      
      patches.emplace_back(std::make_unique<Family>(np,
                                                    patchRelease["maleReleases"],
                                                    patchRelease["femaleReleases"],
                                                    patchRelease["eggReleases"],
                                                    patchRelease["matedFemaleReleases"]));
      
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
    
    if(verbose_){Rcpp::Rcout <<  "Begin repetition " << rep+parameters::instance().get_run_id()<< "\n";};

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
    // maybe schedule as static? 
    #pragma omp parallel for default(shared) private(sHold) schedule(auto)
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
    //  patches could be different size, and all things are written, so putting 
    //  this with dynamic allocation, but not sure if that's beneficial
    #pragma omp parallel for default(shared) schedule(dynamic)
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
    Progress pb(tMax,verbose_);
    while(parameters::instance().get_t_now() <= tMax){

      // test for user interrupt
      if(checkInterrupt()) return;

      // Independent daily operations
      //  use dynamic schedule?
      //  Could be better once patches aren't all the same size
      #pragma omp parallel for default(shared) private(myThread) schedule(dynamic)
      for(size_t np=0; np<numPatches; ++np){
        #ifdef _OPENMP
        // get unique thread for prng
        myThread = omp_get_thread_num();
        #endif
        
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
      //  Same thoughts as for init_output above
      if( (parameters::instance().get_t_now() % parameters::instance().get_samp_time()) == 0 ){
        #pragma omp parallel for default(shared) schedule(dynamic)
        for(size_t np=0; np<numPatches; ++np){
          patches[np]->oneDay_writeOutput(M_output[patches[np]->get_patchID()], F_output[patches[np]->get_patchID()]);
        }
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
      M_output[np].close();
      F_output[np].close();
    }
    
    
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
    BigBrother::instance().reset();
    #pragma omp parallel for default(shared) schedule(dynamic)
    for(size_t np=0; np<numPatches; ++np){
      patches[np]->reset_Patch();
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



