///////////////////////////////////////////////////////////////////////////////
//                           ____  __          ______          
//                ____ ___  / __ \/ /__  _  __/ ____/___  ____ 
//               / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
//              / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
//             /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/ 
//                                              /_/   /_/      
///////////////////////////////////////////////////////////////////////////////

#include "4_Parameters.hpp"

/* constructor & destructor */
parameters::parameters(){};
parameters::~parameters(){};

/* utility methods */
parameters& parameters::instance(){
    static parameters instance;
    return instance;
};

/* set parameters: necessarily ugly, cannot use initializer lists because its not a constructor */
void parameters::set_parameters(/* simulation fields */
                                const int& n_patch_, const int& sim_time_,
                                const int& samp_time_, const int& run_id_,
                                /* biological parameters */
                                const std::vector<int>& stage_time_, const double& beta_, const dVec& mu_,
                                /* patch-specific derived parameters */
                                const dVec& alpha_, const iVec& larva_eq_, const iVec& adult_pop_eq_,
                                // migration
                                const dMat& male_migration_, const dMat& female_migration_,
                                // batch parameters
                                const dVec& batchProbs_, const dMat& sexProbs_, const dMat& moveMat_){



  /* simulation fields */
  n_patch = n_patch_;
  sim_time = sim_time_;
  samp_time = samp_time_;
  t_now = 0; /* starts at 0 because initial condition*/
  male_migration = male_migration_;
  female_migration = female_migration_;
  run_id = run_id_;

  /* biological parameters */
  stage_time = stage_time_;
  stage_sum.resize(stage_time.size());
    std::partial_sum(stage_time.begin(), stage_time.end(), stage_sum.begin());
  beta = beta_;
  mu = mu_;

  /* set patch-specific parameters */
  alpha = alpha_;
  larva_eq = larva_eq_;
  adult_pop_eq = adult_pop_eq_;
  
  /* set batch migration parameters */
  batchProbs = batchProbs_;
  sexProbs = sexProbs_;
  batchLocations = moveMat_;
  
};
