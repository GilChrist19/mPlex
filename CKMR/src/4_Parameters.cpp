///////////////////////////////////////////////////////////////////////////////
//     ________ __ __  _______ 
//    / ____/ //_//  |/  / __ \
//   / /   / ,<  / /|_/ / /_/ /
//  / /___/ /| |/ /  / / _, _/ 
//  \____/_/ |_/_/  /_/_/ |_|  
//   
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
                                const int& n_patch_, const int& sim_time_,const int& run_id_,
                                /* biological parameters */
                                const iVec& stage_time_, const double& beta_, const bool& beta_const_,
                                const dVec& mu_, const int& male_max_age_, const int& female_max_age_,
                                /* patch-specific derived parameters */
                                const arma::Mat<double>& alpha_, const iVec& larva_eq_, const iVec& adult_pop_eq_,
                                // migration
                                const dMat& male_migration_, const dMat& female_migration_,
                                // batch parameters
                                const dVec& batchProbs_, const dMat& sexProbs_, const dMat& moveMat_,
                                // sampling parameters
                                const arma::Cube<unsigned int>& sampDays_, const arma::Cube<double>& sampCov_){



  /* simulation fields */
  n_patch = n_patch_;
  sim_time = sim_time_;
  t_now = 0; /* starts at 0 because initial condition*/
  male_migration = male_migration_;
  female_migration = female_migration_;
  run_id = run_id_;

  /* biological parameters */
  stage_time = stage_time_;
  stage_sum.resize(stage_time.size());
    std::partial_sum(stage_time.begin(), stage_time.end(), stage_sum.begin());
  beta = beta_;
  beta_const = beta_const_;
  mu = mu_;
  male_max_age = male_max_age_;
  female_max_age = female_max_age_;

  /* set patch-specific parameters */
  alpha = alpha_;
  larva_eq = larva_eq_;
  adult_pop_eq = adult_pop_eq_;
  
  /* set batch migration parameters */
  batchProbs = batchProbs_;
  sexProbs = sexProbs_;
  batchLocations = moveMat_;
  
  /* set sampling parameters */
  sampDays = sampDays_;
  sampCov = sampCov_;
  
};
