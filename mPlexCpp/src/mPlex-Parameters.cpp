/*
 #      __  _____________       _       ______
 #     /  |/  / ____/ __ \_____(_)   __/ ____/
 #    / /|_/ / / __/ / / / ___/ / | / / __/
 #   / /  / / /_/ / /_/ / /  / /| |/ / /___
 #  /_/  /_/\____/_____/_/  /_/ |___/_____/
 #
 #  Marshall Lab
 #  January 2018
 #  Global Parameters Singleton
 #
*/

#include "mPlex-Parameters.hpp"

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
                                const int& n_patch_, const int& sim_time_, const double& move_var_, const int& run_id_,
                                /* biological parameters */
                                const std::vector<int>& stage_time_, const double& beta_, const std::vector<double>& mu_,
                                /* patch-specific derived parameters */
                                const std::vector<double>& alpha_, const std::vector<int>& larva_eq_, const std::vector<int>& adult_pop_eq_,
                                // batch parameters
                                const std::vector<double>& batchProbs_, const arma::Cube<double>& sexProbs_, const Rcpp::NumericMatrix& moveMat_){



  /* simulation fields */
  n_patch = n_patch_;
  sim_time = sim_time_;
  t_now = 0; /* starts at 0 because initial condition*/
  move_var = move_var_;
  run_id = run_id_;

  /* biological parameters */
  stage_time = stage_time_;
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
