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

#ifndef PARAMETERS_MPLEX
#define PARAMETERS_MPLEX

#include <tuple>
#include <vector>

#include <RcppArmadillo.h>

class parameters final {
public:
  /* utility methods */
  static parameters&    instance();
  void                  set_parameters(/* simulation fields */
                                       const int& n_patch_, const int& sim_time_, const double& move_var_, const int run_id_,
                                       /* biological parameters */
                                        const std::vector<int>& stage_time_, const double& beta_, const std::vector<double>& mu_,
                                       /* patch-specific derived parameters */
                                       const std::vector<double>& alpha_, const std::vector<double>& larva_eq_, const std::vector<double>& adult_pop_eq_,
                                       // batch parameters
                                       const std::vector<double>& batchProbs_, const arma::Cube<double>& sexProbs_, const Rcpp::NumericMatrix& moveMat_);
  
  void                  set_batch_parameters(const std::vector<double>& batchProbs_, const arma::Cube<double>& sexProbs_, const Rcpp::NumericMatrix& moveMat_);

  
  
  
  
  
  
  /* accessors */

  /* simulation fields */
  int               get_n_patch(){return n_patch;};
  int               get_sim_time(){return sim_time;};
  int               get_t_now(){return t_now;};
  void              reset_t_now(){t_now = 0;};
  void              increment_t_now(){t_now++;};
  double            get_move_var(){return move_var;};
  int               get_run_id(){return run_id;};


  /* biological parameters */
  int               get_stage_time(size_t stage){return stage_time[stage];};
  int               get_stage_sum(size_t stage){return std::accumulate(stage_time.begin(), stage_time[stage], 0)}
  double            get_beta(){return beta;};
  int               get_mu(size_t stage){return mu[stage];};
  
  
  /* patch-specific derived parameters */
  double            get_alpha(const size_t patch){return alpha[patch];};
  int               get_larva_eq(const size_t patch){return larva_eq[patch];};
  int               get_adult_pop_eq(const size_t patch){return adult_pop_eq[patch];};
  
  

  
  
  
  /* batch migration parameters */
  double                get_batchProbs(const size_t patch){return batchProbs[patch];};
  std::vector<double>   get_batchMale(const size_t patch){return sexProbs.tube(patch, 0);};
  std::vector<double>   get_batchFemale(const size_t patch){return sexProbs.tube(patch,1);};
  std::vector<double>   get_batchLocation(const size_t patch){return batchLocations.row(patch);};
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
private:

  /* constructor & destructor */
  parameters();
  ~parameters();

  /* delete all copy & move semantics */
  parameters(const parameters&) = delete;
  parameters& operator=(const parameters&) = delete;
  parameters(parameters&&) = delete;
  parameters& operator=(parameters&&) = delete;

  /* fields */

  /* simulation fields */
  int                   n_patch;
  int                   sim_time;
  int                   t_now;
  double                move_var;
  int                   run_id;
  
  /* biological parameters */
  std::vector<int>      stage_time;
  double                beta;
  std::vector<double>   mu;
  
  /* patch-specific derived parameters */
  std::vector<double>   alpha;
  std::vector<int>      larva_eq;
  std::vector<int>      adult_pop_eq;

  /* batch migration parameters */
  std::vector<double>   batchProbs;
  arma::Cube<double>    sexProbs;
  Rcpp::NumericMatrix   batchLocations;

};

#endif
