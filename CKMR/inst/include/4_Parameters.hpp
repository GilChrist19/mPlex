///////////////////////////////////////////////////////////////////////////////
//     ________ __ __  _______ 
//    / ____/ //_//  |/  / __ \
//   / /   / ,<  / /|_/ / /_/ /
//  / /___/ /| |/ /  / / _, _/ 
//  \____/_/ |_/_/  /_/_/ |_|  
//   
///////////////////////////////////////////////////////////////////////////////

#ifndef PARAMETERS_CKMR
#define PARAMETERS_CKMR


#include <tuple>
#include <vector>
#include <numeric>
#include <RcppArmadillo.h>

using dMat = std::vector<std::vector<double> >;
using dVec = std::vector<double>;
using iMat = std::vector<std::vector<int> >;
using iVec = std::vector<int>;


class parameters final {
public:
  /* utility methods */
  static parameters&    instance();
  void                  set_parameters(/* simulation fields */
                                       const int& n_patch_, const int& sim_time_, const int& run_id_,
                                       /* biological parameters */
                                        const iVec& stage_time_, const double& beta_, const bool& beta_const_,
                                        const dVec& mu_, const int& male_max_age_, const int& female_max_age_,
                                       /* patch-specific derived parameters */
                                       const arma::Mat<double>& alpha_, const iVec& larva_eq_, const iVec& adult_pop_eq_,
                                       // migration
                                       const dMat& male_migration_, const dMat& female_migration_,
                                       // batch parameters
                                       const std::vector<double>& batchProbs_, const dMat& sexProbs_, const dMat& moveMat_,
                                       // sampling parameters
                                       const arma::Cube<unsigned int>& sampDays_, const arma::Cube<double>& sampCov_);
  
  
  /* accessors */

  /* simulation fields */
  int               get_n_patch(){return n_patch;};
  int               get_sim_time(){return sim_time;};
  int               get_t_now(){return t_now;};
  void              reset_t_now(){t_now = 0;};
  void              increment_t_now(){t_now++;};
  dVec              get_male_migration(const size_t& patch){return male_migration[patch];};
  dVec              get_female_migration(const size_t& patch){return female_migration[patch];};
  int               get_run_id(){return run_id;};


  /* biological parameters */
  int               get_stage_time(const size_t& stage){return stage_time[stage];};
  int               get_stage_sum(const size_t& stage){return stage_sum[stage];};
  double            get_beta(){return beta;};
  bool              get_beta_const(){return beta_const;};
  double            get_mu(const size_t& stage){return mu[stage];};
  int               get_male_max_age(){return male_max_age;};
  int               get_female_max_age(){return female_max_age;};
  
  /* patch-specific derived parameters */
  double            get_alpha(const size_t& patch){return alpha.at(t_now,patch);};
  int               get_larva_eq(const size_t& patch){return larva_eq[patch];};
  int               get_adult_pop_eq(const size_t& patch){return adult_pop_eq[patch];};
  
  
  /* batch migration parameters */
  double            get_batchProbs(const size_t& patch){return batchProbs[patch];};
  double            get_batchMale(const size_t& patch){return sexProbs[patch][0];};
  double            get_batchFemale(const size_t& patch){return sexProbs[patch][1];};
  dVec              get_batchLocation(const size_t& patch){return batchLocations[patch];};
  
  /* sampling parameters*/
  bool              get_sampDay(const size_t& stage, const size_t& patch){return sampDays.at(t_now,stage,patch);};
  double            get_sampCov(const size_t& stage, const size_t& patch){return sampCov.at(t_now,stage,patch);};
  
  
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
  int       n_patch;
  int       sim_time;
  int       t_now;
  int       run_id;
  dMat      male_migration;
  dMat      female_migration;
  
  /* biological parameters */
  iVec      stage_time;
  iVec      stage_sum;
  double    beta;
  bool      beta_const;
  dVec      mu;
  int       male_max_age;
  int       female_max_age;
  
  /* patch-specific derived parameters */
  arma::Mat<double>      alpha;
  iVec                   larva_eq;
  iVec                   adult_pop_eq;

  /* batch migration parameters */
  dVec      batchProbs;
  dMat      sexProbs;
  dMat      batchLocations;
  
  /* sampling parameters */
  arma::Cube<unsigned int>      sampDays;
  arma::Cube<double>            sampCov;
  
};

#endif
