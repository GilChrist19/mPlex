/*                  ____  __
 *       ____ ___  / __ \/ /__  _  __
 *      / __ `__ \/ /_/ / / _ \| |/_/
 *     / / / / / / ____/ /  __/>  <
 *    /_/ /_/ /_/_/   /_/\___/_/|_|
 *
 *    Marshall Lab
 *    PRNG threadsafe singleton
 *    January 2018
*/

#include "mPlex-PRNG.hpp"

/* constructor & destructor */
prng::prng(){
  #ifdef DEBUG_MPLEX
  std::cout << "prng being born at " << this << std::endl;
  #endif
};

prng::~prng(){
  #ifdef DEBUG_MPLEX
  std::cout << "prng being killed at " << this << std::endl;
  #endif
};

/* utility methods */
prng* prng::instance(){
    static prng instance;
    return &instance;
};

void prng::set_seed(const uint_least32_t &seed){
  rng.seed(seed);
  runif = std::uniform_real_distribution<double>(0,1);
};

void prng::suicide(){
  prng::~prng();
};

/* continuous random variate sampling */
double prng::get_runif(){
  return runif(rng);
};

double prng::get_rexp(const double &rate){
  std::exponential_distribution<double>rexp(rate);
  return rexp(rng);
};

double prng::get_rlnorm(const double &meanlog, const double &sdlog){
  std::lognormal_distribution<double>rlnorm(meanlog,sdlog);
  return rlnorm(rng);
};

/* discrete random variate sampling */
int prng::get_rpois(const double &lambda){
  std::poisson_distribution<int>rpois(lambda);
  return rpois(rng);
};

int prng::get_rbinom(const int& n, const double& p){
  std::binomial_distribution<int>rbinom(n,p);
  return rbinom(rng);
};

std::vector<int> prng::get_rmultinom(const int &size, const std::vector<double> prob){
  std::vector<int>sample(prob.size(),0);
  std::discrete_distribution<int>rmultinom(prob.begin(),prob.end());
  for(size_t i=0; i<size; i++){
    size_t j = rmultinom(rng);
    sample[j] += 1;
  }
  return sample;
};

/* resample template type T x 'size' times */
template<typename T>
T prng::get_resample(const T &x, const int &size){
  std::uniform_int_distribution<int>runif_int(0,x.size()-1); /* guaranteed unbiased; no modulo arithmetic bias */
  T out;
  for(size_t i=0; i<size; i++){
    size_t j = runif_int(rng);
    out.push_back(x[j]); /* will use type T's copy constructor, so if using classes with move semantics, will need a new function for them */
  };
  return out;
};
