///////////////////////////////////////////////////////////////////////////////
//                           ____  __          ______          
//                ____ ___  / __ \/ /__  _  __/ ____/___  ____ 
//               / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
//              / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
//             /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/ 
//                                              /_/   /_/      
///////////////////////////////////////////////////////////////////////////////

#include "4_PRNG.hpp"

/* constructor & destructor */
prng::prng(){};
prng::~prng(){};

/* utility methods */
prng& prng::instance(){
    static prng instance;
    return instance;
};

void prng::set_seed(const uint_least32_t &seed){
  rng.seed(seed);
  runif = std::uniform_real_distribution<double>(0,1);
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

size_t prng::get_oneSample(const std::vector<double>& prob){
  std::discrete_distribution<size_t> sample(prob.begin(),prob.end());
  return sample(rng);
}

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

std::vector<double> prng::get_rdirichlet(const std::vector<double>& prob){
  
  std::vector<double> sample(prob);
  
  for(auto& sampleIt : sample){
    std::gamma_distribution<double> gamma(sampleIt, 1.0);
    sampleIt = gamma(rng);
  }
  
  double hold = std::accumulate(sample.begin(), sample.end(), 0.0);
  
  for(auto& sampleIt : sample){
    sampleIt /= hold;
  }

  return sample;
}

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
