///////////////////////////////////////////////////////////////////////////////
//                           ____  __          ______          
//                ____ ___  / __ \/ /__  _  __/ ____/___  ____ 
//               / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
//              / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
//             /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/ 
//                                              /_/   /_/      
///////////////////////////////////////////////////////////////////////////////

#include "4_PRNG.hpp"

///////////////////////////////////////////////////////////////////////////////
// constructor & destructor
///////////////////////////////////////////////////////////////////////////////
// prng::prng(const uint_least32_t& seed) : rng(seed){
//   runif = std::uniform_real_distribution<double>(0,1);
// };

prng::prng(const std::array<std::uint64_t, 4>& seed){
  // set seed
  rng.seed(seed);
  // initialize uniform
  runif = std::uniform_real_distribution<double>(0,1);
};
prng::~prng(){};


///////////////////////////////////////////////////////////////////////////////
// continuous random variate sampling
///////////////////////////////////////////////////////////////////////////////
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


///////////////////////////////////////////////////////////////////////////////
// set one sample probs
///////////////////////////////////////////////////////////////////////////////
void prng::set_oneSample(const std::vector<double>& prob){
  sample = std::discrete_distribution<size_t> (prob.begin(),prob.end());
}

size_t prng::get_oneSample(){
  return sample(rng);
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

// bernoulli things
bool prng::get_rBern(const double& p){
  std::bernoulli_distribution berni(p);
  return berni(rng);
}

void prng::set_cBern(const double& p){
  bernoulli = std::bernoulli_distribution(p);
};

bool prng::get_cBern(){
  return bernoulli(rng);
};

std::vector<int> prng::get_rmultinom(const int &size, const std::vector<double>& prob){
  std::vector<int>sample(prob.size(),0);
  
  // if there is only one thing to draw from, return it
  if(prob.size()==0){
    sample[0] = size;
    return sample;
  }
  
  // Draw from distribution
  //  Grows with pop size, need to fix? Does not grow with probs size
  std::discrete_distribution<int>rmultinom(prob.begin(),prob.end());
  for(size_t i=0; i<size; i++){
    size_t j = rmultinom(rng);
    sample[j] += 1;
  }
  return sample;
};

/* Startek, M. (2016). An asymptotically optimal, online algorithm for weighted random sampling with replacement, 1–11. Retrieved from http://arxiv.org/abs/1611.00532 */
std::vector<int> prng::get_rmultinom_online(int size, const std::vector<double>& prob,
                                            const double& switchover){

  std::vector<int> sample(prob.size(),0);

  double pprob = 0.0;
  double cprob = 0.0;
  unsigned int pidx(0);
  double p;
  int nrtaken(0);
    
  while(size > 0){
    
    pprob += prob[pidx];
    while(((pprob - cprob) * size / (1.0 - cprob)) < switchover){
      cprob += get_beta_1_b(size) * (1.0 - cprob);
      while(pprob < cprob)
        pprob += prob[++pidx];
      if(sample.size() == pidx)
        sample[pidx] = 1;
      else
        sample[pidx] += 1;
      size--;
      if(size == 0)
        break;
    } // end inner while
    
    if(size == 0) break;
    
    p = (pprob-cprob)/(1.0-cprob);
    nrtaken = 0;
    
    if(p >= 1.0){
      nrtaken = size;
    } else {
      std::binomial_distribution<int> rbinom(size, p);
      nrtaken = rbinom(rng);
    }
    
    if(nrtaken > 0){
      if(sample.size() == pidx)
        sample[pidx] = nrtaken;
      else
        sample[pidx] += nrtaken;
    }
    
    size -= nrtaken;
    pidx++;
    cprob = pprob;
  } // end outer while

  return sample;
};

std::vector<double> prng::get_rdirichlet(const std::vector<double>& prob){
  
  std::vector<double> sample(prob);
  double hold = 0.0;
  
  for(auto& sampleIt : sample){
    if(sampleIt==0) continue;
    std::gamma_distribution<double> gamma(sampleIt, 1.0);
    sampleIt = gamma(rng);
    // get total sum of draws
    hold += sampleIt;
  }
  
  // invert hold
  hold = 1.0/hold;

  // normalize
  for(auto& sampleIt : sample){
    sampleIt *= hold;
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
