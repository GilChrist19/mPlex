///////////////////////////////////////////////////////////////////////////////
//                           ____  __          ______          
//                ____ ___  / __ \/ /__  _  __/ ____/___  ____ 
//               / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
//              / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
//             /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/ 
//                                              /_/   /_/      
///////////////////////////////////////////////////////////////////////////////

#ifndef PRNG_MPLEX
#define PRNG_MPLEX

///////////////////////////////////////////////////////////////////////////////
// includes and forward declarations
///////////////////////////////////////////////////////////////////////////////

#include <random>
#include <vector>
#include <math.h>

#include "4_xoshiro.hpp"

///////////////////////////////////////////////////////////////////////////////
// class declaration
///////////////////////////////////////////////////////////////////////////////
/* threadsafe prng if one is created for each thread!! */
class prng {
public:

    /* constructor & destructor */
    prng(const std::array<std::uint64_t, 4>& seed);
    virtual ~prng();
//    ~prng() = default;

    /* delete all copy & move semantics */
//    prng(const prng&) = delete;
//    prng& operator=(const prng&) = delete;
//    prng(prng&&) = delete;
//    prng& operator=(prng&&) = delete;


    /* continuous random variate sampling */
    double                                 get_runif();
    double                                 get_rexp(const double& rate);
    double                                 get_rlnorm(const double& meanlog, const double& sdlog);
    double                                 get_beta_1_b(const double& b){return 1.0 - pow(runif(rng), 1.0/b);};
    std::vector<double>                    get_rdirichlet(const std::vector<double>& prob);

    /* discrete random variate sampling */
    int                                    get_rpois(const double& lambda);
    int                                    get_rbinom(const int& n, const double& p);
    bool                                   get_rBern(const double& p);
    void                                   set_cBern(const double& p);
    bool                                   get_cBern();
    void                                   set_oneSample(const std::vector<double>& prob);
    size_t                                 get_oneSample();
    std::vector<int>                       get_rmultinom(const int& size, const std::vector<double>& prob);
    std::vector<int>                       get_rmultinom_online(int size, const std::vector<double>& prob, const double& switchover = 1.0);
    
    /* resample template type T x 'size' times */
    template<typename T>
    T                                      get_resample(const T& x, const int& size);

private:

  //std::mt19937                            rng;
  xoshiro256ss				                    rng;
  std::uniform_real_distribution<double>  runif;
  std::bernoulli_distribution             bernoulli;
  std::discrete_distribution<size_t>      sample;

};

#endif
