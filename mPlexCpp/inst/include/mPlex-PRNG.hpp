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

#ifndef MPLEX_PRNG_
#define MPLEX_PRNG_

#include <random>
#include <vector>
#include <iostream>

/* threadsafe prng singleton */
class prng final {
public:
    /* utility methods */
    static prng*                           instance(); /* get instance */
    void                                   set_seed(const uint_least32_t& seed);
    void                                   suicide();

    /* continuous random variate sampling */
    double                                 get_runif();
    double                                 get_rexp(const double& rate);
    double                                 get_rlnorm(const double& meanlog, const double& sdlog);

    /* discrete random variate sampling */
    int                                    get_rpois(const double& lambda);
    int                                    get_rbinom(const int& n, const double& p);
    std::vector<int>                       get_rmultinom(const int& size, const std::vector<double> prob);

    /* resample template type x 'size' times */
    template<typename T>
    T                                      get_resample(const T& x, const int& size);

private:
  /* constructor & destructor */
  prng();
  ~prng();

  /* delete all copy & move semantics */
  prng(const prng&) = delete;
  prng& operator=(const prng&) = delete;
  prng(prng&&) = delete;
  prng& operator=(prng&&) = delete;

  std::mt19937                            rng;
  std::uniform_real_distribution<double>  runif;
};

#endif
