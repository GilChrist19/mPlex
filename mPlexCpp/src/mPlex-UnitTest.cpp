/*                  ____  __
 *       ____ ___  / __ \/ /__  _  __
 *      / __ `__ \/ /_/ / / _ \| |/_/
 *     / / / / / / ____/ /  __/>  <
 *    /_/ /_/ /_/_/   /_/\___/_/|_|
 *
 *    Marshall Lab
 *    Unit Tests
 *    January 2018
*/

#include <Rcpp.h>

/* library source */
#include <iostream>

/* mPlex source */
#include "mPlex-PRNG.hpp"

//' @export
// [[Rcpp::export]]
void unit_test_PRNG(const uint_least32_t& seed){
  std::cout << "test threaded generation of random variates" << std::endl;
  prng::instance()->set_seed(seed); /* set seed of prng */

  // const size_t nthread = std::thread::hardware_concurrency();
  // std::cout<<"parallel ("<<nthread<<" threads):"<<std::endl;

  prng::instance()->suicide();
};
