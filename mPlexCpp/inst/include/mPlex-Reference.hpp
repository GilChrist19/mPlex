/*
 #      __  _____________       _       ______
 #     /  |/  / ____/ __ \_____(_)   __/ ____/
 #    / /|_/ / / __/ / / / ___/ / | / / __/
 #   / /  / / /_/ / /_/ / /  / /| |/ / /___
 #  /_/  /_/\____/_____/_/  /_/ |___/_____/
 #
 #  Marshall Lab
 #  January 2018
 #  Cube Singleton
 #
*/

#ifndef REFERENCE_MPLEX
#define REFERENCE_MPLEX

#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>

#include <Rcpp.h>


using dVec = std::vector<double>; 
using sVec = std::vector<std::string>;

using dArVec = std::vector<dVec>;
using sArVec = std::vector<sVec>;


class reference final {
public:


  /* utility methods */
  static reference&   instance();
  
  void                set_reference(const Rcpp::ListOf<double>& s_);
  
  void                set_mendelian(const Rcpp::ListOf<Rcpp::ListOf<dVec> >& probs_, const Rcpp::ListOf<Rcpp::ListOf<sVec> >& alleles_);
  void                set_homing(const Rcpp::ListOf<Rcpp::ListOf<dVec> >& probs_, const Rcpp::ListOf<Rcpp::ListOf<sVec> >& alleles_);
  void                set_cutting(const Rcpp::ListOf<Rcpp::ListOf<dVec> >& probs_, const Rcpp::ListOf<Rcpp::ListOf<sVec> >& alleles_);
  

  
  // HOW DO i SET AND GET S/??????
  // getters
  double      get_s(std::string genType);
  
  dVec        get_mendelian_probs(size_t locus, size_t allele){return mendelian_probs[locus][allele];};
  sVec        get_mendelian_allele(size_t locus, size_t allele){return mendelian_alleles[locus][allele];};
  dVec        get_homing_probs(size_t locus, size_t allele){return homing_probs[locus][allele];};
  sVec        get_homing_allele(size_t locus, size_t allele){return homing_alleles[locus][allele];};
  dVec        get_cutting_probs(size_t locus, size_t allele){return cutting_probs[locus][allele];};
  sVec        get_cutting_allele(size_t locus, size_t allele){return cutting_alleles[locus][allele];};
  
  
  
  
private:
  /* constructor & destructor */
  reference();
  ~reference();

  /* delete all copy & move semantics */
  reference(const reference&) = delete;
  reference& operator=(const reference&) = delete;
  reference(reference&&) = delete;
  reference& operator=(reference&&) = delete;

  /* fields */
  std::unordered_map<std::string, double>   s;
  
  std::vector<dArVec>                       mendelian_probs;
  std::vector<sArVec>                       mendelian_alleles;
  std::vector<dArVec>                       homing_probs;
  std::vector<sArVec>                       homing_alleles;
  std::vector<dArVec>                       cutting_probs;
  std::vector<sArVec>                       cutting_alleles;

};


#endif
