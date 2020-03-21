///////////////////////////////////////////////////////////////////////////////
//     ________ __ __  _______ 
//    / ____/ //_//  |/  / __ \
//   / /   / ,<  / /|_/ / /_/ /
//  / /___/ /| |/ /  / / _, _/ 
//  \____/_/ |_/_/  /_/_/ |_|  
//   
///////////////////////////////////////////////////////////////////////////////

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
  
  void                set_reference(const Rcpp::NumericVector& eta_, const Rcpp::NumericVector& phi_,
                                    const Rcpp::NumericVector& omega_, const Rcpp::NumericVector& xiF_,
                                    const Rcpp::NumericVector& xiM_, const Rcpp::NumericVector& s_);
  
  void                set_mendelian(const Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericVector> >& probs_,
                                    const Rcpp::ListOf<Rcpp::ListOf<Rcpp::StringVector> >& alleles_);
  void                set_homing(const Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericVector> >& probs_,
                                 const Rcpp::ListOf<Rcpp::ListOf<Rcpp::StringVector> >& alleles_);
  void                set_cutting(const Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericVector> >& probs_,
                                  const Rcpp::ListOf<Rcpp::ListOf<Rcpp::StringVector> >& alleles_);
  void                set_alleleTypes(const Rcpp::ListOf<Rcpp::List>& alleleList_){alleleTypes = alleleList_;};

  
  // getters
  double      get_eta(const std::string& genType);
  double      get_phi(const std::string& genType);
  double      get_omega(const std::string& genType);
  double      get_xiF(const std::string& genType);
  double      get_xiM(const std::string& genType);
  double      get_s(const std::string& genType);

  dVec::iterator    get_mendelian_probs_begin(size_t locus, size_t allele){return mendelian_probs[locus][allele].begin();};
  dVec::iterator    get_mendelian_probs_end(size_t locus, size_t allele){return mendelian_probs[locus][allele].end();};
  sVec::iterator    get_mendelian_allele_begin(size_t locus, size_t allele){return mendelian_alleles[locus][allele].begin();};
  sVec::iterator    get_mendelian_allele_end(size_t locus, size_t allele){return mendelian_alleles[locus][allele].end();};
  
  dVec::iterator    get_homing_probs_begin(size_t locus, size_t allele){return homing_probs[locus][allele].begin();};
  dVec::iterator    get_homing_probs_end(size_t locus, size_t allele){return homing_probs[locus][allele].end();};
  sVec::iterator    get_homing_allele_begin(size_t locus, size_t allele){return homing_alleles[locus][allele].begin();};
  sVec::iterator    get_homing_allele_end(size_t locus, size_t allele){return homing_alleles[locus][allele].end();};
  
  dVec::iterator    get_cutting_probs_begin(size_t locus, size_t allele){return cutting_probs[locus][allele].begin();};
  dVec::iterator    get_cutting_probs_end(size_t locus, size_t allele){return cutting_probs[locus][allele].end();};
  sVec::iterator    get_cutting_allele_begin(size_t locus, size_t allele){return cutting_alleles[locus][allele].begin();};
  sVec::iterator    get_cutting_allele_end(size_t locus, size_t allele){return cutting_alleles[locus][allele].end();};
  
  Rcpp::List        get_alleloTypes(size_t patch){return alleleTypes[patch];};
  
  
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
  std::unordered_map<std::string, double>   eta;
  std::unordered_map<std::string, double>   phi;
  std::unordered_map<std::string, double>   omega;
  std::unordered_map<std::string, double>   xiF;
  std::unordered_map<std::string, double>   xiM;
  std::unordered_map<std::string, double>   s;
  std::unordered_map<std::string, double>::iterator itHold;
  
  std::vector<dArVec>                       mendelian_probs;
  std::vector<sArVec>                       mendelian_alleles;
  std::vector<dArVec>                       homing_probs;
  std::vector<sArVec>                       homing_alleles;
  std::vector<dArVec>                       cutting_probs;
  std::vector<sArVec>                       cutting_alleles;
  
  Rcpp::ListOf<Rcpp::List>                  alleleTypes;

};

#endif
