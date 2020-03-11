///////////////////////////////////////////////////////////////////////////////
//                           ____  __          ______          
//                ____ ___  / __ \/ /__  _  __/ ____/___  ____ 
//               / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
//              / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
//             /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/ 
//                                              /_/   /_/      
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
  
  void                set_mendelian(const Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericVector> >& probsF_,
                                    const Rcpp::ListOf<Rcpp::ListOf<Rcpp::StringVector> >& allelesF_,
                                    const Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericVector> >& probsM_,
                                    const Rcpp::ListOf<Rcpp::ListOf<Rcpp::StringVector> >& allelesM_);
  void                set_homing(const Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericVector> >& probsF_,
                                 const Rcpp::ListOf<Rcpp::ListOf<Rcpp::StringVector> >& allelesF_,
                                 const Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericVector> >& probsM_,
                                 const Rcpp::ListOf<Rcpp::ListOf<Rcpp::StringVector> >& allelesM_);
  void                set_cutting(const Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericVector> >& probsF_,
                                  const Rcpp::ListOf<Rcpp::ListOf<Rcpp::StringVector> >& allelesF_,
                                  const Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericVector> >& probsM_,
                                  const Rcpp::ListOf<Rcpp::ListOf<Rcpp::StringVector> >& allelesM_);
  void                set_alleleTypes(const Rcpp::ListOf<Rcpp::List>& alleleList_){alleleTypes = alleleList_;};

  
  // getters
  double      get_eta(const std::string& genType);
  double      get_phi(const std::string& genType);
  double      get_omega(const std::string& genType);
  double      get_xiF(const std::string& genType);
  double      get_xiM(const std::string& genType);
  double      get_s(const std::string& genType);

  dVec::iterator    get_mendelian_probs_begin(size_t sex, size_t locus, size_t allele);
  dVec::iterator    get_mendelian_probs_end(size_t sex, size_t locus, size_t allele);
  sVec::iterator    get_mendelian_allele_begin(size_t sex, size_t locus, size_t allele);
  sVec::iterator    get_mendelian_allele_end(size_t sex, size_t locus, size_t allele);
  
  dVec::iterator    get_homing_probs_begin(size_t sex, size_t locus, size_t allele);
  dVec::iterator    get_homing_probs_end(size_t sex, size_t locus, size_t allele);
  sVec::iterator    get_homing_allele_begin(size_t sex, size_t locus, size_t allele);
  sVec::iterator    get_homing_allele_end(size_t sex, size_t locus, size_t allele);
  
  dVec::iterator    get_cutting_probs_begin(size_t sex, size_t locus, size_t allele);
  dVec::iterator    get_cutting_probs_end(size_t sex, size_t locus, size_t allele);
  sVec::iterator    get_cutting_allele_begin(size_t sex, size_t locus, size_t allele);
  sVec::iterator    get_cutting_allele_end(size_t sex, size_t locus, size_t allele);
  
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
  
  std::vector<dArVec>                       mendelian_probs_f;
  std::vector<sArVec>                       mendelian_alleles_f;
  std::vector<dArVec>                       mendelian_probs_m;
  std::vector<sArVec>                       mendelian_alleles_m;

  std::vector<dArVec>                       homing_probs_f;
  std::vector<sArVec>                       homing_alleles_f;
  std::vector<dArVec>                       homing_probs_m;
  std::vector<sArVec>                       homing_alleles_m;

  std::vector<dArVec>                       cutting_probs_f;
  std::vector<sArVec>                       cutting_alleles_f;
  std::vector<dArVec>                       cutting_probs_m;
  std::vector<sArVec>                       cutting_alleles_m;

  
  Rcpp::ListOf<Rcpp::List>                  alleleTypes;

};

#endif
