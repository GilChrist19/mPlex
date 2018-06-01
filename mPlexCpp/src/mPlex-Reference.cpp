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

#include "mPlex-Reference.hpp"

/* constructor & destructor */
reference::reference(){};
reference::~reference(){};

/* utility methods */
reference& reference::instance(){
    static reference instance;
    return instance;
};

/**************************************
 * genotype dependent parameters
**************************************/

// s parameter
void reference::set_reference(const Rcpp::NumericVector& s_){

  // get genotypes from s vector
  sVec listNames = s_.names();
  
  // for all of them, fill the map, casting things as key,value pairs
  for(size_t it=0; it<s_.size(); ++it){
    s.insert(std::make_pair<std::string,double>(listNames[it], s_[it]) );
  }
  
};

double reference::get_s(std::string genType){
  
  double hold = 1.0;
  
  // iterator to element if it exists in the map
  std::unordered_map<std::string, double>::iterator it = s.find(genType);
  
  // if it doesn't exist, it returns the end of the map, so check that
  if(it != s.end()){
    hold = it->second;
  }
  
  // return
  return hold;
  
}


/**************************************
 * reproduction parameters
**************************************/

void reference::set_mendelian(const Rcpp::ListOf<Rcpp::ListOf<dVec> >& probs_,
                              const Rcpp::ListOf<Rcpp::ListOf<sVec> >& alleles_){
  
  // There are always 4 alleles possible
  std::vector<int> alleles = {0,1,2,3};
  
  // loop over number of loci
  for(size_t locus=0; locus < probs_.size(); ++locus){
    // loop over possible alleles at that locus
    for(size_t allele : alleles){
      mendelian_probs[locus][allele] = probs_[locus][allele];
      mendelian_alleles[locus][allele] = alleles_[locus][allele];
    }
  }

}
 
 
 void reference::set_homing(const Rcpp::ListOf<Rcpp::ListOf<dVec> >& probs_,
                            const Rcpp::ListOf<Rcpp::ListOf<sVec> >& alleles_){
   
   // There are always 4 alleles possible
   std::vector<int> alleles = {0,1,2,3};
   
   // loop over number of loci
   for(size_t locus=0; locus < probs_.size(); ++locus){
     // loop over possible alleles at that locus
     for(size_t allele : alleles){
       homing_probs[locus][allele] = probs_[locus][allele];
       homing_alleles[locus][allele] = alleles_[locus][allele];
     }
   }
   
}



void reference::set_cutting(const Rcpp::ListOf<Rcpp::ListOf<dVec> >& probs_,
                            const Rcpp::ListOf<Rcpp::ListOf<sVec> >& alleles_){
  
  // There are always 4 alleles possible
  std::vector<int> alleles = {0,1,2,3};
  
  // loop over number of loci
  for(size_t locus=0; locus < probs_.size(); ++locus){
    // loop over possible alleles at that locus
    for(size_t allele : alleles){
      cutting_probs[locus][allele] = probs_[locus][allele];
      cutting_alleles[locus][allele] = alleles_[locus][allele];
    }
  }
  
}
