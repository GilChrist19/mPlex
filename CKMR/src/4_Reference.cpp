///////////////////////////////////////////////////////////////////////////////
//     ________ __ __  _______ 
//    / ____/ //_//  |/  / __ \
//   / /   / ,<  / /|_/ / /_/ /
//  / /___/ /| |/ /  / / _, _/ 
//  \____/_/ |_/_/  /_/_/ |_|  
//   
///////////////////////////////////////////////////////////////////////////////

#include "4_Reference.hpp"

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

// set genotype dependent parameters
void reference::set_reference(const Rcpp::NumericVector& eta_, const Rcpp::NumericVector& phi_,
                              const Rcpp::NumericVector& omega_, const Rcpp::NumericVector& xiF_,
                              const Rcpp::NumericVector& xiM_, const Rcpp::NumericVector& s_){

  // names vector for all things
  sVec listNames;
  
  // set genotype specific mating fitness, eta
  if(eta_.size() !=0 ){
    listNames = Rcpp::as<sVec>(eta_.names());
    for(size_t it=0; it<eta_.size(); ++it){ // for all of them, fill the map, casting things as key,value pairs
      eta.insert(std::make_pair(listNames[it], eta_[it]) );
    }
  }

  // set genotype specific sex ratio at emergence, phi
  if(phi_.size() != 0){
    listNames = Rcpp::as<sVec>(phi_.names());
    for(size_t it=0; it<phi_.size(); ++it){
      phi.insert(std::make_pair(listNames[it], phi_[it]) );
    }
  }

  // set genotype specific multiplicative modifier of adult mortality, omegs
  if(omega_.size() != 0){
    listNames = Rcpp::as<sVec>(omega_.names());
    for(size_t it=0; it<omega_.size(); ++it){
      omega.insert(std::make_pair(listNames[it], omega_[it]) );
    }
  }
  
  // set genotype specific female pupatory success, xiF
  if(xiF_.size() != 0){
    listNames = Rcpp::as<sVec>(xiF_.names());
    for(size_t it=0; it<xiF_.size(); ++it){
      xiF.insert(std::make_pair(listNames[it], xiF_[it]) );
    }
  }
  
  //set genotype specific male pupatory success, xiM
  if(xiM_.size() != 0){
    listNames = Rcpp::as<sVec>(xiM_.names());
    for(size_t it=0; it<xiM_.size(); ++it){
      xiM.insert(std::make_pair(listNames[it], xiM_[it]) );
    }
  }
  
  // set fertility fraction, s
  if(s_.size() != 0){
    listNames = Rcpp::as<sVec>(s_.names());
    for(size_t it=0; it<s_.size(); ++it){
      s.insert(std::make_pair(listNames[it], s_[it]) );
    }
  }
  
};


// get genotype dependent parameters
double reference::get_eta(const std::string& genType){

  // iterator to element if it exists in the map
  itHold = eta.find(genType);
  
  // if it doesn't exist, it returns the end of the map, so check that
  if(itHold != eta.end()){
    return itHold->second;
  } else {
    return 1.0; // default value
  }

}

double reference::get_phi(const std::string& genType){
  
  // iterator to element if it exists in the map
  itHold = phi.find(genType);
  
  // if it doesn't exist, it returns the end of the map, so check that
  if(itHold != phi.end()){
    return itHold->second;
  } else {
    return 0.5; // default value
  }
  
}

double reference::get_omega(const std::string& genType){

  // iterator to element if it exists in the map
  itHold = omega.find(genType);
  
  // if it doesn't exist, it returns the end of the map, so check that
  if(itHold != omega.end()){
    return itHold->second;
  } else {
    return 1.0; // default value
  }
  
}

double reference::get_xiF(const std::string& genType){
  
  // iterator to element if it exists in the map
  itHold = xiF.find(genType);
  
  // if it doesn't exist, it returns the end of the map, so check that
  if(itHold != xiF.end()){
    return itHold->second;
  } else {
    return 1.0; // default value
  }
  
}

double reference::get_xiM(const std::string& genType){
  
  // iterator to element if it exists in the map
  itHold = xiM.find(genType);
  
  // if it doesn't exist, it returns the end of the map, so check that
  if(itHold != xiM.end()){
    return itHold->second;
  } else {
    return 1.0; // default value
  }
  
}

double reference::get_s(const std::string& genType){
  
  // iterator to element if it exists in the map
  itHold = s.find(genType);
  
  // if it doesn't exist, it returns the end of the map, so check that
  if(itHold != s.end()){
    return itHold->second;
  } else {
    return 1.0; // default value
  }
  
}

/**************************************
 * reproduction parameters
**************************************/
void reference::set_mendelian(const Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericVector> >& probs_,
                              const Rcpp::ListOf<Rcpp::ListOf<Rcpp::StringVector> >& alleles_){

  // set size of outer vector
  mendelian_probs.resize(probs_.size());
  mendelian_alleles.resize(alleles_.size());;
   
  // holder objects to make life easier in the loop
  std::vector<std::string> holdA;
  
  // loop over number of loci
  for(size_t locus=0; locus < probs_.size(); ++locus){
    
    // size inner vector
    mendelian_probs[locus].resize(4);
    mendelian_alleles[locus].resize(4);
    
    // loop over possible alleles at that locus
    for(size_t allele : {0,1,2,3}){
      // convert things in holders
      holdA = Rcpp::as<sVec>(alleles_[locus][allele]);
      
      // set things in final reference
      mendelian_probs[locus][allele].insert(mendelian_probs[locus][allele].end(),
                                            probs_[locus][allele].begin(),
                                            probs_[locus][allele].end());

      mendelian_alleles[locus][allele].insert(mendelian_alleles[locus][allele].end(),
                                              holdA.begin(),
                                              holdA.end());
    }
  }

}
 
void reference::set_homing(const Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericVector> >& probs_,
                            const Rcpp::ListOf<Rcpp::ListOf<Rcpp::StringVector> >& alleles_){
  
  // set size of outer vector
  homing_probs.resize(probs_.size());
  homing_alleles.resize(alleles_.size());;
  
  // holder objects to make life easier in the loop
  std::vector<std::string> holdA;
  
  
   // loop over number of loci
   for(size_t locus=0; locus < probs_.size(); ++locus){
     
     // size inner vector
     homing_probs[locus].resize(4);
     homing_alleles[locus].resize(4);
     
     // loop over possible alleles at that locus
     for(size_t allele : {0,1,2,3}){
       // convert things in holders
       holdA = Rcpp::as<sVec>(alleles_[locus][allele]);
       
       // set things in final reference
       homing_probs[locus][allele].insert(homing_probs[locus][allele].end(),
                                             probs_[locus][allele].begin(),
                                             probs_[locus][allele].end());
       
       homing_alleles[locus][allele].insert(homing_alleles[locus][allele].end(),
                                               holdA.begin(),
                                               holdA.end());
     }
   }
   
}

void reference::set_cutting(const Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericVector> >& probs_,
                            const Rcpp::ListOf<Rcpp::ListOf<Rcpp::StringVector> >& alleles_){
  
  // set size of outer vector
  cutting_probs.resize(probs_.size());
  cutting_alleles.resize(alleles_.size());;
  
  // holder objects to make life easier in the loop
  std::vector<std::string> holdA;
  
  // loop over number of loci
  for(size_t locus=0; locus < probs_.size(); ++locus){
    
    // size inner vector
    cutting_probs[locus].resize(4);
    cutting_alleles[locus].resize(4);
    
    // loop over possible alleles at that locus
    for(size_t allele : {0,1,2,3}){
      // convert things in holders
      holdA = Rcpp::as<sVec>(alleles_[locus][allele]);
      
      // set things in final reference
      cutting_probs[locus][allele].insert(cutting_probs[locus][allele].end(),
                                         probs_[locus][allele].begin(),
                                         probs_[locus][allele].end());
      
      cutting_alleles[locus][allele].insert(cutting_alleles[locus][allele].end(),
                                           holdA.begin(),
                                           holdA.end());
    }
  }
  
}
