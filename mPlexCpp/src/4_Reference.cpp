///////////////////////////////////////////////////////////////////////////////
//                           ____  __          ______          
//                ____ ___  / __ \/ /__  _  __/ ____/___  ____ 
//               / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
//              / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
//             /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/ 
//                                              /_/   /_/      
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
// set mendelian inheritance
void reference::set_mendelian(const Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericVector> >& probsF_,
                              const Rcpp::ListOf<Rcpp::ListOf<Rcpp::StringVector> >& allelesF_,
                              const Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericVector> >& probsM_,
                              const Rcpp::ListOf<Rcpp::ListOf<Rcpp::StringVector> >& allelesM_){

  // set size of outer vector
  //  All should be the same size
  mendelian_probs_f.resize(probsF_.size());
  mendelian_alleles_f.resize(allelesF_.size());
  mendelian_probs_m.resize(probsF_.size());
  mendelian_alleles_m.resize(allelesF_.size());
   
  // holder objects to make life easier in the loop
  std::vector<std::string> holdA;
  
  // loop over number of loci
  for(size_t locus=0; locus < probsF_.size(); ++locus){
    
    // size inner vector
    mendelian_probs_f[locus].resize(4);
    mendelian_alleles_f[locus].resize(4);
    mendelian_probs_m[locus].resize(4);
    mendelian_alleles_m[locus].resize(4);
    
    // loop over possible alleles at that locus
    for(size_t allele : {0,1,2,3}){
      /**********
       * Female
       **********/
      // convert things in holders
      holdA = Rcpp::as<sVec>(allelesF_[locus][allele]);
      
      // set things in final reference
      mendelian_probs_f[locus][allele].insert(mendelian_probs_f[locus][allele].end(),
                                              probsF_[locus][allele].begin(),
                                              probsF_[locus][allele].end());

      mendelian_alleles_f[locus][allele].insert(mendelian_alleles_f[locus][allele].end(),
                                                holdA.begin(),
                                                holdA.end());
      
      /**********
       * Male
       **********/
      // convert things in holders
      holdA = Rcpp::as<sVec>(allelesM_[locus][allele]);
      
      // set things in final reference
      mendelian_probs_m[locus][allele].insert(mendelian_probs_m[locus][allele].end(),
                                              probsM_[locus][allele].begin(),
                                              probsM_[locus][allele].end());
      
      mendelian_alleles_m[locus][allele].insert(mendelian_alleles_m[locus][allele].end(),
                                                holdA.begin(),
                                                holdA.end());
      
    }// end loop over alleles
  }// end loop over loci 

}

// getters for mendelian inheritance
dVec::iterator reference::get_mendelian_probs_begin(size_t sex, size_t locus, size_t allele){
  // 0 is female, 1 is male
  switch(sex){
    case 0:
      return mendelian_probs_f[locus][allele].begin();
    case 1:
      return mendelian_probs_m[locus][allele].begin();
  }
}
  
dVec::iterator reference::get_mendelian_probs_end(size_t sex, size_t locus, size_t allele){
  // 0 is female, 1 is male
  switch(sex){
  case 0:
    return mendelian_probs_f[locus][allele].end();
  case 1:
    return mendelian_probs_m[locus][allele].end();
  }
}

sVec::iterator reference::get_mendelian_allele_begin(size_t sex, size_t locus, size_t allele){
  // 0 is female, 1 is male
  switch(sex){
  case 0:
    return mendelian_alleles_f[locus][allele].begin();
  case 1:
    return mendelian_alleles_m[locus][allele].begin();
  }
}
  
sVec::iterator reference::get_mendelian_allele_end(size_t sex, size_t locus, size_t allele){
  // 0 is female, 1 is male
  switch(sex){
  case 0:
    return mendelian_alleles_f[locus][allele].end();
  case 1:
    return mendelian_alleles_m[locus][allele].end();
  }
}
 

// set homing inheritance
void reference::set_homing(const Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericVector> >& probsF_,
                           const Rcpp::ListOf<Rcpp::ListOf<Rcpp::StringVector> >& allelesF_,
                           const Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericVector> >& probsM_,
                           const Rcpp::ListOf<Rcpp::ListOf<Rcpp::StringVector> >& allelesM_){
  
  // set size of outer vector
  //  All should be the same size
  homing_probs_f.resize(probsF_.size());
  homing_alleles_f.resize(allelesF_.size());
  homing_probs_m.resize(probsF_.size());
  homing_alleles_m.resize(allelesF_.size());
  
  // holder objects to make life easier in the loop
  std::vector<std::string> holdA;
  
  
  // loop over number of loci
  for(size_t locus=0; locus < probsF_.size(); ++locus){
    
    // size inner vector
    homing_probs_f[locus].resize(4);
    homing_alleles_f[locus].resize(4);
    homing_probs_m[locus].resize(4);
    homing_alleles_m[locus].resize(4);
    
    // loop over possible alleles at that locus
    for(size_t allele : {0,1,2,3}){
      /**********
       * Female
       **********/
      // convert things in holders
      holdA = Rcpp::as<sVec>(allelesF_[locus][allele]);
      
      // set things in final reference
      homing_probs_f[locus][allele].insert(homing_probs_f[locus][allele].end(),
                                           probsF_[locus][allele].begin(),
                                           probsF_[locus][allele].end());
      
      homing_alleles_f[locus][allele].insert(homing_alleles_f[locus][allele].end(),
                                                holdA.begin(),
                                                holdA.end());
      
      /**********
       * Male
       **********/
      // convert things in holders
      holdA = Rcpp::as<sVec>(allelesM_[locus][allele]);
      
      // set things in final reference
      homing_probs_m[locus][allele].insert(homing_probs_m[locus][allele].end(),
                                           probsM_[locus][allele].begin(),
                                           probsM_[locus][allele].end());
      
      homing_alleles_m[locus][allele].insert(homing_alleles_m[locus][allele].end(),
                                             holdA.begin(),
                                             holdA.end());
      
    }// end loop over alleles
  }// end loop over loci 
   
}

// getters for homing inheritance
dVec::iterator reference::get_homing_probs_begin(size_t sex, size_t locus, size_t allele){
  // 0 is female, 1 is male
  switch(sex){
  case 0:
    return homing_probs_f[locus][allele].begin();
  case 1:
    return homing_probs_m[locus][allele].begin();
  }
}

dVec::iterator reference::get_homing_probs_end(size_t sex, size_t locus, size_t allele){
  // 0 is female, 1 is male
  switch(sex){
  case 0:
    return homing_probs_f[locus][allele].end();
  case 1:
    return homing_probs_m[locus][allele].end();
  }
}
  
sVec::iterator reference::get_homing_allele_begin(size_t sex, size_t locus, size_t allele){
  // 0 is female, 1 is male
  switch(sex){
  case 0:
    return homing_alleles_f[locus][allele].begin();
  case 1:
    return homing_alleles_m[locus][allele].begin();
  }
}
  
sVec::iterator reference::get_homing_allele_end(size_t sex, size_t locus, size_t allele){
  // 0 is female, 1 is male
  switch(sex){
  case 0:
    return homing_alleles_f[locus][allele].end();
  case 1:
    return homing_alleles_m[locus][allele].end();
  }
}
  

// set cutting inheritance
void reference::set_cutting(const Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericVector> >& probsF_,
                            const Rcpp::ListOf<Rcpp::ListOf<Rcpp::StringVector> >& allelesF_,
                            const Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericVector> >& probsM_,
                            const Rcpp::ListOf<Rcpp::ListOf<Rcpp::StringVector> >& allelesM_){
  
  // set size of outer vector
  //  All should be the same size
  cutting_probs_f.resize(probsF_.size());
  cutting_alleles_f.resize(allelesF_.size());
  cutting_probs_m.resize(probsF_.size());
  cutting_alleles_m.resize(allelesF_.size());
  
  // holder objects to make life easier in the loop
  std::vector<std::string> holdA;
  
  // loop over number of loci
  for(size_t locus=0; locus < probsF_.size(); ++locus){
    
    // size inner vector
    cutting_probs_f[locus].resize(4);
    cutting_alleles_f[locus].resize(4);
    cutting_probs_m[locus].resize(4);
    cutting_alleles_m[locus].resize(4);
    
    // loop over possible alleles at that locus
    for(size_t allele : {0,1,2,3}){
      /**********
       * Female
       **********/
      // convert things in holders
      holdA = Rcpp::as<sVec>(allelesF_[locus][allele]);
      
      // set things in final reference
      cutting_probs_f[locus][allele].insert(cutting_probs_f[locus][allele].end(),
                                              probsF_[locus][allele].begin(),
                                              probsF_[locus][allele].end());
      
      cutting_alleles_f[locus][allele].insert(cutting_alleles_f[locus][allele].end(),
                                                holdA.begin(),
                                                holdA.end());
      
      /**********
       * Male
       **********/
      // convert things in holders
      holdA = Rcpp::as<sVec>(allelesM_[locus][allele]);
      
      // set things in final reference
      cutting_probs_m[locus][allele].insert(cutting_probs_m[locus][allele].end(),
                                              probsM_[locus][allele].begin(),
                                              probsM_[locus][allele].end());
      
      cutting_alleles_m[locus][allele].insert(cutting_alleles_m[locus][allele].end(),
                                                holdA.begin(),
                                                holdA.end());
      
    }// end loop over alleles
  }// end loop over loci 
  
}

// getters for cutting inheritance
dVec::iterator reference::get_cutting_probs_begin(size_t sex, size_t locus, size_t allele){
  // 0 is female, 1 is male
  switch(sex){
  case 0:
    return cutting_probs_f[locus][allele].begin();
  case 1:
    return cutting_probs_m[locus][allele].begin();
  }
}
dVec::iterator reference::get_cutting_probs_end(size_t sex, size_t locus, size_t allele){
  // 0 is female, 1 is male
  switch(sex){
  case 0:
    return cutting_probs_f[locus][allele].end();
  case 1:
    return cutting_probs_m[locus][allele].end();
  }
}
  
sVec::iterator reference::get_cutting_allele_begin(size_t sex, size_t locus, size_t allele){
  // 0 is female, 1 is male
  switch(sex){
  case 0:
    return cutting_alleles_f[locus][allele].begin();
  case 1:
    return cutting_alleles_m[locus][allele].begin();
  }
}
  
sVec::iterator reference::get_cutting_allele_end(size_t sex, size_t locus, size_t allele){
  // 0 is female, 1 is male
  switch(sex){
  case 0:
    return cutting_alleles_f[locus][allele].end();
  case 1:
    return cutting_alleles_m[locus][allele].end();
  }
}
  
