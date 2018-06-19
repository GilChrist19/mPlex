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
  
  Rcpp::Rcout << "\nMade it into set_reference!\n";
  
  // set genotype specific mating fitness, eta
  if(eta_.size() !=0 ){
    listNames = Rcpp::as<sVec>(eta_.names());
    for(size_t it=0; it<eta_.size(); ++it){ // for all of them, fill the map, casting things as key,value pairs
      eta.insert(std::make_pair(listNames[it], eta_[it]) );
    }
  }

  Rcpp::Rcout << "\tPassed eta\n";
  
  // set genotype specific sex ratio at emergence, phi
  if(phi_.size() != 0){
    listNames = Rcpp::as<sVec>(phi_.names());
    for(size_t it=0; it<phi_.size(); ++it){
      phi.insert(std::make_pair(listNames[it], phi_[it]) );
    }
  }

  Rcpp::Rcout << "\tPassed phi\n";  
  
  // set genotype specific multiplicative modifier of adult mortality, omegs
  if(omega_.size() != 0){
    listNames = Rcpp::as<sVec>(omega_.names());
    for(size_t it=0; it<omega_.size(); ++it){
      omega.insert(std::make_pair(listNames[it], omega_[it]) );
    }
  }
  
  Rcpp::Rcout << "\tPassed omega\n";

  // set genotype specific female pupatory success, xiF
  if(xiF_.size() != 0){
    listNames = Rcpp::as<sVec>(xiF_.names());
    for(size_t it=0; it<xiF_.size(); ++it){
      xiF.insert(std::make_pair(listNames[it], xiF_[it]) );
    }
  }
  
  Rcpp::Rcout << "\tPassed xiF\n";

  //set genotype specific male pupatory success, xiM
  if(xiM_.size() != 0){
    listNames = Rcpp::as<sVec>(xiM_.names());
    for(size_t it=0; it<xiM_.size(); ++it){
      xiM.insert(std::make_pair(listNames[it], xiM_[it]) );
    }
  }
  
  Rcpp::Rcout << "\tPassed xiM\n";

  // set fertility fraction, s
  if(s_.size() != 0){
    listNames = Rcpp::as<sVec>(s_.names());
    for(size_t it=0; it<s_.size(); ++it){
      s.insert(std::make_pair(listNames[it], s_[it]) );
    }
  }
  
  Rcpp::Rcout << "\tPassed s\n";

};


// get genotype dependent parameters
double reference::get_eta(std::string genType){
  
  double hold = 1.0;
  
  // iterator to element if it exists in the map
  std::unordered_map<std::string, double>::iterator it = eta.find(genType);
  
  // if it doesn't exist, it returns the end of the map, so check that
  if(it != eta.end()){
    hold = it->second;
  }
  
  // return
  return hold;
}

double reference::get_phi(std::string genType){
  
  double hold = 0.5;
  
  // iterator to element if it exists in the map
  std::unordered_map<std::string, double>::iterator it = phi.find(genType);
  
  // if it doesn't exist, it returns the end of the map, so check that
  if(it != phi.end()){
    hold = it->second;
  }
  
  // return
  return hold;
}

double reference::get_omega(std::string genType){
  
  double hold = 1.0;
  
  // iterator to element if it exists in the map
  std::unordered_map<std::string, double>::iterator it = omega.find(genType);
  
  // if it doesn't exist, it returns the end of the map, so check that
  if(it != omega.end()){
    hold = it->second;
  }
  
  // return
  return hold;
}

double reference::get_xiF(std::string genType){
  
  double hold = 1.0;
  
  // iterator to element if it exists in the map
  std::unordered_map<std::string, double>::iterator it = xiF.find(genType);
  
  // if it doesn't exist, it returns the end of the map, so check that
  if(it != xiF.end()){
    hold = it->second;
  }
  
  // return
  return hold;
}

double reference::get_xiM(std::string genType){
  
  double hold = 1.0;
  
  // iterator to element if it exists in the map
  std::unordered_map<std::string, double>::iterator it = xiM.find(genType);
  
  // if it doesn't exist, it returns the end of the map, so check that
  if(it != xiM.end()){
    hold = it->second;
  }
  
  // return
  return hold;
}

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

void reference::set_mendelian(const Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericVector> >& probs_,
                              const Rcpp::ListOf<Rcpp::ListOf<Rcpp::StringVector> >& alleles_){
  
  Rcpp::Rcout << "Made it into set_mendelian\n";
  
  Rcpp::Rcout << "probs size is " << probs_.size() << std::endl;
  
  // There are always 4 alleles possible
  std::vector<int> alleles = {0,1,2,3};
  
  homing_probs.resize(probs_.size()); // size of outer vector
  homing_alleles.resize(alleles_.size());; // size outer vector
   
   
  Rcpp::Rcout << "Setup homing probs and alleles\n";
  
  Rcpp::Rcout << "dim of homing probs is: " << homing_probs.size() << std::endl;
  Rcpp::Rcout << homing_probs[0].size() << std::endl;
//  Rcpp::Rcout << homing_probs[0][0].size() << std::endl;
  
  
  //   
  // using dVec = std::vector<double>; 
  // using sVec = std::vector<std::string>;
  // 
  // using dArVec = std::vector<dVec>;
  // using sArVec = std::vector<sVec>;
  //   
  // (probs_.size(), std::vector<dVec>(4));
  
  Rcpp::NumericVector holdP;
  Rcpp::StringVector holdA;
  
  
  // homing_probs.reserve(probs_.size());
  // homing_alleles.reserve(probs_.size());
  
  // loop over number of loci
  for(size_t locus=0; locus < probs_.size(); ++locus){
    
    Rcpp::Rcout << "In loop iter: "<<locus<<std::endl;
    
    // size inner vector
    homing_probs[locus].resize(4);
    homing_alleles[locus].resize(4);
    
    Rcpp::Rcout << "dim of homing probs is: " << homing_probs.size() << std::endl;
    Rcpp::Rcout << homing_probs[0].size() << std::endl;
    Rcpp::Rcout << homing_probs[0][0].size() << std::endl;
    
    // // create place to put things
    //homing_probs.push_back();
    // homing_alleles.push_back(std::vector<sVec>);
    

    
    
    // loop over possible alleles at that locus
    for(size_t allele : alleles){
      
      Rcpp::Rcout << "Inner loop iter: " << allele << std::endl;
      
      
      holdP = Rcpp::as<dVec>(probs_[locus][allele]);
      holdA = Rcpp::as<sVec>(alleles_[locus][allele]);
      
      //mendelian_probs[locus][allele].reserve(holdP.size());
      //mendelian_alleles[locus][allele].reserve(holdP.size());
      
      // mendelian_probs[locus][allele] = Rcpp::as<dVec>(probs_[locus][allele]);
      // mendelian_alleles[locus][allele] = Rcpp::as<sVec>(alleles_[locus][allele]);
      
      for(auto it=0; it < holdP.size(); ++ it){
        
        mendelian_probs[locus][allele].emplace_back(holdP[it]);
        mendelian_alleles[locus][allele].emplace_back(holdA[it]);
        
        
      }

      
      
      // mendelian_probs[locus][allele].insert(mendelian_probs[locus][allele].end(),
      //                                       holdP.begin(),
      //                                       holdP.end());
      // 
      // mendelian_alleles[locus][allele].insert(mendelian_alleles[locus][allele].end(),
      //                                         holdA.begin(),
      //                                         holdA.end());

      
      
      
      // mendelian_probs[locus][allele] = Rcpp::as<dVec>(probs_[locus][allele]);
      // mendelian_alleles[locus][allele] = Rcpp::as<sVec>(alleles_[locus][allele]);
    }
  }
  
  Rcpp::Rcout << "Finished mendelian stuffs" << std::endl;

}
 
 
void reference::set_homing(const Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericVector> >& probs_,
                            const Rcpp::ListOf<Rcpp::ListOf<Rcpp::StringVector> >& alleles_){
   
   // There are always 4 alleles possible
   std::vector<int> alleles = {0,1,2,3};
   
   // loop over number of loci
   for(size_t locus=0; locus < probs_.size(); ++locus){
     // loop over possible alleles at that locus
     for(size_t allele : alleles){
       homing_probs[locus].push_back( Rcpp::as<dVec>(probs_[locus][allele]) );
       homing_alleles[locus].push_back( Rcpp::as<sVec>(alleles_[locus][allele]) );
     }
   }
   
}



void reference::set_cutting(const Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericVector> >& probs_,
                            const Rcpp::ListOf<Rcpp::ListOf<Rcpp::StringVector> >& alleles_){
  
  // There are always 4 alleles possible
  std::vector<int> alleles = {0,1,2,3};
  
  // loop over number of loci
  for(size_t locus=0; locus < probs_.size(); ++locus){
    // loop over possible alleles at that locus
    for(size_t allele : alleles){
      cutting_probs[locus][allele] = Rcpp::as<dVec>(probs_[locus][allele]);
      cutting_alleles[locus][allele] = Rcpp::as<sVec>(alleles_[locus][allele]);
    }
  }
  
}
