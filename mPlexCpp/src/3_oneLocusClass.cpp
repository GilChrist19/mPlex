///////////////////////////////////////////////////////////////////////////////
//                           ____  __          ______          
//                ____ ___  / __ \/ /__  _  __/ ____/___  ____ 
//               / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
//              / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
//             /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/ 
//                                              /_/   /_/      
///////////////////////////////////////////////////////////////////////////////

#include "3_DriveDefinitions.hpp"


/******************************************************************************
 * Constructor & Destructor
 ******************************************************************************/
oneLocus::oneLocus(const int& patchID_,
                   const Rcpp::List& maleReleases_,
                   const Rcpp::List& femaleReleases_,
                   const Rcpp::List& eggReleases_,
                   const Rcpp::List& matedFemaleReleases_) : Patch::Patch(patchID_,
                                                                          maleReleases_,
                                                                          femaleReleases_,
                                                                          eggReleases_,
                                                                          matedFemaleReleases_)
{
 
  /****************
  * SET POPULATIONS
  ****************/
  // set initial genotypes and their probabilities
  twoLociGenotypes(reference::instance().get_alleloTypes(patchID), genos0, probs0);
  
  fillPopulation(patchID, genos0, probs0,
                 eggs, larva, pupa,
                 adult_male, adult_female, unmated_female,
                 CreateMosquitoes);
  
  
  // Reproduction setup
  numLoci = reference::instance().get_alleloTypes(patchID).size();
 
};

oneLocus::~oneLocus(){};

/******************************************************************************
 * Default (compiler-generated) move semantics
 ******************************************************************************/
oneLocus::oneLocus(oneLocus&& d) = default;
oneLocus& oneLocus::operator=(oneLocus&& d) = default;

/******************************************************************************
 * Reset
 ******************************************************************************/
void oneLocus::reset_Patch(){
  
  /****************
   * RESET POPULATIONS
   ****************/
  eggs.clear();
  larva.clear();
  pupa.clear();
  adult_male.clear();
  adult_female.clear();
  unmated_female.clear();
  
  /****************
   * SET POPULATIONS
   ****************/
  fillPopulation(patchID, genos0, probs0,
                 eggs, larva, pupa,
                 adult_male, adult_female, unmated_female,
                 CreateMosquitoes);
  
  /****************
   * RESET RELEASES
   ****************/
  releaseM = releaseM0;
  releaseF = releaseF0;
  releaseE = releaseE0;
  releaseMF = releaseMF0;
}

/******************************************************************************
 * Lay eggs
 ******************************************************************************/
void oneLocus::oneDay_layEggs(prng& myPRNG){

  for(auto female : adult_female){

    // calculate genotypes and probs of offspring
    MultiplexOffspring_oLocus(female.get_genotype(), female.get_mate());

    // get number of new offspring based on genotype and poisson randomness
    index = myPRNG.get_rpois(parameters::instance().get_beta()
                                              * reference::instance().get_s(female.get_genotype()));

    // pull eggs over offspring probs
    newEggs = myPRNG.get_rmultinom_online(index, holdProbs1);

    // create new eggs
    for(size_t eggIndex=0; eggIndex<newEggs.size(); ++eggIndex){
      for(size_t it=0; it<newEggs[eggIndex]; ++it){
        eggs.emplace_back(Mosquito(0, holdGens1[eggIndex]));
      } // end loop over number of eggs per genotype
    } // end loop over newEggs vector
  } // end loop over females
  
}

/******************************************************************************
 * FUNCTIONS
 ******************************************************************************/

/**************************************
 * SETUP
 **************************************/
void twoLociGenotypes(const Rcpp::ListOf<Rcpp::List>& aTypes, sVec& retGenos, dVec& retProbs){
  
  // This only good for things that have 2 loci, but multiple alleles

  // number of loci is the number of sublists in aTypes
  size_t numLoci = aTypes.size();

  // holders for things in loop
  sVec genotypes, holdGens, alleles;
  dVec genoProbs, holdProbs,probs;
  
  
  // fill first index
  probs = Rcpp::as<dVec>( aTypes[0]["probs"] );
  alleles = Rcpp::as<sVec>( aTypes[0]["alleles"] );
  
  for(size_t index=0; index < alleles.size(); ++index){
    genotypes.push_back(alleles[index]);
    genoProbs.push_back(probs[index]);
  }
  
  
  // loop over rest of them
  for(int index=1; index<numLoci; ++index){
    
    // get probabilities
    probs = Rcpp::as<dVec>( aTypes[index]["probs"] );
    alleles = Rcpp::as<sVec>( aTypes[index]["alleles"] );
    
    // loop over current elements in return, and elements in next vector
    for(size_t current=0; current<genotypes.size(); ++current){
      for(size_t newOne=0; newOne<alleles.size(); ++newOne){
        // combine strings and probabilities
        holdGens.push_back(genotypes[current] + alleles[newOne]);
        holdProbs.push_back(genoProbs[current] * probs[newOne]);
      } // end loop over new things
    } // end loop over current things
    
    // set returns
    genotypes = holdGens;
    genoProbs = holdProbs;

    // clear holders
    holdGens.clear();
    holdProbs.clear();
    
  } // end loop over loci
  
  
  // combine possible alleles at each locus
  twoAlleleLocus(genoProbs, genotypes, holdProbs, holdGens);
  
  
  // normalize and remove zeros
  double normalizer = std::accumulate(holdProbs.begin(), holdProbs.end(),0.0);
  for(size_t index=0; index<holdGens.size(); ++index){
    // if not zero, normalize and keep
    if(holdProbs[index] != 0){
      retGenos.push_back(holdGens[index]);
      retProbs.push_back(holdProbs[index]/normalizer);
    }
  } // end normalize and remove zeros
  
}


/**************************************
 * GENERATING FUNCTION
 **************************************/
void oneLocus::MultiplexOffspring_oLocus(const std::string& fGen, const std::string& mGen){




  // get number of alleles, this is reused every time and resets here
  //numAlleles = fGen.size()/2;
  

  /*****************************************************************************/
  // Score Each Allele
  /*****************************************************************************/
  // these get reused, and set here
  fScore = false;
  mScore = false;
  index = 0;

  //loop over loci, separate alleles and score
  while((!fScore || !mScore) && index < fGen.size()){
    // female score
    if(fGen[index] == 'H') {fScore = true;}
    // male score
    if(mGen[index] == 'H') {mScore = true;}
    // increment index
    ++index;
  }
  
  
  /*****************************************************************************/
  // End Score Each Allele
  /*****************************************************************************/
  
  /*****************************************************************************/
  // Determine Next-Gen alleles
  /*****************************************************************************/
  // reset things that are reused
  index=0;

  holdGens1.clear();
  holdProbs1.clear();

  fAllele.clear();
  fProbs.clear();
  mAllele.clear();
  mProbs.clear();


  // FEMALES
  if(fScore){
    // TRUE - there is homing
    // loop over alleles, only two
    for(size_t allele=0; allele<2; ++allele, index+=numLoci){
      // loop over each target in allele
      for(size_t i=0; i<numLoci; ++i){
        // fill allele and probs
        if(fGen[i+index] == 'W'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_homing_allele_begin(i,0),
                            reference::instance().get_homing_allele_end(i,0));
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_homing_probs_begin(i,0),
                            reference::instance().get_homing_probs_end(i,0));
        } else if(fGen[i+index] == 'H'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_homing_allele_begin(i,1),
                            reference::instance().get_homing_allele_end(i,1));
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_homing_probs_begin(i,1),
                            reference::instance().get_homing_probs_end(i,1));
        } else if(fGen[i+index] == 'R'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_homing_allele_begin(i,2),
                            reference::instance().get_homing_allele_end(i,2));
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_homing_probs_begin(i,2),
                           reference::instance().get_homing_probs_end(i,2));
        } else if(fGen[i+index] == 'S'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_homing_allele_begin(i,3),
                            reference::instance().get_homing_allele_end(i,3));
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_homing_probs_begin(i,3),
                           reference::instance().get_homing_probs_end(i,3));
        }
        
        /********************************
        // All Combinations of Target Sites
        *********************************/
        // make combinations
        for(size_t first=0; first < holdGens1.size(); ++first){
          for(size_t second=0; second < holdGens2.size(); ++second){
            holdGens3.push_back(holdGens1[first] + holdGens2[second]);
            holdProbs3.push_back(holdProbs1[first] * holdProbs2[second]);
          }
        }

        // set combinations in return matrix
        //  Is there a better way to set the first one?
        if(holdGens1.size() == 0){
          holdGens1 = holdGens2;
          holdProbs1 = holdProbs2;
        } else {
          holdGens1 = holdGens3;
          holdProbs1 = holdProbs3;
        }

        // clear hold values
        holdGens2.clear();
        holdProbs2.clear();
        holdGens3.clear();
        holdProbs3.clear();
        /********************************
        // End All Combinations of Target Sites
        *********************************/
      } // end loop over target sites at each allele

      /********************************
      // "Unlist"
      *********************************/
      // put finished set of alleles from each locus into a final vector, they
      //  can't mix with each other
      //  There could be duplicates here, do I care?
      fAllele.insert(fAllele.end(), holdGens1.begin(), holdGens1.end());
      fProbs.insert(fProbs.end(), holdProbs1.begin(), holdProbs1.end());

      // clear holder
      holdGens1.clear();
      holdProbs1.clear();

      /********************************
      // End "Unlist"
      *********************************/
    } // end allele loop

  } else {
    // FALSE - there is not homing
    // loop over alleles, only two
    for(size_t allele=0; allele<2; ++allele, index+=numLoci){
      // loop over each target in allele
      for(size_t i=0; i<numLoci; ++i){
        // fill allele and probs
        if(fGen[i+index] == 'W'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_mendelian_allele_begin(i,0),
                                 reference::instance().get_mendelian_allele_end(i,0));
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_mendelian_probs_begin(i,0),
                                reference::instance().get_mendelian_probs_end(i,0));
        } else if(fGen[i+index] == 'H'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_mendelian_allele_begin(i,1),
                                 reference::instance().get_mendelian_allele_end(i,1));
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_mendelian_probs_begin(i,1),
                                reference::instance().get_mendelian_probs_end(i,1));
        } else if(fGen[i+index] == 'R'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_mendelian_allele_begin(i,2),
                                 reference::instance().get_mendelian_allele_end(i,2));
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_mendelian_probs_begin(i,2),
                                reference::instance().get_mendelian_probs_end(i,2));
        } else if(fGen[i+index] == 'S'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_mendelian_allele_begin(i,3),
                                 reference::instance().get_mendelian_allele_end(i,3));
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_mendelian_probs_begin(i,3),
                                reference::instance().get_mendelian_probs_end(i,3));
        }
        
        /********************************
        // All Combinations of Target Sites
        *********************************/
        // make combinations
        for(size_t first=0; first < holdGens1.size(); ++first){
          for(size_t second=0; second < holdGens2.size(); ++second){
            holdGens3.push_back(holdGens1[first] + holdGens2[second]);
            holdProbs3.push_back(holdProbs1[first] * holdProbs2[second]);
          }
        }
        
        // set combinations in return matrix
        //  Is there a better way to set the first one?
        if(holdGens1.size() == 0){
          holdGens1 = holdGens2;
          holdProbs1 = holdProbs2;
        } else {
          holdGens1 = holdGens3;
          holdProbs1 = holdProbs3;
        }
        
        // clear hold values
        holdGens2.clear();
        holdProbs2.clear();
        holdGens3.clear();
        holdProbs3.clear();
        /********************************
        // End All Combinations of Target Sites
        *********************************/
      } // end loop over target sites at each allele

      /********************************
      // "Unlist"
      *********************************/
      // put finished set of alleles from each locus into a final vector, they
      //  can't mix with each other
      //  There could be duplicates here, do I care?
      fAllele.insert(fAllele.end(), holdGens1.begin(), holdGens1.end());
      fProbs.insert(fProbs.end(), holdProbs1.begin(), holdProbs1.end());

      // clear holder
      holdGens1.clear();
      holdProbs1.clear();

      /********************************
      // End "Unlist"
      *********************************/
    } // end allele loop

  } // end female if statement

  // reset index
  index=0;

  // MALES
  if(mScore){
    // TRUE - there is homing
    // loop over alleles, only two
    for(size_t allele=0; allele<2; ++allele, index+=numLoci){
      // loop over each target in allele
      for(size_t i=0; i<numLoci; ++i){
        // fill allele and probs
        if(mGen[i+index] == 'W'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_homing_allele_begin(i,0),
                                 reference::instance().get_homing_allele_end(i,0));
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_homing_probs_begin(i,0),
                                reference::instance().get_homing_probs_end(i,0));
        } else if(mGen[i+index] == 'H'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_homing_allele_begin(i,1),
                                 reference::instance().get_homing_allele_end(i,1));
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_homing_probs_begin(i,1),
                                reference::instance().get_homing_probs_end(i,1));
        } else if(mGen[i+index] == 'R'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_homing_allele_begin(i,2),
                                 reference::instance().get_homing_allele_end(i,2));
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_homing_probs_begin(i,2),
                                reference::instance().get_homing_probs_end(i,2));
        } else if(mGen[i+index] == 'S'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_homing_allele_begin(i,3),
                                 reference::instance().get_homing_allele_end(i,3));
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_homing_probs_begin(i,3),
                                reference::instance().get_homing_probs_end(i,3));
        }

        /********************************
        // All Combinations of Target Sites
        *********************************/
        // make combinations
        for(size_t first=0; first < holdGens1.size(); ++first){
          for(size_t second=0; second < holdGens2.size(); ++second){
            holdGens3.push_back(holdGens1[first] + holdGens2[second]);
            holdProbs3.push_back(holdProbs1[first] * holdProbs2[second]);
          }
        }

        // set combinations in return matrix
        //  Is there a better way to set the first one?
        if(holdGens1.size() == 0){
          holdGens1 = holdGens2;
          holdProbs1 = holdProbs2;
        } else {
          holdGens1 = holdGens3;
          holdProbs1 = holdProbs3;
        }

        // clear hold values
        holdGens2.clear();
        holdProbs2.clear();
        holdGens3.clear();
        holdProbs3.clear();
        /********************************
        // End All Combinations of Target Sites
        *********************************/
      } // end loop over target sites at each allele

      /********************************
      // "Unlist"
      *********************************/
      // put finished set of alleles from each locus into a final vector, they
      //  can't mix with each other
      //  There could be duplicates here, do I care?
      mAllele.insert(mAllele.end(), holdGens1.begin(), holdGens1.end());
      mProbs.insert(mProbs.end(), holdProbs1.begin(), holdProbs1.end());

      // clear holder
      holdGens1.clear();
      holdProbs1.clear();

      /********************************
      // End "Unlist"
      *********************************/
    } // end allele loop

  } else {

    // FALSE - there is no homing
    // loop over alleles, only two
    for(size_t allele=0; allele<2; ++allele, index+=numLoci){
      // loop over each target in allele
      for(size_t i=0; i<numLoci; ++i){
        // fill allele and probs
        if(mGen[i+index] == 'W'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_mendelian_allele_begin(i,0),
                                 reference::instance().get_mendelian_allele_end(i,0));
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_mendelian_probs_begin(i,0),
                                reference::instance().get_mendelian_probs_end(i,0));
        } else if(mGen[i+index] == 'H'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_mendelian_allele_begin(i,1),
                                 reference::instance().get_mendelian_allele_end(i,1));
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_mendelian_probs_begin(i,1),
                                reference::instance().get_mendelian_probs_end(i,1));
        } else if(mGen[i+index] == 'R'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_mendelian_allele_begin(i,2),
                                 reference::instance().get_mendelian_allele_end(i,2));
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_mendelian_probs_begin(i,2),
                                reference::instance().get_mendelian_probs_end(i,2));
        } else if(mGen[i+index] == 'S'){
          holdGens2.insert(holdGens2.end(), reference::instance().get_mendelian_allele_begin(i,3),
                                 reference::instance().get_mendelian_allele_end(i,3));
          holdProbs2.insert(holdProbs2.end(), reference::instance().get_mendelian_probs_begin(i,3),
                                reference::instance().get_mendelian_probs_end(i,3));
        }

        /********************************
        // All Combinations of Target Sites
        *********************************/
        // make combinations
        for(size_t first=0; first < holdGens1.size(); ++first){
          for(size_t second=0; second < holdGens2.size(); ++second){
            holdGens3.push_back(holdGens1[first] + holdGens2[second]);
            holdProbs3.push_back(holdProbs1[first] * holdProbs2[second]);
          }
        }

        // set combinations in return matrix
        //  Is there a better way to set the first one?
        if(holdGens1.size() == 0){
          holdGens1 = holdGens2;
          holdProbs1 = holdProbs2;
        } else {
          holdGens1 = holdGens3;
          holdProbs1 = holdProbs3;
        }

        // clear hold values
        holdGens2.clear();
        holdProbs2.clear();
        holdGens3.clear();
        holdProbs3.clear();
        /********************************
        // End All Combinations of Target Sites
        *********************************/
      } // end loop over target sites at each allele

      /********************************
      // "Unlist"
      *********************************/
      // put finished set of alleles from each locus into a final vector, they
      //  can't mix with each other
      //  There could be duplicates here, do I care?
      mAllele.insert(mAllele.end(), holdGens1.begin(), holdGens1.end());
      mProbs.insert(mProbs.end(), holdProbs1.begin(), holdProbs1.end());

      // clear holder
      holdGens1.clear();
      holdProbs1.clear();

      /********************************
      // End "Unlist"
      *********************************/
    } // end allele loop

  } // end male if statement
  /*****************************************************************************/
  // End Next-Gen alleles
  /*****************************************************************************/

  /*****************************************************************************/
  // Cartesian Product of male/female loci
  /*****************************************************************************/
  // These get reused, but two don't need cleared
  // holdAllele
  duplicates.clear();
  // std::unordered_map<std::string, double>::iterator value;


  // all combinations of female/male alleles and probs
  for(size_t fem=0; fem<fAllele.size(); ++fem){
    for(size_t mal=0; mal<mAllele.size(); ++mal){

      // combine alleles
      holdAllele = fAllele[fem] + mAllele[mal];

      // Here we aggregate non-unique values
      // check if it is in map
      value = duplicates.find(holdAllele);
      if(value == duplicates.end()){
        // not in map, add it
        duplicates.insert(std::make_pair(holdAllele, fProbs[fem] * mProbs[mal]));
      } else {
        // is in map, combine with existing value
        value->second += fProbs[fem]*mProbs[mal];
      }

    } // end male loop
  } //  end female loop

  // clear for reuse
  holdGens1.clear();
  holdProbs1.clear();

  // Take unique values in map, return to matrix
  for(auto elem : duplicates){
    holdGens1.push_back(elem.first);
    holdProbs1.push_back(elem.second);
  }

  /*****************************************************************************/
  // End Cartesian Product of male/female loci
  /*****************************************************************************/

  
    // test normalize
  double normalizer = std::accumulate(holdProbs1.begin(), holdProbs1.end(),0.0);
  
  for(auto& it : holdProbs1){
    it /= normalizer;
  }
  

} // end function








