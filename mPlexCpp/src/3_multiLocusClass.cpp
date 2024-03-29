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
multiLocus::multiLocus(const int& patchID_,
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
   twoAlleleGenotypes(reference::instance().get_alleloTypes(patchID), genos0, probs0);
   
   fillPopulation(patchID, genos0, probs0, 
                  eggs, larva, pupa,
                  adult_male, adult_female, unmated_female,
                  CreateMosquitoes);
   
   // Reproduction setup
   numAlleles = reference::instance().get_alleloTypes(patchID).size();
   fProbs.resize(numAlleles);
   mProbs.resize(numAlleles);
   fAllele.resize(numAlleles);
   mAllele.resize(numAlleles);

};

multiLocus::~multiLocus(){};

/******************************************************************************
 * Default (compiler-generated) move semantics
 ******************************************************************************/
multiLocus::multiLocus(multiLocus&& d) = default;
multiLocus& multiLocus::operator=(multiLocus&& d) = default;

/******************************************************************************
 * Reset
 ******************************************************************************/
void multiLocus::reset_Patch(){
  
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
void multiLocus::oneDay_layEggs(prng& myPRNG){
  
  for(auto female : adult_female){
    
    // calculate genotypes and probs of offspring
    MultiplexOffspring_mLoci(female.get_genotype(), female.get_mate());
    
    // get number of new offspring based on genotype and poisson randomness
    index = myPRNG.get_rpois(parameters::instance().get_beta()
                                              * reference::instance().get_s(female.get_genotype()));
    
    // pull eggs over offspring probs
    newEggs = myPRNG.get_rmultinom_online(index, finalProbs);
    
    // create new eggs
    for(size_t eggIndex=0; eggIndex<newEggs.size(); ++eggIndex){
      for(size_t it=0; it<newEggs[eggIndex]; ++it){
        eggs.emplace_back(Mosquito(0, finalGenotypes[eggIndex]));
      } // end loop over number of eggs per genotype
    } // end loop over newEggs vector
    
  } // end loop over females
  
}

/******************************************************************************
 * FUNCTIONS
 ******************************************************************************/

/**************************************
 * GENERATING FUNCTION
 **************************************/
void multiLocus::MultiplexOffspring_mLoci(const std::string& fGen, const std::string& mGen){



  //// list of objects I create every time ////
  // vector<bool> fScore
  // vector<bool> mScore






  // get number of alleles
  // don't need to reset because it's redefined every time
  // numAlleles = fGen.size()/2;

  /*****************************************************************************/
  // Score Each Allele
  /*****************************************************************************/
  std::vector<bool> fScore(numAlleles, false), mScore(numAlleles, false);

  //loop over loci, separate alleles and score
  index=0;
  for(size_t i=0; i<fGen.size(); i+=2, ++index){
    // female score
    if( (fGen[i] == 'H') || (fGen[i+1] == 'H')) {fScore[index] = true;}
    // male score
    if( (mGen[i] == 'H') || (mGen[i+1] == 'H')) {mScore[index] = true;}

  } // end scoring loop
  /*****************************************************************************/
  // End Score Each Allele
  /*****************************************************************************/

  /*****************************************************************************/
  // Determine Next-Gen alleles
  /*****************************************************************************/
  // these get reused, so clear them first
  for(index=0; index < numAlleles; ++index){
    fProbs[index].clear();
    mProbs[index].clear();
    fAllele[index].clear();
    mAllele[index].clear();
  }


  // loop over all loci
  index=0;
  for(size_t i=0; i<numAlleles; ++i, index+=2){

    // FEMALES
    if( fScore[i] ){
      // homing
      // loop over alleles, all loci diploid
      for(size_t j=0; j<2; ++j){
        // fill allele and probs
        if(fGen[index+j] == 'W'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_homing_allele_begin(0,i,0),
                            reference::instance().get_homing_allele_end(0,i,0));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_homing_probs_begin(0,i,0),
                           reference::instance().get_homing_probs_end(0,i,0));
        } else if(fGen[index+j] == 'H'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_homing_allele_begin(0,i,1),
                            reference::instance().get_homing_allele_end(0,i,1));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_homing_probs_begin(0,i,1),
                           reference::instance().get_homing_probs_end(0,i,1));
        } else if(fGen[index+j] == 'R'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_homing_allele_begin(0,i,2),
                            reference::instance().get_homing_allele_end(0,i,2));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_homing_probs_begin(0,i,2),
                           reference::instance().get_homing_probs_end(0,i,2));
        } else if(fGen[index+j] == 'S'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_homing_allele_begin(0,i,3),
                            reference::instance().get_homing_allele_end(0,i,3));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_homing_probs_begin(0,i,3),
                           reference::instance().get_homing_probs_end(0,i,3));
        }

      } // end allele loop
    } else {
      // no homing
      // loop over alleles, all loci diploid
      for(size_t j=0; j<2; ++j){
        // fill allele and probs
        if(fGen[index+j] == 'W'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_mendelian_allele_begin(0,i,0),
                            reference::instance().get_mendelian_allele_end(0,i,0));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_mendelian_probs_begin(0,i,0),
                           reference::instance().get_mendelian_probs_end(0,i,0));
        } else if(fGen[index+j] == 'H'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_mendelian_allele_begin(0,i,1),
                            reference::instance().get_mendelian_allele_end(0,i,1));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_mendelian_probs_begin(0,i,1),
                           reference::instance().get_mendelian_probs_end(0,i,1));
        } else if(fGen[index+j] == 'R'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_mendelian_allele_begin(0,i,2),
                            reference::instance().get_mendelian_allele_end(0,i,2));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_mendelian_probs_begin(0,i,2),
                           reference::instance().get_mendelian_probs_end(0,i,2));
        } else if(fGen[index+j] == 'S'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_mendelian_allele_begin(0,i,3),
                            reference::instance().get_mendelian_allele_end(0,i,3));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_mendelian_probs_begin(0,i,3),
                           reference::instance().get_mendelian_probs_end(0,i,3));
        }

      } // end allele loop
    } // end females


    // MALES
    if( mScore[i] ){
      // homing
      // loop over alleles, all loci diploid
      for(size_t j=0; j<2; ++j){
        // fill allele and probs
        if(mGen[index+j] == 'W'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_homing_allele_begin(1,i,0),
                            reference::instance().get_homing_allele_end(1,i,0));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_homing_probs_begin(1,i,0),
                           reference::instance().get_homing_probs_end(1,i,0));
        } else if(mGen[index+j] == 'H'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_homing_allele_begin(1,i,1),
                            reference::instance().get_homing_allele_end(1,i,1));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_homing_probs_begin(1,i,1),
                           reference::instance().get_homing_probs_end(1,i,1));
        } else if(mGen[index+j] == 'R'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_homing_allele_begin(1,i,2),
                            reference::instance().get_homing_allele_end(1,i,2));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_homing_probs_begin(1,i,2),
                           reference::instance().get_homing_probs_end(1,i,2));
        } else if(mGen[index+j] == 'S'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_homing_allele_begin(1,i,3),
                            reference::instance().get_homing_allele_end(1,i,3));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_homing_probs_begin(1,i,3),
                           reference::instance().get_homing_probs_end(1,i,3));
        }

      } // end allele loop
    } else {
      // no homing
      // loop over alleles, all loci diploid
      for(size_t j=0; j<2; ++j){
        // fill allele and probs
        if(mGen[index+j] == 'W'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_mendelian_allele_begin(1,i,0),
                            reference::instance().get_mendelian_allele_end(1,i,0));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_mendelian_probs_begin(1,i,0),
                           reference::instance().get_mendelian_probs_end(1,i,0));
        } else if(mGen[index+j] == 'H'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_mendelian_allele_begin(1,i,1),
                            reference::instance().get_mendelian_allele_end(1,i,1));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_mendelian_probs_begin(1,i,1),
                           reference::instance().get_mendelian_probs_end(1,i,1));
        } else if(mGen[index+j] == 'R'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_mendelian_allele_begin(1,i,2),
                            reference::instance().get_mendelian_allele_end(1,i,2));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_mendelian_probs_begin(1,i,2),
                           reference::instance().get_mendelian_probs_end(1,i,2));
        } else if(mGen[index+j] == 'S'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_mendelian_allele_begin(1,i,3),
                            reference::instance().get_mendelian_allele_end(1,i,3));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_mendelian_probs_begin(1,i,3),
                           reference::instance().get_mendelian_probs_end(1,i,3));
        }

      } // end allele loop
    } // end males

  } // end loci loop

  /*****************************************************************************/
  // End Next-Gen alleles
  /*****************************************************************************/

  /*****************************************************************************/
  // All Combinations of Male/Female for each Loci
  /*****************************************************************************/

  // Aggregate duplicates before storing
  // these get reused, so save them
  duplicates.clear();
  // value
  // holdAllele is defined in class, always gets reset


  // loop over alleles
  for(size_t i=0; i<numAlleles; ++i){

    // all combinations of female/male alleles, and probs
    for(size_t fem=0; fem<fAllele[i].size(); ++fem){
      for(size_t mal=0; mal<mAllele[i].size(); ++mal){
        // combine and sort allele
        holdAllele = fAllele[i][fem] + mAllele[i][mal];
        std::sort(holdAllele.begin(), holdAllele.end() );

        // Here we aggregate non-unique values
        // check if it is in map
        value = duplicates.find(holdAllele);
        if(value == duplicates.end()){
          // not in map, add it
          duplicates.insert(std::make_pair(holdAllele, fProbs[i][fem]*mProbs[i][mal]));
        } else {
          // is in map, combine with existing value
          value->second += fProbs[i][fem]*mProbs[i][mal];
        }

      } // end male loop
    } // end female loop

    // clear fAllele and fProbs for re-use
    fAllele[i].clear();
    fProbs[i].clear();

    // Take unique values in map, return to matrix
    for(auto elem : duplicates){
      fAllele[i].push_back(elem.first);
      fProbs[i].push_back(elem.second);
    }

    // clear map
    duplicates.clear();
  } // end

  /*****************************************************************************/
  // End All Combinations of MAle/Female for Each Loci
  /*****************************************************************************/

  /*****************************************************************************/
  // Cartesian Product of All Loci
  /*****************************************************************************/
  // these get reused, clear them
  finalGenotypes.clear();
  finalProbs.clear();

  holdGens.clear();
  holdProbs.clear();
  
  // fill first index
  for(index=0; index<fAllele[0].size(); ++index){
    finalGenotypes.push_back(fAllele[0][index]);
    finalProbs.push_back(fProbs[0][index]);
  }

  // loop over rest of them
  for(index=1; index<numAlleles; ++index){
    // loop over current elements in return, and elements in next vector
    for(size_t current=0; current<finalGenotypes.size(); ++current){
      for(size_t newOne=0; newOne<fAllele[index].size(); ++newOne){
        // combine strings and probabilities
        holdGens.push_back(finalGenotypes[current] + fAllele[index][newOne]);
        holdProbs.push_back(finalProbs[current] * fProbs[index][newOne]);
      } // end loop over new things
    } // end loop over current things

    // set returns
    finalGenotypes = holdGens;
    finalProbs = holdProbs;

    // clear holders
    holdGens.clear();
    holdProbs.clear();
  } // end loop over loci

  /*****************************************************************************/
  // End Cartesian Product of All Loci
  /*****************************************************************************/

  
  
  // test normalize
  double normalizer = std::accumulate(finalProbs.begin(), finalProbs.end(),0.0);
  
  for(auto& it : finalProbs){
    it /= normalizer;
  }

} // end function










