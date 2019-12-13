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
Daisy::Daisy(const int& patchID_,
             const Rcpp::List& maleReleases_,
             const Rcpp::List& femaleReleases_,
             const Rcpp::List& eggReleases_) : Patch::Patch(patchID_,
                                                             maleReleases_,
                                                             femaleReleases_,
                                                             eggReleases_)
{
             
  /****************
  * SET POPULATIONS
  ****************/
  fillPopulation(patchID, eggs, larva, pupa,
                 adult_male, adult_female, unmated_female,
                 CreateMosquitoes2Allele);

  
  // Reproduction setup
  numAlleles = reference::instance().get_alleloTypes(patchID).size();
  fProbs.resize(numAlleles);
  mProbs.resize(numAlleles);
  fAllele.resize(numAlleles);
  mAllele.resize(numAlleles);
  
};

Daisy::~Daisy(){};

/******************************************************************************
 * Default (compiler-generated) move semantics
******************************************************************************/
Daisy::Daisy(Daisy&& d) = default;
Daisy& Daisy::operator=(Daisy&& d) = default;

/******************************************************************************
 * Reset
******************************************************************************/
void Daisy::reset_Patch(){
  
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
  fillPopulation(patchID, eggs, larva, pupa,
                 adult_male, adult_female, unmated_female,
                 CreateMosquitoes2Allele);
  
  /****************
   * RESET RELEASES
   ****************/
  releaseM = releaseM0;
  releaseF = releaseF0;
  releaseE = releaseE0;
}


/******************************************************************************
 * Lay eggs
******************************************************************************/
void Daisy::oneDay_layEggs(){

  for(auto female : adult_female){

    // calculate genotypes and probs of offspring
    DaisyOffspring(female.get_genotype(), female.get_mate());
    
    // get number of new offspring based on genotype and poisson randomness
    index = prng::instance().get_rpois(parameters::instance().get_beta()
                                              * reference::instance().get_s(female.get_genotype()));

    // pull eggs over offspring probs
    newEggs = prng::instance().get_rmultinom_online(index, finalProbs);

    
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
 * SETUP
 **************************************/
std::vector<double> markovDist(const double& life, const int& time){
  
  // setup matrix to solve
  arma::Mat<double> markovMat(time, time, arma::fill::eye);
  
  // create and fill offDiagonal vector
  arma::Col<double> offDiag(time-1);
  offDiag.fill((1.0 - life) - 1.0);
  
  // put off diagonal in matrix
  markovMat.diag(1) = offDiag;
  
  // invert matrix
  markovMat = markovMat.i();
  
  // get normalized vector of larval ratios
  arma::Row<double> solVec(markovMat.row(0)/arma::sum(markovMat.row(0)));
  
  // store as standard vector, both for return and because arma::Row doesn't 
  //  have some of the functions I need
  std::vector<double> retVec(solVec.begin(), solVec.end());
  
  return(retVec);
}

std::vector<double> popDist(const double& mu, const double& alpha, const int& larvalEQ, const iVec& timeAq) {
  // This comes from MGDrivE, where Sean solved the initial conditions as a CTMC
  
  /***********************************/
  // setup different death rates
  /***********************************/
  double epLife(1.0 - mu);
  double lLife(std::pow(alpha/(alpha + larvalEQ), 1.0/timeAq[1]) * epLife );
  
  /***********************************/
  // Solve larval distribution
  /***********************************/
  std::vector<double> hold(markovDist(lLife, timeAq[1]));
  
  /***********************************/
  // Egg distribution at the front
  /***********************************/
  for(int i = timeAq[0]; i > 0; i--){
    hold.insert(hold.begin(), hold.front()/epLife);
  }
  
  /***********************************/
  // Pupal distribution at the end
  /***********************************/
  hold.push_back(hold.back() * lLife);
  for(int i = timeAq[2]-1; i > 0; i--){
    hold.push_back(hold.back() * epLife);
  }
  
  // return normalized vector of initial conditions
  return(hold);
}

// this one is for Daisy, oneLocus, and multiLocus classes
void fillPopulation(const int& patchID, popVec& eggVec, popVec& larvaVec, popVec& pupaVec,
                    popVec& aMaleVec, popVec& aFemaleVec, popVec& unFemaleVec,
                    void (*populationFill)(const int&, const int&, const dVec&,
                          const Rcpp::ListOf<Rcpp::List>&, popVec&)
                      ){
  
  // holder objects, these are for distribution function
  int minAge, maxAge;
  dVec distHold;
  
  // solve aquatic distribution
  dVec aDist = popDist(parameters::instance().get_mu(0),
                       parameters::instance().get_alpha(patchID),
                       parameters::instance().get_larva_eq(patchID),
                       {parameters::instance().get_stage_time(0),
                        parameters::instance().get_stage_time(1),
                        parameters::instance().get_stage_time(2)});
  
  
  /*********/
  // eggs
  /*********/
  minAge = 0;
  maxAge = parameters::instance().get_stage_time(0) - 1;
  
  // set age distribution vector
  distHold.resize(maxAge+1);
  std::copy(aDist.begin(), aDist.begin() + maxAge+1, distHold.begin());
  
  // reserve estimated population size
  eggVec.reserve(round(1.1*parameters::instance().get_larva_eq(patchID) * accumulate(distHold.begin(), distHold.end(), 0.0)));
  
  // fill initial population
  populationFill(parameters::instance().get_larva_eq(patchID),
                 minAge, distHold, reference::instance().get_alleloTypes(patchID), eggVec);
  
  
  /*********/
  // larva
  /*********/
  minAge = parameters::instance().get_stage_time(0);
  maxAge = parameters::instance().get_stage_sum(1)-1;
  
  // set age distribution vector
  distHold.resize(maxAge+1 - minAge);
  std::copy(aDist.begin() + minAge, aDist.begin() + maxAge + 1, distHold.begin());
  
  // reserve estimated population size
  larvaVec.reserve(round(1.1*parameters::instance().get_larva_eq(patchID) * accumulate(distHold.begin(), distHold.end(), 0.0)));
  
  // fill initial population
  populationFill(parameters::instance().get_larva_eq(patchID),
                 minAge, distHold, reference::instance().get_alleloTypes(patchID), larvaVec);

  
  /*********/
  // pupa
  /*********/
  minAge = parameters::instance().get_stage_sum(1);
  maxAge = parameters::instance().get_stage_sum(2) - 1;
  
  // set age distribution vector
  distHold.resize(maxAge+1 - minAge);
  std::copy(aDist.begin() + minAge, aDist.begin() + maxAge + 1, distHold.begin());
  
  // reserve estimated population size
  pupaVec.reserve(round(1.1*parameters::instance().get_larva_eq(patchID) * accumulate(distHold.begin(), distHold.end(), 0.0)));
  
  // fill initial population
  populationFill(parameters::instance().get_larva_eq(patchID),
                 minAge, distHold, reference::instance().get_alleloTypes(patchID), pupaVec);
  
  
  /***********************************/
  // Solve adult distribution
  /***********************************/
  minAge = parameters::instance().get_stage_sum(2);
  
  // set age distribution vector
  aDist = markovDist(1.0 - parameters::instance().get_mu(3), parameters::instance().get_stage_time(3) * 3);
  
  // reserve estimated population size
  int popSize(round(1.1 * parameters::instance().get_adult_pop_eq(patchID) * accumulate(aDist.begin(), aDist.end(), 0.0)));
  aMaleVec.reserve(popSize);
  aFemaleVec.reserve(popSize);
  unFemaleVec.reserve(popSize);
  
  // fill initial population
  populationFill(parameters::instance().get_adult_pop_eq(patchID)/2,
                 minAge, aDist, reference::instance().get_alleloTypes(patchID), aMaleVec);
  populationFill(parameters::instance().get_adult_pop_eq(patchID)/2,
                 minAge, aDist, reference::instance().get_alleloTypes(patchID), unFemaleVec);
  
}

// this one is solely for family class
void fillPopulation(const int& patchID, popVec& eggVec, popVec& larvaVec, popVec& pupaVec,
                    popVec& aMaleVec, popVec& aFemaleVec, popVec& unFemaleVec,
                    void (*populationFill)(const int&, const int&, const dVec&, popVec&)
){
  
  // holder objects, these are for distribution function
  int minAge, maxAge;
  dVec distHold;
  
  // solve aquatic distribution
  dVec aDist = popDist(parameters::instance().get_mu(0),
                       parameters::instance().get_alpha(patchID),
                       parameters::instance().get_larva_eq(patchID),
                       {parameters::instance().get_stage_time(0),
                        parameters::instance().get_stage_time(1),
                        parameters::instance().get_stage_time(2)});
  
  
  /*********/
  // eggs
  /*********/
  minAge = 0;
  maxAge = parameters::instance().get_stage_time(0) - 1;
  
  // set age distribution vector
  distHold.resize(maxAge+1);
  std::copy(aDist.begin(), aDist.begin() + maxAge+1, distHold.begin());
  
  // reserve estimated population size
  eggVec.reserve(round(1.1*parameters::instance().get_larva_eq(patchID) * accumulate(distHold.begin(), distHold.end(), 0.0)));
  
  // fill initial population
  populationFill(parameters::instance().get_larva_eq(patchID),
                 minAge, distHold, eggVec);
  
  
  /*********/
  // larva
  /*********/
  minAge = parameters::instance().get_stage_time(0);
  maxAge = parameters::instance().get_stage_sum(1)-1;
  
  // set age distribution vector
  distHold.resize(maxAge+1 - minAge);
  std::copy(aDist.begin() + minAge, aDist.begin() + maxAge + 1, distHold.begin());
  
  // reserve estimated population size
  larvaVec.reserve(round(1.1*parameters::instance().get_larva_eq(patchID) * accumulate(distHold.begin(), distHold.end(), 0.0)));
  
  // fill initial population
  populationFill(parameters::instance().get_larva_eq(patchID),
                 minAge, distHold, larvaVec);
  
  
  /*********/
  // pupa
  /*********/
  minAge = parameters::instance().get_stage_sum(1);
  maxAge = parameters::instance().get_stage_sum(2) - 1;
  
  // set age distribution vector
  distHold.resize(maxAge+1 - minAge);
  std::copy(aDist.begin() + minAge, aDist.begin() + maxAge + 1, distHold.begin());
  
  // reserve estimated population size
  pupaVec.reserve(round(1.1*parameters::instance().get_larva_eq(patchID) * accumulate(distHold.begin(), distHold.end(), 0.0)));
  
  // fill initial population
  populationFill(parameters::instance().get_larva_eq(patchID),
                 minAge, distHold, pupaVec);
  
  
  /***********************************/
  // Solve adult distribution
  /***********************************/
  minAge = parameters::instance().get_stage_sum(2);
  
  // set age distribution vector
  aDist = markovDist(1.0 - parameters::instance().get_mu(3), parameters::instance().get_stage_time(3) * 3);
  
  // reserve estimated population size
  int popSize(round(1.1 * parameters::instance().get_adult_pop_eq(patchID) * accumulate(aDist.begin(), aDist.end(), 0.0)));
  aMaleVec.reserve(popSize);
  aFemaleVec.reserve(popSize);
  unFemaleVec.reserve(popSize);
  
  // fill initial population
  populationFill(parameters::instance().get_adult_pop_eq(patchID)/2,
                 minAge, aDist, aMaleVec);
  populationFill(parameters::instance().get_adult_pop_eq(patchID)/2,
                 minAge, aDist, unFemaleVec);
  
}


// This is also used in the multiLocusClass setup and reset functions
void CreateMosquitoes2Allele(const int& Leq, const int& minAge, const dVec& ageDist,
                             const Rcpp::ListOf<Rcpp::List>& aTypes, popVec& returnPop){
  
  // This only good for things that have 2 alleles per locus

  // number of loci is the number of sublists in aTypes
  size_t numLoci = aTypes.size();
  
  // holders for things in loop
  std::vector<double> probs;
  std::vector<std::string> alleles;
  std::vector<std::string> genHold(2);
  std::string genotype;
  int age;
  
  // loop over each age that mosquitoes can be
  for(size_t age = 0; age < ageDist.size(); age++){
    // loop over number of mosquitoes to create
    for(size_t num=0; num < round(Leq * ageDist[age]); num++){
      /*********/
      // set genotype
      /*********/
      // loop over each locus
      for(size_t loci=0; loci < numLoci; ++ loci){
        
        probs = Rcpp::as<std::vector<double> >( aTypes[loci]["probs"] );
        alleles = Rcpp::as<std::vector<std::string> >( aTypes[loci]["alleles"] );
        
        // generate each allele
        prng::instance().set_oneSample(probs);
        for(size_t gen=0; gen < 2; gen++){
          genHold[gen] = alleles[prng::instance().get_oneSample()];
        }
        
        // sort alleles
        std::sort(genHold.begin(), genHold.end());
        
        // combine alleles into final genotype
        for(auto const& s : genHold){ genotype += s; }
      } // end loop over loci
      
      
      // add new mosquito to population
      returnPop.emplace_back(Mosquito(minAge + age, genotype));
      
      // clear genotype for next one
      genotype.clear();
      
    } // end loop over number of mosquitoes
  } // end loop over age distribution

}



/**************************************
 * GENERATING FUNCTION
 **************************************/
void Daisy::DaisyOffspring(const std::string& fGen, const std::string& mGen){
  
  
  //// list of objects I create every time ////
  // vector<bool> fScore
  // vector<bool> mScore

  
  


  

  
  
  
  
  
  
  // get number of alleles
  // gets redefined every time
  // numAlleles = fGen.size()/2;
  
  /*****************************************************************************/
  // Score Each Allele
  /*****************************************************************************/
  std::vector<bool> fScore(numAlleles+1, false), mScore(numAlleles+1, false);
  
  //loop over loci, separate alleles and score
  index=1;
  for(size_t i=0; i<fGen.size(); i+=2, ++index){
    // female score
    if( (fGen[i] == 'H') || (fGen[i+1] == 'H')) {fScore[index] = true;}
    // male score
    if( (mGen[i] == 'H') || (mGen[i+1] == 'H')) {mScore[index] = true;}
    
  } // end scoring loop
  /*****************************************************************************/
  //End Split and Score
  /*****************************************************************************/
  
  
  //*****************************************************************************/
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
    if( (!fScore[i] && !fScore[i+1]) || (!fScore[i] && fScore[i+1]) ){
      //FF  or FT case
      // loop over alleles, all loci diploid
      for(size_t j=0; j<2; ++j){
        // fill allele and probs
        if(fGen[index+j] == 'W'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_mendelian_allele_begin(i,0),
                            reference::instance().get_mendelian_allele_end(i,0));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_mendelian_probs_begin(i,0),
                           reference::instance().get_mendelian_probs_end(i,0));
        } else if(fGen[index+j] == 'H'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_mendelian_allele_begin(i,1),
                            reference::instance().get_mendelian_allele_end(i,1));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_mendelian_probs_begin(i,1),
                           reference::instance().get_mendelian_probs_end(i,1));
        } else if(fGen[index+j] == 'R'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_mendelian_allele_begin(i,2),
                            reference::instance().get_mendelian_allele_end(i,2));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_mendelian_probs_begin(i,2),
                           reference::instance().get_mendelian_probs_end(i,2));
        } else if(fGen[index+j] == 'S'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_mendelian_allele_begin(i,3),
                            reference::instance().get_mendelian_allele_end(i,3));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_mendelian_probs_begin(i,3),
                           reference::instance().get_mendelian_probs_end(i,3));
        }
      } // end allele loop
      
    } else if(fScore[i] && !fScore[i+1]){
      //TF case
      // loop over alleles, all loci diploid
      for(size_t j=0; j<2; ++j){
        // fill allele and probs
        if(fGen[index+j] == 'W'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_cutting_allele_begin(i,0),
                            reference::instance().get_cutting_allele_end(i,0));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_cutting_probs_begin(i,0),
                           reference::instance().get_cutting_probs_end(i,0));
        } else if(fGen[index+j] == 'H'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_cutting_allele_begin(i,1),
                            reference::instance().get_cutting_allele_end(i,1));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_cutting_probs_begin(i,1),
                           reference::instance().get_cutting_probs_end(i,1));
        } else if(fGen[index+j] == 'R'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_cutting_allele_begin(i,2),
                            reference::instance().get_cutting_allele_end(i,2));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_cutting_probs_begin(i,2),
                           reference::instance().get_cutting_probs_end(i,2));
        } else if(fGen[index+j] == 'S'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_cutting_allele_begin(i,3),
                            reference::instance().get_cutting_allele_end(i,3));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_cutting_probs_begin(i,3),
                           reference::instance().get_cutting_probs_end(i,3));
        }
      } // end allele loop
      
    } else if(fScore[i] && fScore[i+1]){
      //TT case
      // loop over alleles, all loci diploid
      for(size_t j=0; j<2; ++j){
        // fill allele and probs
        if(fGen[index+j] == 'W'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_homing_allele_begin(i,0),
                            reference::instance().get_homing_allele_end(i,0));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_homing_probs_begin(i,0),
                           reference::instance().get_homing_probs_end(i,0));
        } else if(fGen[index+j] == 'H'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_homing_allele_begin(i,1),
                            reference::instance().get_homing_allele_end(i,1));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_homing_probs_begin(i,1),
                           reference::instance().get_homing_probs_end(i,1));
        } else if(fGen[index+j] == 'R'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_homing_allele_begin(i,2),
                            reference::instance().get_homing_allele_end(i,2));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_homing_probs_begin(i,2),
                           reference::instance().get_homing_probs_end(i,2));
        } else if(fGen[index+j] == 'S'){
          fAllele[i].insert(fAllele[i].end(), reference::instance().get_homing_allele_begin(i,3),
                            reference::instance().get_homing_allele_end(i,3));
          fProbs[i].insert(fProbs[i].end(), reference::instance().get_homing_probs_begin(i,3),
                           reference::instance().get_homing_probs_end(i,3));
        }
      } // end allele loop 
    } // end females
    
    // MALES
    if( (!mScore[i] && !mScore[i+1]) || (!mScore[i] && mScore[i+1]) ){
      //FF  or FT case
      // loop over alleles, all loci diploid
      for(size_t j=0; j<2; ++j){
        // fill allele and probs
        if(mGen[index+j] == 'W'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_mendelian_allele_begin(i,0),
                            reference::instance().get_mendelian_allele_end(i,0));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_mendelian_probs_begin(i,0),
                           reference::instance().get_mendelian_probs_end(i,0));
        } else if(mGen[index+j] == 'H'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_mendelian_allele_begin(i,1),
                            reference::instance().get_mendelian_allele_end(i,1));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_mendelian_probs_begin(i,1),
                           reference::instance().get_mendelian_probs_end(i,1));
        } else if(mGen[index+j] == 'R'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_mendelian_allele_begin(i,2),
                            reference::instance().get_mendelian_allele_end(i,2));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_mendelian_probs_begin(i,2),
                           reference::instance().get_mendelian_probs_end(i,2));
        } else if(mGen[index+j] == 'S'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_mendelian_allele_begin(i,3),
                            reference::instance().get_mendelian_allele_end(i,3));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_mendelian_probs_begin(i,3),
                           reference::instance().get_mendelian_probs_end(i,3));
        }
      } // end allele loop
      
    } else if(mScore[i] && !mScore[i+1]){
      //TF case
      // loop over alleles, all loci diploid
      for(size_t j=0; j<2; ++j){
        // fill allele and probs
        if(mGen[index+j] == 'W'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_cutting_allele_begin(i,0),
                            reference::instance().get_cutting_allele_end(i,0));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_cutting_probs_begin(i,0),
                           reference::instance().get_cutting_probs_end(i,0));
        } else if(mGen[index+j] == 'H'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_cutting_allele_begin(i,1),
                            reference::instance().get_cutting_allele_end(i,1));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_cutting_probs_begin(i,1),
                           reference::instance().get_cutting_probs_end(i,1));
        } else if(mGen[index+j] == 'R'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_cutting_allele_begin(i,2),
                            reference::instance().get_cutting_allele_end(i,2));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_cutting_probs_begin(i,2),
                           reference::instance().get_cutting_probs_end(i,2));
        } else if(mGen[index+j] == 'S'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_cutting_allele_begin(i,3),
                            reference::instance().get_cutting_allele_end(i,3));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_cutting_probs_begin(i,3),
                           reference::instance().get_cutting_probs_end(i,3));
        }
      } // end allele loop
      
      
    } else if(mScore[i] && mScore[i+1]){
      //TT case
      // loop over alleles, all loci diploid
      for(size_t j=0; j<2; ++j){
        // fill allele and probs
        if(mGen[index+j] == 'W'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_homing_allele_begin(i,0),
                            reference::instance().get_homing_allele_end(i,0));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_homing_probs_begin(i,0),
                           reference::instance().get_homing_probs_end(i,0));
        } else if(mGen[index+j] == 'H'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_homing_allele_begin(i,1),
                            reference::instance().get_homing_allele_end(i,1));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_homing_probs_begin(i,1),
                           reference::instance().get_homing_probs_end(i,1));
        } else if(mGen[index+j] == 'R'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_homing_allele_begin(i,2),
                            reference::instance().get_homing_allele_end(i,2));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_homing_probs_begin(i,2),
                           reference::instance().get_homing_probs_end(i,2));
        } else if(mGen[index+j] == 'S'){
          mAllele[i].insert(mAllele[i].end(), reference::instance().get_homing_allele_begin(i,3),
                            reference::instance().get_homing_allele_end(i,3));
          mProbs[i].insert(mProbs[i].end(), reference::instance().get_homing_probs_begin(i,3),
                           reference::instance().get_homing_probs_end(i,3));
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
  // value gets redefined, don't clear
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



