#include <Rcpp.h>
using namespace Rcpp;

//' Daisy Drive Offspring
//'
//' Create a list of offspring genotypes and their probabilities
//'
//' @usage DaisyOffspring(fGen, mGen, reference)
//'
//' @param fGen Female genotype
//' @param mGen Male genotype
//' @param reference Offspring reference list
//'
//' @details Using the reference generated by \code{\link{MakeReference_DaisyDrive}},
//' this function expands the possible genotypes of the offspring and the ratios
//' that they occur. Similar to \code{\link{MultiplexOffspring_mLoci}} and
//' \code{\link{MultiplexOffspring_oLocus}}.
//'
//' @return List(Alleles, Probabilities)
//'
//' @examples
//' ref <- MakeReferenceDaisy(H = 0.9, R = 0, S = 0, d = .001)
//' fGen <- "WW"
//' mGen <- "WW"
//'
//' DaisyOffspring(fGen, mGen, ref)
//'
//' @export
// [[Rcpp::export]]
List DaisyOffspring_C(StringVector fGen, StringVector mGen, List& reference) {

  //Safety check. Maybe don't need.
  if(fGen[0].size() != mGen[0].size()){
    stop("Genotypes not the same length.");
  }

  //get number of alleles. Divide by two because diploid
  int numAlleles = fGen[0].size()/2;

  //setup matrix for individual parent alleles
  StringMatrix momAlleles(numAlleles, 2);
  StringMatrix dadAlleles(numAlleles, 2);

  //Scoring if they have homing allele
  LogicalVector fScore(numAlleles+1);
  LogicalVector mScore(numAlleles+1);

  /*****************************************************************************/
  //Split Loci Within Alleles and Score Each Allele
  /*****************************************************************************/

  //loop over loci, separate alleles and score
  for(int i=0; i<numAlleles; i++){
    //loop over alleles. All diploid
    for(int j=0; j<2; j++){
      momAlleles(i,j) = fGen[0][2*i+j];
      dadAlleles(i,j) = mGen[0][2*i+j];
    }//end locus loop

    //score if homing allele is present - female
    //Scoring: first score is always false. Last score is ignored.
    // The [] operators should allow ignoreing out-of-bounds on the final score.
    if(momAlleles(i,0) == 'H' | momAlleles(i,1) == 'H'){
      fScore[i+1]=1;
    }
    //score male
    if(dadAlleles(i,0) == 'H' | dadAlleles(i,1) == 'H'){
      fScore[i+1]=1;
    }
  }//end sublist and scoring loop over loci
  /*****************************************************************************/
  //End Split and Score
  /*****************************************************************************/

  /*****************************************************************************/
  //Determine Next-Gen alleles
  /*****************************************************************************/

  //setup male/female allele/probs lists
  List fAllele(numAlleles);
  List fProbs(numAlleles);
  List mAllele(numAlleles);
  List mProbs(numAlleles);

  //holder because Rcpp can't do lists of lists.
  List AlleleHold(2);
  List ProbsHold(2);


  //loop over all loci
  for(int i=0; i<numAlleles; i++){

    ////////////FEMALES////////////
    //3 if statements for 3 cases
    if((fScore[i]==FALSE && fScore[i+1]==FALSE) ||(fScore[i]==FALSE && fScore[i+1]==TRUE)){
      //FF  or FT case

      //loop over alleles. All diploid
      for(int j=0; j<2; j++){

        //Fill allele with letter and probs
        if(momAlleles(i,j) == 'W'){
          AlleleHold(j) = Rcpp::as<CharacterVector>(Rcpp::as<List>(Rcpp::as<List>(reference["mendelianAlleles"])[i])[0]);
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["W"];
        } else if(momAlleles(i,j) == 'H'){
          AlleleHold(j) = Rcpp::as<CharacterVector>(Rcpp::as<List>(Rcpp::as<List>(reference["mendelianAlleles"])[i])[1]);
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["H"];
        } else if(momAlleles(i,j) == 'R'){
          AlleleHold(j) = Rcpp::as<CharacterVector>(Rcpp::as<List>(Rcpp::as<List>(reference["mendelianAlleles"])[i])[2]);
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["R"];
        } else if(momAlleles(i,j) == 'S'){
          AlleleHold(j) = Rcpp::as<CharacterVector>(Rcpp::as<List>(Rcpp::as<List>(reference["mendelianAlleles"])[i])[3]);
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["S"];
        }

      }//end allele loop
    } else if(fScore[i]==TRUE && fScore[i+1]==FALSE){
      //TF case

      //loop over alleles at the locus. Everything is diploid here.
      for(int j=0; j<2; j++){

        //Fill allele with letter and probs
        if(momAlleles(i,j) == 'W'){
          AlleleHold(j) = Rcpp::as<CharacterVector>(Rcpp::as<List>(Rcpp::as<List>(reference["cuttingAlleles"])[i])[0]);
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["cutting"])[i])["W"];
        } else if(momAlleles(i,j) == 'H'){
          AlleleHold(j) = Rcpp::as<CharacterVector>(Rcpp::as<List>(Rcpp::as<List>(reference["cuttingAlleles"])[i])[1]);
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["cutting"])[i])["H"];
        } else if(momAlleles(i,j) == 'R'){
          AlleleHold(j) = Rcpp::as<CharacterVector>(Rcpp::as<List>(Rcpp::as<List>(reference["cuttingAlleles"])[i])[2]);
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["cutting"])[i])["R"];
        } else if(momAlleles(i,j) == 'S'){
          AlleleHold(j) = Rcpp::as<CharacterVector>(Rcpp::as<List>(Rcpp::as<List>(reference["cuttingAlleles"])[i])[3]);
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["cutting"])[i])["S"];
        }

      }//end allele loop
    } else if(fScore[i]==TRUE && fScore[i+1]==TRUE){
      //TT case

      //loop over alleles at the locus. Everything is diploid.
      for(int j=0; j<2; j++){

        //Fill allele with letter and probs
        if(momAlleles(i,j) == 'W'){
          AlleleHold(j) = Rcpp::as<CharacterVector>(Rcpp::as<List>(Rcpp::as<List>(reference["homingAlleles"])[i])[0]);
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["W"];
        } else if(momAlleles(i,j) == 'H'){
          AlleleHold(j) = Rcpp::as<CharacterVector>(Rcpp::as<List>(Rcpp::as<List>(reference["homingAlleles"])[i])[1]);
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["H"];
        } else if(momAlleles(i,j) == 'R'){
          AlleleHold(j) = Rcpp::as<CharacterVector>(Rcpp::as<List>(Rcpp::as<List>(reference["homingAlleles"])[i])[2]);
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["R"];
        } else if(momAlleles(i,j) == 'S'){
          AlleleHold(j) = Rcpp::as<CharacterVector>(Rcpp::as<List>(Rcpp::as<List>(reference["homingAlleles"])[i])[3]);
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["S"];
        }

      }//end allele loop
    }//end female if statement

    //Set things in list
    fAllele(i) = AlleleHold;
    fProbs(i) = ProbsHold;

    //Reset holder
    AlleleHold = List(2);
    ProbsHold = List(2);

    ////////////MALES////////////
    if((mScore[i]==FALSE && mScore[i+1]==FALSE) ||(mScore[i]==FALSE && mScore[i+1]==TRUE)){
      //FF or FT case

      //Loop over alleles at the locus. All diploid here.
      for(int j=0; j<2; j++){

        //Fill with allele and probs
        if(dadAlleles(i,j) == 'W'){
          AlleleHold(j) = Rcpp::as<CharacterVector>(Rcpp::as<List>(Rcpp::as<List>(reference["mendelianAlleles"])[i])[0]);
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["W"];
        } else if(dadAlleles(i,j) == 'H'){
          AlleleHold(j) = Rcpp::as<CharacterVector>(Rcpp::as<List>(Rcpp::as<List>(reference["mendelianAlleles"])[i])[1]);
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["H"];
        } else if(dadAlleles(i,j) == 'R'){
          AlleleHold(j) = Rcpp::as<CharacterVector>(Rcpp::as<List>(Rcpp::as<List>(reference["mendelianAlleles"])[i])[2]);
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["R"];
        } else if(dadAlleles(i,j) == 'S'){
          AlleleHold(j) = Rcpp::as<CharacterVector>(Rcpp::as<List>(Rcpp::as<List>(reference["mendelianAlleles"])[i])[3]);
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["S"];
        }

      }//end allele loop
    } else if(mScore[i]==TRUE && mScore[i+1]==FALSE){
      //TF case

      //loop over allele at this locus. All diploid
      for(int j=0; j<2; j++){

        //Fill with allele and probs
        if(dadAlleles(i,j) == 'W'){
          AlleleHold(j) = Rcpp::as<CharacterVector>(Rcpp::as<List>(Rcpp::as<List>(reference["cuttingAlleles"])[i])[0]);
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["cutting"])[i])["W"];
        } else if(dadAlleles(i,j) == 'H'){
          AlleleHold(j) = Rcpp::as<CharacterVector>(Rcpp::as<List>(Rcpp::as<List>(reference["cuttingAlleles"])[i])[1]);
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["cutting"])[i])["H"];
        } else if(dadAlleles(i,j) == 'R'){
          AlleleHold(j) = Rcpp::as<CharacterVector>(Rcpp::as<List>(Rcpp::as<List>(reference["cuttingAlleles"])[i])[2]);
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["cutting"])[i])["R"];
        } else if(dadAlleles(i,j) == 'S'){
          AlleleHold(j) = Rcpp::as<CharacterVector>(Rcpp::as<List>(Rcpp::as<List>(reference["cuttingAlleles"])[i])[3]);
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["cutting"])[i])["S"];
        }

      }//end allele loop
    } else if(mScore[i]==TRUE && mScore[i+1]==TRUE){
      //TT case

      //Loop over alleles at the locus. All diploid here.
      for(int j=0; j<2; j++){

        //Fill allele with letter and probs
        if(dadAlleles(i,j) == 'W'){
          AlleleHold(j) = Rcpp::as<CharacterVector>(Rcpp::as<List>(Rcpp::as<List>(reference["homingAlleles"])[i])[0]);
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["W"];
        } else if(dadAlleles(i,j) == 'H'){
          AlleleHold(j) = Rcpp::as<CharacterVector>(Rcpp::as<List>(Rcpp::as<List>(reference["homingAlleles"])[i])[1]);
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["H"];
        } else if(dadAlleles(i,j) == 'R'){
          AlleleHold(j) = Rcpp::as<CharacterVector>(Rcpp::as<List>(Rcpp::as<List>(reference["homingAlleles"])[i])[2]);
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["R"];
        } else if(dadAlleles(i,j) == 'S'){
          AlleleHold(j) = Rcpp::as<CharacterVector>(Rcpp::as<List>(Rcpp::as<List>(reference["homingAlleles"])[i])[3]);
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["S"];
        }

      }//end allele loop
    }//end male if statement

    //Set things in list
    mAllele(i) = AlleleHold;
    mProbs(i) = ProbsHold;

    //Reset holder
    AlleleHold = List(2);
    ProbsHold = List(2);

  }//end loci loop
  /*****************************************************************************/
  //End Alleles and Probs at Each Locus
  /*****************************************************************************/

  /*****************************************************************************/
  //Unlist Loci
  /*****************************************************************************/

  //unlist loop!
  for(int i=0; i<numAlleles; i++){

    int lenF = 0;
    int lenM = 0;

    //loop over depth to get total length of vector
    for(int j = 0; j<2; j++){
      lenF += Rcpp::as<StringVector>(Rcpp::as<List>(fAllele[i])[j]).length();
      lenM += Rcpp::as<StringVector>(Rcpp::as<List>(mAllele[i])[j]).length();
    }

    //set index to hold position and output vectors
    int fIndex = 0;
    StringVector fAHold(lenF);
    NumericVector fPHold(lenF);

    int mIndex = 0;
    StringVector mAHold(lenM);
    NumericVector mPHold(lenM);

    //loop over depth, concatenate alleles
    for(int j = 0; j<2; j++){

      //Holder for things
      StringVector holdAf = Rcpp::as<List>(fAllele[i])[j];
      NumericVector holdPf = Rcpp::as<List>(fProbs[i])[j];

      StringVector holdAm = Rcpp::as<List>(mAllele[i])[j];
      NumericVector holdPm = Rcpp::as<List>(mProbs[i])[j];

      //copy things into new outputs
      std::copy(holdAf.begin(), holdAf.end(), fAHold.begin()+fIndex);
      std::copy(holdPf.begin(), holdPf.end(), fPHold.begin()+fIndex);

      std::copy(holdAm.begin(), holdAm.end(), mAHold.begin()+mIndex);
      std::copy(holdPm.begin(), holdPm.end(), mPHold.begin()+mIndex);

      //increment index
      fIndex += holdAf.size();
      mIndex += holdAm.size();

    }//end allele loop

    //Set things in the list
    fAllele[i] = fAHold;
    fProbs[i] = fPHold;

    mAllele[i] = mAHold;
    mProbs[i] = mPHold;

  }//end loci loop
  /*****************************************************************************/
  //All Combinations of Male/Female for each Loci
  /*****************************************************************************/

  List lociAList(numAlleles);
  List lociPList(numAlleles);

  StringVector charHold(2);
  int lenF = 0;
  int lenM = 0;

  for(int i=0; i<numAlleles; i++){

    lenF = Rcpp::as<StringVector>(fAllele[i]).length();
    lenM = Rcpp::as<StringVector>(mAllele[i]).length();

    StringVector holdAll( lenF*lenM );
    NumericVector holdProb( lenF*lenM);

    //This is expand.grid(), then combining(pasting or multiplication), then sorting
    for(int y=0; y<lenF; y++){
      for(int z=0; z<lenM; z++){

        //combine and sort each possible allele at loci i
        charHold[0] = Rcpp::as<CharacterVector>(fAllele[i])[y];
        charHold[1] = Rcpp::as<CharacterVector>(mAllele[i])[z];
        holdAll[lenM*y+z] = collapse(charHold.sort());

        //combine probs via multiplication
        holdProb[lenM*y+z] = Rcpp::as<NumericVector>(fProbs[i])[y] * Rcpp::as<NumericVector>(mProbs[i])[z];

      }//end male choices loop
    }//end female choices loop

    ////////////AGGREGATE FUNCTION////////////
    //keep unique allele types
    lociAList[i] = unique(holdAll);
    LogicalVector matchedAlleles(holdAll.length());
    NumericVector combinedProbs(Rcpp::as<StringVector>(lociAList[i]).length());

    //This is the aggregate loop
    for(int y=0; y<Rcpp::as<StringVector>(lociAList[i]).length(); y++){

      //get locations of matching allele types
      for(int z=0; z<holdAll.length(); z++){
        matchedAlleles[z] = (holdAll[z] == Rcpp::as<StringVector>(lociAList[i])[y]);
      }

      //sum the value of each of them
      combinedProbs[y] = sum(as<NumericVector>(holdProb[matchedAlleles]));
    }//end aggregate loop

    lociPList[i] = combinedProbs;
    ////////////END AGGREGATE FUNCTION////////////

  }//end loop over alleles
  /*****************************************************************************/
  //End All Combinations of MAle/Female for Each Loci
  /*****************************************************************************/

  /*****************************************************************************/
  //Final Combinations of Loci for each Allele
  /*****************************************************************************/

  //final things in program
  //total length of output
  int depth = Rcpp::as<CharacterVector>(lociAList[0]).length();

  //vectors of cumulative and sublist lengths
  IntegerVector cumLen(numAlleles);
  IntegerVector subLen(numAlleles);

  //first sublist is length depth
  subLen[0] = depth;

  //get values from rest of list.
  for(int i=1; i<numAlleles; i++){
    cumLen[i] = depth;
    subLen[i] = Rcpp::as<CharacterVector>(lociAList[i]).length();
    depth *= subLen[i];
  }

  //setup output and holders for them.
  CharacterVector outAList(depth);
  CharacterVector oALHolder(depth);

  NumericVector outPList(depth);
  NumericVector oPLHolder(depth);

  //fill lists initially
  for(int i=0; i<subLen[0]; i++){
    //Alleles
    outAList[i] = Rcpp::as<CharacterVector>(lociAList[0])[i];
    oALHolder[i] = outAList[i];
    //Probs
    outPList[i] = Rcpp::as<NumericVector>(lociPList[0])[i];
    oPLHolder[i] = outPList[i];
  }


  for(int i=1; i<numAlleles; i++){
    for(int j=0; j<cumLen[i]; j++){
      for(int k=0; k<subLen[i]; k++){
        outAList[j*subLen[i]+k] = collapse(CharacterVector::create(oALHolder[j],Rcpp::as<CharacterVector>(lociAList[i])[k]));
        outPList[j*subLen[i]+k] = oPLHolder[j]*Rcpp::as<NumericVector>(lociPList[i])[k];
      }
    }

    //copy for next iteration
    std::copy(outAList.begin(), outAList.end(), oALHolder.begin());
    std::copy(outPList.begin(), outPList.end(), oPLHolder.begin());

  }//end shenanigans



  return List::create(
    _["Alleles"]  = outAList,
    _["Probabilities"] = outPList/sum(outPList)
  ) ;
}