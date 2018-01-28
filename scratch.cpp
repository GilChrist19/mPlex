#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
List TEST(StringVector fGen, StringVector mGen, List& reference) {

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
  LogicalVector fScore(numAlleles);
  LogicalVector mScore(numAlleles);

  //loop over loci, separate alleles and score
  for(int i=0; i<numAlleles; i++){
    //loop over alleles. All diploid
    for(int j=0; j<2; j++){
      momAlleles(i,j) = fGen[0][2*i+j];
      dadAlleles(i,j) = mGen[0][2*i+j];
    }//end locus loop

    //score if homing allele is present - female
    if(momAlleles(i,0) == 'H' | momAlleles(i,1) == 'H'){
      fScore[i]=1;
    }
    //score male
    if(dadAlleles(i,0) == 'H' | dadAlleles(i,1) == 'H'){
      fScore[i]=1;
    }
  }//end sublist and scoring loop over loci

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

    //Females
    if(fScore[i]){
      //if homing
      //loop over alleles. All diploid
      for(int j=0; j<2; j++){

        //Fill allele with letter and probs
        if(momAlleles(i,j) == 'W'){
          AlleleHold(j) = StringVector::create("W", "H", "R", "S");
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["W"];
        } else if(momAlleles(i,j) == 'H'){
          AlleleHold(j) = StringVector::create("H", "S");
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["H"];
        } else if(momAlleles(i,j) == 'R'){
          AlleleHold(j) = "R";
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["R"];
        } else if(momAlleles(i,j) == 'S'){
          AlleleHold(j) = "S";
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["S"];
        }

      }//end allele loop
    } else{
      //If there is not homing
      //loop over alleles at the locus. Everything is diploid here.
      for(int j=0; j<2; j++){

        //Fill allele with letter and probs
        if(momAlleles(i,j) == 'W'){
          AlleleHold(j) = StringVector::create("W", "S");
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["W"];
        } else if(momAlleles(i,j) == 'H'){
          AlleleHold(j) = StringVector::create("H", "S");
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["H"];
        } else if(momAlleles(i,j) == 'R'){
          AlleleHold(j) = "R";
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["R"];
        } else if(momAlleles(i,j) == 'S'){
          AlleleHold(j) = "S";
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["S"];
        }

      }//end allele loop
    }//end female if statement

    //Set things in list
    fAllele(i) = AlleleHold;
    fProbs(i) = ProbsHold;

    //Reset holder
    AlleleHold = List(2);
    ProbsHold = List(2);

    //Males!
    if(mScore[i]){
      //If there is homing
      //Loop over alleles at the locus. All diploid here.
      for(int j=0; j<2; j++){

        //Fill allele with letter and probs
        if(dadAlleles(i,j) == 'W'){
          AlleleHold(j) = StringVector::create("W", "H", "R", "S");
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["W"];
        } else if(dadAlleles(i,j) == 'H'){
          AlleleHold(j) = StringVector::create("H", "S");
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["H"];
        } else if(dadAlleles(i,j) == 'R'){
          AlleleHold(j) = "R";
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["R"];
        } else if(dadAlleles(i,j) == 'S'){
          AlleleHold(j) = "S";
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["S"];
        }

      }//end allele loop
    } else{
      //if no homing
      //loop over allele at this locus. All diploid
      for(int j=0; j<2; j++){

        //Fill with allele and probs
        if(dadAlleles(i,j) == 'W'){
          AlleleHold(j) = StringVector::create("W", "S");
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["W"];
        } else if(dadAlleles(i,j) == 'H'){
          AlleleHold(j) = StringVector::create("H", "S");
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["H"];
        } else if(dadAlleles(i,j) == 'R'){
          AlleleHold(j) = "R";
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["R"];
        } else if(dadAlleles(i,j) == 'S'){
          AlleleHold(j) = "S";
          ProbsHold(j) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["S"];
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

    //FEMALES
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

    //BEGIN AGGREGATE FUNCTION
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
    //END AGGREGATE FUNCTION

  }//end loop over alleles



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


//This is a version of expand.grid()

// [[Rcpp::export]]
CharacterVector ListTest(List myList){

  //Length of total list
  int myListLen = myList.length();
  //Totel length of multiplied sublists
  int depth = Rcpp::as<CharacterVector>(myList[0]).length();

  //vectors of cumulative and sublist lengths
  IntegerVector cumLen(myListLen);
  IntegerVector subLen(myListLen);

  //first sublist is length depth
  subLen[0] = depth;

  //get values from rest of list.
  for(int i=1; i<myList.length(); i++){
    cumLen[i] = depth;
    subLen[i] = Rcpp::as<CharacterVector>(myList[i]).length();
    depth *= subLen[i];
  }

  CharacterVector holder(depth);
  CharacterVector output(depth);


  //fill lists initially
  for(int i=0; i<subLen[0]; i++){
    output[i] = Rcpp::as<CharacterVector>(myList[0])[i];
    holder[i] = output[i];
  }



  for(int i=1; i<myListLen; i++){
    for(int j=0; j<cumLen[i]; j++){
      for(int k=0; k<subLen[i]; k++){
        output[j*subLen[i]+k] = collapse(CharacterVector::create(holder[j],Rcpp::as<CharacterVector>(myList[i])[k]));
      }

    }

    std::copy(output.begin(), output.end(), holder.begin());

  }//end shenanigans





  return(output);
}



























