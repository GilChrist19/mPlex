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
List oLocus(StringVector fGen, StringVector mGen, List& reference){

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
  bool fScore = false;
  bool mScore = false;

/*****************************************************************************/
//Split Loci Within Alleles and Score Each Allele
/*****************************************************************************/

  //loop over alleles. Everything diploid
  for(int j=0; j<2; j++){
    //loop over loci, variable number
    for(int i=0; i<numAlleles; i++){
      momAlleles(i,j) = fGen[0][numAlleles*j+i];
      dadAlleles(i,j) = mGen[0][numAlleles*j+i];

      //Score male/female at each locus. Just need one H anywhere.
      if(momAlleles(i,j) == 'H'){
        fScore = true;
      }
      if(dadAlleles(i,j) == 'H'){
        mScore = true;
      }

    }
  }
/*****************************************************************************/
//End Split and Score
/*****************************************************************************/

/*****************************************************************************/
//Determine Next-Gen alleles
/*****************************************************************************/

  //setup male/female allele/probs lists
  List fAllele(2);
  List fProbs(2);
  List mAllele(2);
  List mProbs(2);

  //holder because Rcpp can't do lists of lists.
  List AlleleHold(numAlleles);
  List ProbsHold(numAlleles);

  ////////////FEMALES////////////
  if(fScore){
    //TRUE - homing allele present
    //loop over alleles. All diploid
    for(int j=0; j<2; j++){
      //loop over loci.
      for(int i=0; i<numAlleles; i++){
        //Fill target with letter and probs
        if(momAlleles(i,j) == 'W'){
          AlleleHold(i) = StringVector::create("W", "H", "R", "S");
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["W"];
        } else if(momAlleles(i,j) == 'H'){
          AlleleHold(i) = StringVector::create("H", "S");
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["H"];
        } else if(momAlleles(i,j) == 'R'){
          AlleleHold(i) = "R";
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["R"];
        } else if(momAlleles(i,j) == 'S'){
          AlleleHold(i) = "S";
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["S"];
        }

      }//end loci loop

      //Set things in list
      fAllele(j) = AlleleHold;
      fProbs(j) = ProbsHold;

      //Reset holder
      AlleleHold = List(numAlleles);
      ProbsHold = List(numAlleles);

    }//end allele loop

  } else{
    //FALSE - no homing present
    //loop over alleles. All diploid
    for(int j=0; j<2; j++){
      //loop over loci.
      for(int i=0; i<numAlleles; i++){
        //Fill target with letter and probs
        if(momAlleles(i,j) == 'W'){
          AlleleHold(i) = StringVector::create("W", "S");
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["W"];
        } else if(momAlleles(i,j) == 'H'){
          AlleleHold(i) = StringVector::create("H", "S");
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["H"];
        } else if(momAlleles(i,j) == 'R'){
          AlleleHold(i) = "R";
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["R"];
        } else if(momAlleles(i,j) == 'S'){
          AlleleHold(i) = "S";
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["S"];
        }

      }//end loci loop

      //Set things in list
      fAllele(j) = AlleleHold;
      fProbs(j) = ProbsHold;

      //Reset holder
      AlleleHold = List(numAlleles);
      ProbsHold = List(numAlleles);

    }//end allele loop
  }//end female if statement



  ////////////MALES////////////
  if(mScore){
    //TRUE - homing allele present
    //loop over alleles. All diploid
    for(int j=0; j<2; j++){
      //loop over loci.
      for(int i=0; i<numAlleles; i++){
        //Fill target with letter and probs
        if(dadAlleles(i,j) == 'W'){
          AlleleHold(i) = StringVector::create("W", "H", "R", "S");
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["W"];
        } else if(dadAlleles(i,j) == 'H'){
          AlleleHold(i) = StringVector::create("H", "S");
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["H"];
        } else if(dadAlleles(i,j) == 'R'){
          AlleleHold(i) = "R";
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["R"];
        } else if(dadAlleles(i,j) == 'S'){
          AlleleHold(i) = "S";
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["S"];
        }

      }//end loci loop

      //Set things in list
      mAllele(j) = AlleleHold;
      mProbs(j) = ProbsHold;

      //Reset holder
      AlleleHold = List(numAlleles);
      ProbsHold = List(numAlleles);

    }//end allele loop

  } else{
    //FALSE - no homing present
    //loop over alleles. All diploid
    for(int j=0; j<2; j++){
      //loop over loci.
      for(int i=0; i<numAlleles; i++){
        //Fill target with letter and probs
        if(dadAlleles(i,j) == 'W'){
          AlleleHold(i) = StringVector::create("W", "S");
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["W"];
        } else if(dadAlleles(i,j) == 'H'){
          AlleleHold(i) = StringVector::create("H", "S");
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["H"];
        } else if(dadAlleles(i,j) == 'R'){
          AlleleHold(i) = "R";
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["R"];
        } else if(dadAlleles(i,j) == 'S'){
          AlleleHold(i) = "S";
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["S"];
        }

      }//end loci loop

      //Set things in list
      mAllele(j) = AlleleHold;
      mProbs(j) = ProbsHold;

      //Reset holder
      AlleleHold = List(numAlleles);
      ProbsHold = List(numAlleles);

    }//end allele loop
  }//end female if statement

/*****************************************************************************/
//End Alleles and Probs at Each Locus
/*****************************************************************************/

/*****************************************************************************/
//All Combinations of Loci for Each Allele
/*****************************************************************************/

  //setup return lists
  List fAllLoci(2);
  List fProbsLoci(2);
  List mAllLoci(2);
  List mProbsLoci(2);

  //Things used in loop
  //total length of output
  //vectors of cumulative and sublist lengths
  int depth = 0;
  IntegerVector cumLen(numAlleles);
  IntegerVector subLen(numAlleles);

  for(int numA=0; numA<2; numA++){

    ////////////FEMALES////////////
    //total length of output
    depth = Rcpp::as<CharacterVector>(Rcpp::as<List>(fAllele[numA])[0]).length();

    //first sublist is length depth
    subLen[0] = depth;

    //get values from rest of list.
    for(int i=1; i<numAlleles; i++){
      cumLen[i] = depth;
      subLen[i] = Rcpp::as<CharacterVector>(Rcpp::as<List>(fAllele[numA])[i]).length();
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
      outAList[i] = Rcpp::as<CharacterVector>(Rcpp::as<List>(fAllele[numA])[0])[i];
      oALHolder[i] = outAList[i];
      //Probs
      outPList[i] = Rcpp::as<NumericVector>(Rcpp::as<List>(fProbs[numA])[0])[i];
      oPLHolder[i] = outPList[i];
    }

    for(int i=1; i<numAlleles; i++){
      for(int j=0; j<cumLen[i]; j++){
        for(int k=0; k<subLen[i]; k++){
          outAList[j*subLen[i]+k] = collapse(CharacterVector::create(oALHolder[j],
                                                Rcpp::as<CharacterVector>(Rcpp::as<List>(fAllele[numA])[i])[k]));
          outPList[j*subLen[i]+k] = oPLHolder[j]*Rcpp::as<NumericVector>(Rcpp::as<List>(fProbs[numA])[i])[k];
        }
      }

      //copy for next iteration
      std::copy(outAList.begin(), outAList.end(), oALHolder.begin());
      std::copy(outPList.begin(), outPList.end(), oPLHolder.begin());
    }//end loop loci

    fAllLoci[numA] = outAList;
    fProbsLoci[numA] = outPList;



    ////////////MALES////////////
    //total length of output
    depth = Rcpp::as<CharacterVector>(Rcpp::as<List>(mAllele[numA])[0]).length();

    //first sublist is length depth
    subLen[0] = depth;

    //get values from rest of list.
    for(int i=1; i<numAlleles; i++){
      cumLen[i] = depth;
      subLen[i] = Rcpp::as<CharacterVector>(Rcpp::as<List>(mAllele[numA])[i]).length();
      depth *= subLen[i];
    }

    //setup output and holders for them.
    outAList = CharacterVector(depth);
    oALHolder = CharacterVector(depth);

    outPList = NumericVector(depth);
    oPLHolder = NumericVector(depth);

    //fill lists initially
    for(int i=0; i<subLen[0]; i++){
      //Alleles
      outAList[i] = Rcpp::as<CharacterVector>(Rcpp::as<List>(mAllele[numA])[0])[i];
      oALHolder[i] = outAList[i];
      //Probs
      outPList[i] = Rcpp::as<NumericVector>(Rcpp::as<List>(mProbs[numA])[0])[i];
      oPLHolder[i] = outPList[i];
    }

    for(int i=1; i<numAlleles; i++){
      for(int j=0; j<cumLen[i]; j++){
        for(int k=0; k<subLen[i]; k++){
          outAList[j*subLen[i]+k] = collapse(CharacterVector::create(oALHolder[j],
                                                                     Rcpp::as<CharacterVector>(Rcpp::as<List>(mAllele[numA])[i])[k]));
          outPList[j*subLen[i]+k] = oPLHolder[j]*Rcpp::as<NumericVector>(Rcpp::as<List>(mProbs[numA])[i])[k];
        }
      }

      //copy for next iteration
      std::copy(outAList.begin(), outAList.end(), oALHolder.begin());
      std::copy(outPList.begin(), outPList.end(), oPLHolder.begin());
    }//end loop loci

    mAllLoci[numA] = outAList;
    mProbsLoci[numA] = outPList;

  }//end loop over 2 alleles
/*****************************************************************************/
//End All Combinations of Loci for Each Allele
/*****************************************************************************/

/*****************************************************************************/
//Unlist Loci
/*****************************************************************************/

  //get lengths for things
  subLen[0] = Rcpp::as<CharacterVector>(fAllLoci[0]).length();
  subLen[1] = Rcpp::as<CharacterVector>(fAllLoci[1]).length();

  cumLen[0] = Rcpp::as<CharacterVector>(mAllLoci[0]).length();
  cumLen[1] = Rcpp::as<CharacterVector>(mAllLoci[1]).length();

  //create new vectors for things
  CharacterVector fLociAll = CharacterVector(subLen[0]+subLen[1]);
  NumericVector fProbsAll = NumericVector(subLen[0]+subLen[1]);

  CharacterVector mLociAll = CharacterVector(cumLen[0]+cumLen[1]);
  NumericVector mProbsAll = NumericVector(cumLen[0]+cumLen[1]);

  //copy things to unlist
  std::copy(Rcpp::as<CharacterVector>(fAllLoci[0]).begin(),
            Rcpp::as<CharacterVector>(fAllLoci[0]).end(),
            fLociAll.begin());
  std::copy(Rcpp::as<CharacterVector>(fAllLoci[1]).begin(),
            Rcpp::as<CharacterVector>(fAllLoci[1]).end(),
            fLociAll.begin()+subLen[0]);
  std::copy(Rcpp::as<NumericVector>(fProbsLoci[0]).begin(),
            Rcpp::as<NumericVector>(fProbsLoci[0]).end(),
            fProbsAll.begin());
  std::copy(Rcpp::as<NumericVector>(fProbsLoci[1]).begin(),
            Rcpp::as<NumericVector>(fProbsLoci[1]).end(),
            fProbsAll.begin()+subLen[0]);

  std::copy(Rcpp::as<CharacterVector>(mAllLoci[0]).begin(),
            Rcpp::as<CharacterVector>(mAllLoci[0]).end(),
            mLociAll.begin());
  std::copy(Rcpp::as<CharacterVector>(mAllLoci[1]).begin(),
            Rcpp::as<CharacterVector>(mAllLoci[1]).end(),
            mLociAll.begin()+cumLen[0]);
  std::copy(Rcpp::as<NumericVector>(mProbsLoci[0]).begin(),
            Rcpp::as<NumericVector>(mProbsLoci[0]).end(),
            mProbsAll.begin());
  std::copy(Rcpp::as<NumericVector>(mProbsLoci[1]).begin(),
            Rcpp::as<NumericVector>(mProbsLoci[1]).end(),
            mProbsAll.begin()+cumLen[0]);

/*****************************************************************************/
//End Unlist
/*****************************************************************************/

/*****************************************************************************/
//All combinations of genotypes and probabilities
/*****************************************************************************/

  //setup output
  depth = fLociAll.length()*mLociAll.length();
  CharacterVector outAList(depth);
  NumericVector outPList(depth);


  depth = mLociAll.length();


  for(int i=0; i<fLociAll.length(); i++){
    for(int j=0; j<depth; j++){

      //Female-Male allele/probs
      outAList[i*depth+j] = collapse(CharacterVector::create(fLociAll[i], mLociAll[j]));
      outPList[i*depth+j] = fProbsAll[i]*mProbsAll[j];

    }


  }

/*****************************************************************************/
//End All Combinations
/*****************************************************************************/

/*****************************************************************************/
//Aggregate and Return
/*****************************************************************************/

  //Get unique genotypes
  fLociAll = unique(outAList);

  //set values for use in loop or return
  depth = fLociAll.length();
  LogicalVector Matches(outAList.length());
  //LogicalVector NotZeroMatches(depth);
  fProbsAll = NumericVector(depth);

  //Aggregate loop
  for(int i=0; i<depth; i++){

    //Get locations of matching alleles
    for(int j=0; j<outAList.length(); j++){
      Matches[j] = (fLociAll[i]==outAList[j]);
    }

    //Sum probabilities of matching alleles
    fProbsAll[i] = sum(Rcpp::as<NumericVector>(outPList[Matches]));

    //Keep track of what's not zero
    //if(fProbsAll[i]!=0){NotZeroMatches[i]=TRUE;}

  }//End Aggregate Loop



    return List::create(
      _["Alleles"]  = fLociAll,
      _["Probabilities"] = fProbsAll/sum(fProbsAll)//Rcpp::as<NumericVector>(fProbsAll[NotZeroMatches])/sum(fProbsAll)
    ) ;
}


// [[Rcpp::export]]
List oLocus_nonZero(StringVector fGen, StringVector mGen, List& reference){

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
  bool fScore = false;
  bool mScore = false;

  /*****************************************************************************/
  //Split Loci Within Alleles and Score Each Allele
  /*****************************************************************************/

  //loop over alleles. Everything diploid
  for(int j=0; j<2; j++){
    //loop over loci, variable number
    for(int i=0; i<numAlleles; i++){
      momAlleles(i,j) = fGen[0][numAlleles*j+i];
      dadAlleles(i,j) = mGen[0][numAlleles*j+i];

      //Score male/female at each locus. Just need one H anywhere.
      if(momAlleles(i,j) == 'H'){
        fScore = true;
      }
      if(dadAlleles(i,j) == 'H'){
        mScore = true;
      }

    }
  }
  /*****************************************************************************/
  //End Split and Score
  /*****************************************************************************/

  /*****************************************************************************/
  //Determine Next-Gen alleles
  /*****************************************************************************/

  //setup male/female allele/probs lists
  List fAllele(2);
  List fProbs(2);
  List mAllele(2);
  List mProbs(2);

  //holder because Rcpp can't do lists of lists.
  List AlleleHold(numAlleles);
  List ProbsHold(numAlleles);

  ////////////FEMALES////////////
  if(fScore){
    //TRUE - homing allele present
    //loop over alleles. All diploid
    for(int j=0; j<2; j++){
      //loop over loci.
      for(int i=0; i<numAlleles; i++){
        //Fill target with letter and probs
        if(momAlleles(i,j) == 'W'){
          AlleleHold(i) = StringVector::create("W", "H", "R", "S");
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["W"];
        } else if(momAlleles(i,j) == 'H'){
          AlleleHold(i) = StringVector::create("H", "S");
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["H"];
        } else if(momAlleles(i,j) == 'R'){
          AlleleHold(i) = "R";
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["R"];
        } else if(momAlleles(i,j) == 'S'){
          AlleleHold(i) = "S";
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["S"];
        }

      }//end loci loop

      //Set things in list
      fAllele(j) = AlleleHold;
      fProbs(j) = ProbsHold;

      //Reset holder
      AlleleHold = List(numAlleles);
      ProbsHold = List(numAlleles);

    }//end allele loop

  } else{
    //FALSE - no homing present
    //loop over alleles. All diploid
    for(int j=0; j<2; j++){
      //loop over loci.
      for(int i=0; i<numAlleles; i++){
        //Fill target with letter and probs
        if(momAlleles(i,j) == 'W'){
          AlleleHold(i) = StringVector::create("W", "S");
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["W"];
        } else if(momAlleles(i,j) == 'H'){
          AlleleHold(i) = StringVector::create("H", "S");
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["H"];
        } else if(momAlleles(i,j) == 'R'){
          AlleleHold(i) = "R";
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["R"];
        } else if(momAlleles(i,j) == 'S'){
          AlleleHold(i) = "S";
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["S"];
        }

      }//end loci loop

      //Set things in list
      fAllele(j) = AlleleHold;
      fProbs(j) = ProbsHold;

      //Reset holder
      AlleleHold = List(numAlleles);
      ProbsHold = List(numAlleles);

    }//end allele loop
  }//end female if statement



  ////////////MALES////////////
  if(mScore){
    //TRUE - homing allele present
    //loop over alleles. All diploid
    for(int j=0; j<2; j++){
      //loop over loci.
      for(int i=0; i<numAlleles; i++){
        //Fill target with letter and probs
        if(dadAlleles(i,j) == 'W'){
          AlleleHold(i) = StringVector::create("W", "H", "R", "S");
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["W"];
        } else if(dadAlleles(i,j) == 'H'){
          AlleleHold(i) = StringVector::create("H", "S");
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["H"];
        } else if(dadAlleles(i,j) == 'R'){
          AlleleHold(i) = "R";
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["R"];
        } else if(dadAlleles(i,j) == 'S'){
          AlleleHold(i) = "S";
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["homing"])[i])["S"];
        }

      }//end loci loop

      //Set things in list
      mAllele(j) = AlleleHold;
      mProbs(j) = ProbsHold;

      //Reset holder
      AlleleHold = List(numAlleles);
      ProbsHold = List(numAlleles);

    }//end allele loop

  } else{
    //FALSE - no homing present
    //loop over alleles. All diploid
    for(int j=0; j<2; j++){
      //loop over loci.
      for(int i=0; i<numAlleles; i++){
        //Fill target with letter and probs
        if(dadAlleles(i,j) == 'W'){
          AlleleHold(i) = StringVector::create("W", "S");
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["W"];
        } else if(dadAlleles(i,j) == 'H'){
          AlleleHold(i) = StringVector::create("H", "S");
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["H"];
        } else if(dadAlleles(i,j) == 'R'){
          AlleleHold(i) = "R";
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["R"];
        } else if(dadAlleles(i,j) == 'S'){
          AlleleHold(i) = "S";
          ProbsHold(i) = Rcpp::as<List>(Rcpp::as<List>(reference["mendelian"])[i])["S"];
        }

      }//end loci loop

      //Set things in list
      mAllele(j) = AlleleHold;
      mProbs(j) = ProbsHold;

      //Reset holder
      AlleleHold = List(numAlleles);
      ProbsHold = List(numAlleles);

    }//end allele loop
  }//end female if statement

  /*****************************************************************************/
  //End Alleles and Probs at Each Locus
  /*****************************************************************************/

  /*****************************************************************************/
  //All Combinations of Loci for Each Allele
  /*****************************************************************************/

  //setup return lists
  List fAllLoci(2);
  List fProbsLoci(2);
  List mAllLoci(2);
  List mProbsLoci(2);

  //Things used in loop
  //total length of output
  //vectors of cumulative and sublist lengths
  int depth = 0;
  IntegerVector cumLen(numAlleles);
  IntegerVector subLen(numAlleles);

  for(int numA=0; numA<2; numA++){

    ////////////FEMALES////////////
    //total length of output
    depth = Rcpp::as<CharacterVector>(Rcpp::as<List>(fAllele[numA])[0]).length();

    //first sublist is length depth
    subLen[0] = depth;

    //get values from rest of list.
    for(int i=1; i<numAlleles; i++){
      cumLen[i] = depth;
      subLen[i] = Rcpp::as<CharacterVector>(Rcpp::as<List>(fAllele[numA])[i]).length();
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
      outAList[i] = Rcpp::as<CharacterVector>(Rcpp::as<List>(fAllele[numA])[0])[i];
      oALHolder[i] = outAList[i];
      //Probs
      outPList[i] = Rcpp::as<NumericVector>(Rcpp::as<List>(fProbs[numA])[0])[i];
      oPLHolder[i] = outPList[i];
    }

    for(int i=1; i<numAlleles; i++){
      for(int j=0; j<cumLen[i]; j++){
        for(int k=0; k<subLen[i]; k++){
          outAList[j*subLen[i]+k] = collapse(CharacterVector::create(oALHolder[j],
                                                                     Rcpp::as<CharacterVector>(Rcpp::as<List>(fAllele[numA])[i])[k]));
          outPList[j*subLen[i]+k] = oPLHolder[j]*Rcpp::as<NumericVector>(Rcpp::as<List>(fProbs[numA])[i])[k];
        }
      }

      //copy for next iteration
      std::copy(outAList.begin(), outAList.end(), oALHolder.begin());
      std::copy(outPList.begin(), outPList.end(), oPLHolder.begin());
    }//end loop loci

    fAllLoci[numA] = outAList;
    fProbsLoci[numA] = outPList;



    ////////////MALES////////////
    //total length of output
    depth = Rcpp::as<CharacterVector>(Rcpp::as<List>(mAllele[numA])[0]).length();

    //first sublist is length depth
    subLen[0] = depth;

    //get values from rest of list.
    for(int i=1; i<numAlleles; i++){
      cumLen[i] = depth;
      subLen[i] = Rcpp::as<CharacterVector>(Rcpp::as<List>(mAllele[numA])[i]).length();
      depth *= subLen[i];
    }

    //setup output and holders for them.
    outAList = CharacterVector(depth);
    oALHolder = CharacterVector(depth);

    outPList = NumericVector(depth);
    oPLHolder = NumericVector(depth);

    //fill lists initially
    for(int i=0; i<subLen[0]; i++){
      //Alleles
      outAList[i] = Rcpp::as<CharacterVector>(Rcpp::as<List>(mAllele[numA])[0])[i];
      oALHolder[i] = outAList[i];
      //Probs
      outPList[i] = Rcpp::as<NumericVector>(Rcpp::as<List>(mProbs[numA])[0])[i];
      oPLHolder[i] = outPList[i];
    }

    for(int i=1; i<numAlleles; i++){
      for(int j=0; j<cumLen[i]; j++){
        for(int k=0; k<subLen[i]; k++){
          outAList[j*subLen[i]+k] = collapse(CharacterVector::create(oALHolder[j],
                                                                     Rcpp::as<CharacterVector>(Rcpp::as<List>(mAllele[numA])[i])[k]));
          outPList[j*subLen[i]+k] = oPLHolder[j]*Rcpp::as<NumericVector>(Rcpp::as<List>(mProbs[numA])[i])[k];
        }
      }

      //copy for next iteration
      std::copy(outAList.begin(), outAList.end(), oALHolder.begin());
      std::copy(outPList.begin(), outPList.end(), oPLHolder.begin());
    }//end loop loci

    mAllLoci[numA] = outAList;
    mProbsLoci[numA] = outPList;

  }//end loop over 2 alleles
  /*****************************************************************************/
  //End All Combinations of Loci for Each Allele
  /*****************************************************************************/

  /*****************************************************************************/
  //Unlist Loci
  /*****************************************************************************/

  //get lengths for things
  subLen[0] = Rcpp::as<CharacterVector>(fAllLoci[0]).length();
  subLen[1] = Rcpp::as<CharacterVector>(fAllLoci[1]).length();

  cumLen[0] = Rcpp::as<CharacterVector>(mAllLoci[0]).length();
  cumLen[1] = Rcpp::as<CharacterVector>(mAllLoci[1]).length();

  //create new vectors for things
  CharacterVector fLociAll = CharacterVector(subLen[0]+subLen[1]);
  NumericVector fProbsAll = NumericVector(subLen[0]+subLen[1]);

  CharacterVector mLociAll = CharacterVector(cumLen[0]+cumLen[1]);
  NumericVector mProbsAll = NumericVector(cumLen[0]+cumLen[1]);

  //copy things to unlist
  std::copy(Rcpp::as<CharacterVector>(fAllLoci[0]).begin(),
            Rcpp::as<CharacterVector>(fAllLoci[0]).end(),
            fLociAll.begin());
  std::copy(Rcpp::as<CharacterVector>(fAllLoci[1]).begin(),
            Rcpp::as<CharacterVector>(fAllLoci[1]).end(),
            fLociAll.begin()+subLen[0]);
  std::copy(Rcpp::as<NumericVector>(fProbsLoci[0]).begin(),
            Rcpp::as<NumericVector>(fProbsLoci[0]).end(),
            fProbsAll.begin());
  std::copy(Rcpp::as<NumericVector>(fProbsLoci[1]).begin(),
            Rcpp::as<NumericVector>(fProbsLoci[1]).end(),
            fProbsAll.begin()+subLen[0]);

  std::copy(Rcpp::as<CharacterVector>(mAllLoci[0]).begin(),
            Rcpp::as<CharacterVector>(mAllLoci[0]).end(),
            mLociAll.begin());
  std::copy(Rcpp::as<CharacterVector>(mAllLoci[1]).begin(),
            Rcpp::as<CharacterVector>(mAllLoci[1]).end(),
            mLociAll.begin()+cumLen[0]);
  std::copy(Rcpp::as<NumericVector>(mProbsLoci[0]).begin(),
            Rcpp::as<NumericVector>(mProbsLoci[0]).end(),
            mProbsAll.begin());
  std::copy(Rcpp::as<NumericVector>(mProbsLoci[1]).begin(),
            Rcpp::as<NumericVector>(mProbsLoci[1]).end(),
            mProbsAll.begin()+cumLen[0]);

  /*****************************************************************************/
  //End Unlist
  /*****************************************************************************/

  /*****************************************************************************/
  //All combinations of genotypes and probabilities
  /*****************************************************************************/

  //setup output
  depth = fLociAll.length()*mLociAll.length();
  CharacterVector outAList(depth);
  NumericVector outPList(depth);


  depth = mLociAll.length();


  for(int i=0; i<fLociAll.length(); i++){
    for(int j=0; j<depth; j++){

      //Female-Male allele/probs
      outAList[i*depth+j] = collapse(CharacterVector::create(fLociAll[i], mLociAll[j]));
      outPList[i*depth+j] = fProbsAll[i]*mProbsAll[j];

    }


  }

  /*****************************************************************************/
  //End All Combinations
  /*****************************************************************************/

  /*****************************************************************************/
  //Aggregate and Return
  /*****************************************************************************/

  //Get unique genotypes
  fLociAll = unique(outAList);

  //set values for use in loop or return
  depth = fLociAll.length();
  LogicalVector Matches(outAList.length());
  LogicalVector NotZeroMatches(depth);
  fProbsAll = NumericVector(depth);

  //Aggregate loop
  for(int i=0; i<depth; i++){

    //Get locations of matching alleles
    for(int j=0; j<outAList.length(); j++){
      Matches[j] = (fLociAll[i]==outAList[j]);
    }

    //Sum probabilities of matching alleles
    fProbsAll[i] = sum(Rcpp::as<NumericVector>(outPList[Matches]));

    //Keep track of what's not zero
    if(fProbsAll[i]!=0){NotZeroMatches[i]=TRUE;}

  }//End Aggregate Loop



  return List::create(
    _["Alleles"]  = fLociAll[NotZeroMatches],
    _["Probabilities"] = Rcpp::as<NumericVector>(fProbsAll[NotZeroMatches])/sum(fProbsAll)
  ) ;
}




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



























