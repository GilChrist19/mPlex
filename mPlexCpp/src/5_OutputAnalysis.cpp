///////////////////////////////////////////////////////////////////////////////
//                           ____  __          ______          
//                ____ ___  / __ \/ /__  _  __/ ____/___  ____ 
//               / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
//              / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
//             /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/ 
//                                              /_/   /_/      
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Namespace and external functions
///////////////////////////////////////////////////////////////////////////////

#include "5_Reader.hpp"
#include <unordered_map>

#include <Rcpp.h>
using namespace Rcpp;

//This is to access the print function.
// It is not available outside of Rcpp!
extern "C" SEXP fwriteMain(SEXP MAT,   //matrix test
                          SEXP filename_Arg,
                          SEXP sep_Arg,
                          SEXP eol_Arg,
                          SEXP dec_Arg,
                          SEXP buffMB_Arg); // [1-1024] default 8MB

///////////////////////////////////////////////////////////////////////////////
// Simulation aggregation
///////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
void simAgg(const std::vector<std::vector<std::string> >& readFiles_,
            const std::vector<std::vector<std::string> >& writeFiles_,
            const std::string& largeFile_,
            const int& simTime_, const int& sampTime_, const int& maxRows_,
            const Rcpp::List& genKey_){
  
  ////////////////////
  // Hash Table
  ////////////////////
  /* This section builds the hash table that will be used to collapse individuals 
   * into genotypes of interest. The list used is generated in another function and 
   * printed as a csv.
   */
  
  // hash table
  std::unordered_map<std::string, int> genKey;
  int numGroups = Rcpp::as<Rcpp::IntegerVector>(genKey_[1]).size();
  std::unordered_map<std::string, int>::iterator findIt, endIt;
  
  // fill hash table
  for(size_t i=0; i<numGroups; ++i){
    genKey.insert(std::pair<std::string, int>(Rcpp::as<std::vector<std::string> >(genKey_[0])[i],
                                              Rcpp::as<std::vector<int> >(genKey_[1])[i]));
  }
  
  // get iterator to end
  numGroups = Rcpp::unique(Rcpp::as<Rcpp::IntegerVector>(genKey_[1])).size();
  endIt = genKey.end();
  
  
  ////////////////////
  // Initialize reading stuff
  ////////////////////
  MPLEXReader fileReader(largeFile_);
  std::vector<int> startVec({18,23}); // The first is male, second is female, corresponds to Time,Age,Genotype,Mate\n
  
  
  
  ////////////////////
  // Initialize data holders and return stuff
  ////////////////////
  // in data
  std::vector<SimpleMos> inData;
  inData.reserve(maxRows_*sizeof(SimpleMos)); // reserve space.
  
  // out data
  Rcpp::IntegerMatrix outData(simTime_/sampTime_ + 1, numGroups+2);
  Rcpp::CharacterVector cNames(Rcpp::unique(Rcpp::as<Rcpp::IntegerVector>(genKey_[1])).sort());
  cNames.push_front("Time");
  cNames.push_back("Other");
  Rcpp::colnames(outData) = cNames;
  
  // other stuff
  int oGroup(numGroups+1), group(0);
  Rcpp::IntegerVector tVec(seq(0,simTime_/sampTime_) * sampTime_ );
  
  
  ////////////////////
  // Loop over read file list
  ////////////////////
  numGroups = readFiles_[0].size();
  for(size_t sexInt : {0,1}){
    ////////////////////
    // Loop over male/female vector files
    ////////////////////
    for(size_t whichFile=0; whichFile < numGroups; ++whichFile){
      
      // reset data holders
      outData.fill(0);
      inData.clear();
      
      // read in file
      fileReader.readSimOut(readFiles_[sexInt][whichFile], startVec[sexInt], sexInt, inData);
      
      ////////////////////
      // Loop over elements from file
      ////////////////////
      for(const auto& it : inData){
        // set default group
        group = oGroup;
        // look to see if key is in map
        findIt = genKey.find(it.gen);
        // check that key is found
        if(findIt != endIt){
          group = findIt->second;
        }
        // store in out data
        outData(it.time / sampTime_, group) += 1;
      } // end loop over in data
      
      // fill in simTime
      outData(_,0)=tVec;
      
      // write to file
      fwriteMain(wrap(outData), wrap(writeFiles_[sexInt][whichFile]), wrap(","), wrap("\n"),
                 wrap("."), wrap(8));
      
    } // end loop over files in sex vector
  } // end loop over male/female file vectors
  
}










