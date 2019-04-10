#include "5_Reader.hpp"


#include <Rcpp.h>
using namespace Rcpp;






// [[Rcpp::export]]
void testRead(std::string& female_, std::string& male_, const int& size_){
  
  Rcpp::Rcout << "In testRead\n";
  
  // initialize reader
  MPLEXReader fileReader(female_);
  
  Rcpp::Rcout<< "Initialized reader\n";
  
  // initialize vector
  std::vector<SimpleMos> inData;
  
  
  // read male
  fileReader.readSimOut(male_, 18, 0, inData);
  
  // print male to screen
  Rcpp::Rcout << "MALES\n";
  for(auto it : inData){
    Rcpp::Rcout << "Time: " << it.time << "\t Gen: " << it.gen << std::endl; 
  }
  
  
  
  //read female
  inData.clear();
  fileReader.readSimOut(female_, 23, 1, inData);
  
  
  // print female to screen
  Rcpp::Rcout << "FEMALES\n";
  for(auto it : inData){
    Rcpp::Rcout << "Time: " << it.time << "\t Gen: " << it.gen << std::endl; 
  }
  
  
}

