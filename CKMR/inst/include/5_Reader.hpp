///////////////////////////////////////////////////////////////////////////////
//     ________ __ __  _______ 
//    / ____/ //_//  |/  / __ \
//   / /   / ,<  / /|_/ / /_/ /
//  / /___/ /| |/ /  / / _, _/ 
//  \____/_/ |_/_/  /_/_/ |_|  
//   
///////////////////////////////////////////////////////////////////////////////

#ifndef CKMR_READ
#define CKMR_READ

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

//#include <Rcpp.h> // because used with armadillo later. can swap without armadillo
//#include <RcppArmadillo.h>

///////////////////////////////////////////////////////////////////////////////
//Struct Definition
///////////////////////////////////////////////////////////////////////////////

struct SimpleMos{
  int time;
  std::string gen;
  
  SimpleMos(int time_, std::string gen_) : time(time_), gen(gen_){}
};

///////////////////////////////////////////////////////////////////////////////
// Class Definition
///////////////////////////////////////////////////////////////////////////////

class MPLEXReader {
public:
  
  // constructor & destructor
  MPLEXReader(const std::string& _bigFile){
    
    // set size for holdString, use largest file (determined and given from elsewhere)
    holdStream.open(_bigFile);  // open stream
    holdStream.seekg(0, std::ios::end);     // skip to end of file
    fileSize = holdStream.tellg();          // get size of file
    holdString.reserve(fileSize);           // reserve space
    holdStream.close();
    
    // initialize holders for conversions
    baseDBL = 0.0;
    expDBL = 0.0;
    baseINT = 0;
    
  };
  
  ~MPLEXReader(){};
  
  // main functions
  // void readFileDBL(const std::string& _file, const int& start, const int& nrow,
  //                  const int& ncol, Rcpp::NumericMatrix& dataPlace);
  
  void readSimOut(const std::string& _file, const int& start, const int& switch_,
                  std::vector<SimpleMos>& dataPlace);
  
  // auxiliary functions
  double str2dbl(const char *p); // https://tinodidriksen.com/uploads/code/cpp/speed-string-to-double.cpp
  int str2int(const char *p);   // https://tinodidriksen.com/uploads/code/cpp/speed-string-to-int.cpp
  
  
private:
  // strings to hold input files
  std::string holdString;

  // input file stream
  std::ifstream holdStream;
  
  // file sizes
  size_t fileSize;
  
  // holder things for conversions
  double baseDBL;
  double expDBL;
  int    baseINT;
  // size_t outerLoop;
  // size_t innerLoop;
  
  std::string::iterator curChar, endChar;
  std::string oneNum, twoNum;
  
};


///////////////////////////////////////
//Function Definitions
///////////////////////////////////////

// convert string to double, no signs
double MPLEXReader::str2dbl(const char *p) {
  baseDBL = 0.0;

  while (*p >= '0' && *p <= '9') {
    baseDBL = (baseDBL*10.0) + (*p - '0');
    ++p;
  }
  if (*p == '.') {
    expDBL = 0.0;
    baseINT = 0;
    ++p;
    while (*p >= '0' && *p <= '9') {
      expDBL = (expDBL*10.0) + (*p - '0');
      ++p;
      ++baseINT;
    }
    baseDBL += expDBL / std::pow(10.0, baseINT); 
  }

  return baseDBL;
}

// convert string to int, no signs
int MPLEXReader::str2int(const char *p) {
  baseINT = 0;

  while (*p >= '0' && *p <= '9') {
    baseINT = (baseINT*10) + (*p - '0');
    ++p;
  }

  return baseINT;
}

// // read doubles from file into an already created matrix
// void MPLEXReader::readFileDBL(const std::string& _file, const int& start, const int& nrow,
//                             const int& ncol, Rcpp::NumericMatrix& dataPlace){
//   
//   holdString.clear();
//   
//   // open file
//   holdStream.open(_file);
//   // get amount to read
//   holdStream.seekg(0, std::ios::end);
//   fileSize = holdStream.tellg(); // size of file
//   holdStream.seekg(0); // back to beginning
// 
//   // read in everything
//   holdStream.read(&holdString[0], fileSize);
//   // set iterator at beginning
//   curChar = holdString.begin();
//   // skip first line
//   curChar += start;
//   
//   
//   // loop over matrix, fill with things
//   for(outerLoop=0; outerLoop < nrow; ++outerLoop){
//     for(innerLoop=0; innerLoop < ncol; ++innerLoop){
//       
//       // clear string holder
//       oneNum.clear();
//       
//       // get all characters part of this number
//       while(*&*curChar != ',' && *&*curChar != '\n'){
//         oneNum += *curChar;
//         ++curChar;
//       }
//       // iterate over the comma or EOL
//       ++curChar;
//       
//       // convert to double and put in matrix
//       dataPlace(outerLoop, innerLoop) = str2dbl(&oneNum[0]);
//       
//     } // end loop over columns
//   } // end loop over rows
//   
//   // close file
//   holdStream.close();
//   
// } // end double reader


// read simulation output, append to vector of objects, read by male or female
void MPLEXReader::readSimOut(const std::string& file_, const int& start_, const int& switch_,
                              std::vector<SimpleMos>& dataPlace_){
  

  // open file
  holdStream.open(file_, std::ios_base::binary);
  
  
  // get amount to read
  holdStream.seekg(0, std::ios::end);
  fileSize = holdStream.tellg(); // size of file
  holdStream.seekg(0); // back to beginning
  
  // resize buffer to current file size
  holdString.resize(fileSize);
  
  // read in everything
  holdStream.read(&holdString[0], fileSize);
  // set iterators
  curChar = holdString.begin();
  endChar = holdString.end();
  // skip first line
  curChar += start_;
  
  // switch statement for male vs female files
  switch (switch_){
    case 0: // male file
      while( curChar != endChar){
        
        // clear string holder
        oneNum.clear();
        twoNum.clear();
        
        // get all characters in time
        while(*&*curChar != ','){
          oneNum += *curChar;
          ++curChar;
        }
        // iterate over the comma
        ++curChar;
        
        // skip over age
        while(*&*curChar != ','){
          ++curChar;
        }
        // iterate over the comma
        ++curChar;
        
        // get all characters in genotype
        while(*&*curChar != '\n'){
          twoNum += *curChar;
          ++curChar;
        }
        // iterate over the EOL
        ++curChar;
        
        // add new mosquito to the vector
        dataPlace_.emplace_back(str2int(&oneNum[0]), twoNum);
      } // end loop over file
      break;
  case 1: //female file  
    while(curChar != endChar){
      
      // clear string holder
      oneNum.clear();
      twoNum.clear();
      
      // get all characters in time
      while(*&*curChar != ','){
        oneNum += *curChar;
        ++curChar;
      }
      // iterate over the comma
      ++curChar;
      
      // skip over age
      while(*&*curChar != ','){
        ++curChar;
      }
      // iterate over the comma
      ++curChar;
      
      // get all characters in genotype
      while(*&*curChar != ','){
        twoNum += *curChar;
        ++curChar;
      }
      // iterate over the EOL
      ++curChar;
      
      // skip over mate
      while(*&*curChar != '\n'){
        ++curChar;
      }
      // iterate over EOL
      ++curChar;
      
      // add new mosquito to the vector
      dataPlace_.emplace_back(str2int(&oneNum[0]), twoNum);
      
    } // end loop over file
    break;
    
  } // end switch
  
  // close file
  holdStream.close();
  
} // end int reader



#endif
