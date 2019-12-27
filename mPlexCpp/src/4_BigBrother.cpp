///////////////////////////////////////////////////////////////////////////////
//                           ____  __          ______          
//                ____ ___  / __ \/ /__  _  __/ ____/___  ____ 
//               / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
//              / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
//             /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/ 
//                                              /_/   /_/      
///////////////////////////////////////////////////////////////////////////////

#include "4_BigBrother.hpp"

// constructor & destructor
BigBrother::BigBrother(){
  idMem=0;
};
BigBrother::~BigBrother(){};

// utility method
BigBrother& BigBrother::instance(){
  static BigBrother instance;
  return instance;
}

// getters/setters
std::string BigBrother::get_ID(){
  
  idMem++; // use 0 as a default for setting stuff
  return(std::to_string(idMem));

};








