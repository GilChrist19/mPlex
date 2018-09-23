///////////////////////////////////////////////////////////////////////////////
//                           ____  __          ______          
//                ____ ___  / __ \/ /__  _  __/ ____/___  ____ 
//               / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
//              / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
//             /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/ 
//                                              /_/   /_/      
///////////////////////////////////////////////////////////////////////////////

#include "1_Mosquito.hpp"

/******************************************************************************
 * Mosquito Class
******************************************************************************/

// constructor & destructor
Mosquito::Mosquito(const int& _age, const std::string& _genotype) : age(_age),genotype(_genotype){};
Mosquito::~Mosquito(){};

// print function
std::string Mosquito::print_male(){
  return std::to_string(age) + "," + genotype + "\n";
}

std::string Mosquito::print_female(){
  return std::to_string(age) + "," + genotype + "," + mate + "\n";
}


/******************************************************************************
 * Mosquito Familial Class
******************************************************************************/

// constructor and destructor
MosFam::MosFam(const int& _age,
               const std::string& _genotype,
               const std::string& _momID,
               const std::string& _dadID) : Mosquito::Mosquito(_age, _genotype), momID(_momID), dadID(_dadID){};
MosFam::~MosFam(){};

// print functions
std::string MosFam::print_male(){
  return std::to_string(age) + "," + genotype + "," + momID + "," + dadID + "\n";
}

std::string MosFam::print_female(){
  return std::to_string(age) + "," + genotype + "," + momID + "," + dadID + "," + 
    mate + "\n";
}


