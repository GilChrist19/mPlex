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
Mosquito::Mosquito(const int& _age, const std::string& _genotype, const std::string& _mateGenotype) : 
  age(_age),genotype(_genotype),mate(_mateGenotype){};
Mosquito::~Mosquito(){};

// print function
std::string Mosquito::print_male(){
  return std::to_string(age) + "," + genotype + "\n";
}

std::string Mosquito::print_female(){
  return std::to_string(age) + "," + genotype + "," + mate + "\n";
}

// print functions
std::string Mosquito::print_maleFam(){
  return std::to_string(age) + "," + genotype + "," + momID + "," + dadID + "\n";
}

std::string Mosquito::print_femaleFam(){
  return std::to_string(age) + "," + genotype + "," + momID + "," + dadID + "," + 
    mate + "\n";
}
