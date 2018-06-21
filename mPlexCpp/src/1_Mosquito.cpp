///////////////////////////////////////////////////////////////////////////////
//                           ____  __          ______          
//                ____ ___  / __ \/ /__  _  __/ ____/___  ____ 
//               / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
//              / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
//             /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/ 
//                                              /_/   /_/      
///////////////////////////////////////////////////////////////////////////////

#include "1_Mosquito.hpp"

// constructor & destructor
Mosquito::Mosquito(const int& _age, const std::string& _genotype) : age(_age),genotype(_genotype){};
Mosquito::~Mosquito(){};

// compiler generated move constructor & assignment operator
//Mosquito::Mosquito(Mosquito&& h) = default;
//Mosquito& Mosquito::operator=(Mosquito&& h) = default;

// print function
std::string Mosquito::print_male(){
  return std::to_string(age) + "," + genotype + "\n";
}

std::string Mosquito::print_female(){
  return std::to_string(age) + "," + genotype + "," + mate + "\n";
}
