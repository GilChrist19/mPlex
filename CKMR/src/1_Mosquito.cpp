///////////////////////////////////////////////////////////////////////////////
//     ________ __ __  _______ 
//    / ____/ //_//  |/  / __ \
//   / /   / ,<  / /|_/ / /_/ /
//  / /___/ /| |/ /  / / _, _/ 
//  \____/_/ |_/_/  /_/_/ |_|  
//   
///////////////////////////////////////////////////////////////////////////////

#include "1_Mosquito.hpp"

/******************************************************************************
 * Mosquito Class
******************************************************************************/

// constructor & destructor
Mosquito::Mosquito(const int& age_, const std::string& myID_,
                   const std::string& mom_, const std::string& dad_) : 
                   age(age_), myID(myID_), momID(mom_), dadID(dad_){};
Mosquito::~Mosquito(){};

// print functions
std::string Mosquito::print_aquatic(){
  return "," + std::to_string(age) + "," + myID + "," + momID + "," + dadID + "\n";
}

std::string Mosquito::print_male(){
  return "," + std::to_string(age) + "," + myID + "," + momID + "," + dadID + "\n";
}

std::string Mosquito::print_female(){
  return "," + std::to_string(age) + "," + myID + "," + momID + "," + dadID + "," + 
    mate + "\n";
}
