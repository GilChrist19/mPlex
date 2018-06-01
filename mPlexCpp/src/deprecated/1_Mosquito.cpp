/*                  ____  __
 *       ____ ___  / __ \/ /__  _  __
 *      / __ `__ \/ /_/ / / _ \| |/_/
 *     / / / / / / ____/ /  __/>  <
 *    /_/ /_/ /_/_/   /_/\___/_/|_|
 *
 *    Marshall Lab
 *    Mosquito class
 *    January 2018
*/

#include "1_Mosquito.hpp"


/******************************************************************************
 * Base Class 
******************************************************************************/

// constructor & destructor
Mosquito::Mosquito(const int& _age, const std::string& _genotype) : age(_age),genotype(_genotype){};
Mosquito::~Mosquito(){};

// compiler generated move constructor & assignment operator
Mosquito::Mosquito(Mosquito&& h) = default;
Mosquito& Mosquito::operator=(Mosquito&& h) = default;

// print function
std::string Mosquito::print(){
  return std::to_string(age) + "," + genotype + "\n";
}

/******************************************************************************
 * Female Class 
******************************************************************************/
// destructor
FemaleMos::~FemaleMos(){};

// compiler generated move constructor & assignment operator
FemaleMos::FemaleMos(FemaleMos&& h) = default;
FemaleMos& FemaleMos::operator=(FemaleMos&& h) = default;

//print function
std::string FemaleMos::print(){
  return std::to_string(age) + "," + genotype + "," + mate + "\n";
}
