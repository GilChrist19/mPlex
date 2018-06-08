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

#ifndef MOSQUITO_MPLEX
#define MOSQUITO_MPLEX

#include <string>

class Mosquito {
public:

  // constructor & destructor
  Mosquito(const int& _age, const std::string& _genotype);
  virtual ~Mosquito();

  // delete copy constructor and assignment operator
//  Mosquito(const Mosquito& h) = delete;
//  Mosquito& operator=(const Mosquito& h) = delete;

  // compiler generated move constructor & assignment operator
//  Mosquito(Mosquito&& h);
// Mosquito& operator=(Mosquito&& h);

  // getters 
  int                   get_age(){return age;};
  std::string&          get_genotype(){return genotype;};
  std::string&          get_mate(){return mate;};
  
  
  // setter
  void                  age_one_day(){age++;};
  void                  set_mate(std::string newMate){mate = newMate;};
  
  
  // print
  std::string           print_male();
  std::string           print_female();
  
protected:

  int                   age;
  std::string           genotype;
  std::string           mate;

};

#endif