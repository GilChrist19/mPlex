///////////////////////////////////////////////////////////////////////////////
//                           ____  __          ______          
//                ____ ___  / __ \/ /__  _  __/ ____/___  ____ 
//               / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
//              / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
//             /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/ 
//                                              /_/   /_/      
///////////////////////////////////////////////////////////////////////////////

#include <string>

/******************************************************************************
 * Mosquito Class
******************************************************************************/

#ifndef MOSQUITO_MPLEX
#define MOSQUITO_MPLEX

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
  void                  set_parents(std::string mom_, std::string dad_)
                                    {momID = mom_; dadID = dad_;};
  
  // print
  std::string           print_male();
  std::string           print_female();
  std::string           print_maleFam();
  std::string           print_femaleFam();
  
protected:

  int                   age;
  std::string           genotype;
  std::string           mate;
  std::string           momID;
  std::string           dadID;

};

#endif