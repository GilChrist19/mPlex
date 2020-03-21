///////////////////////////////////////////////////////////////////////////////
//     ________ __ __  _______ 
//    / ____/ //_//  |/  / __ \
//   / /   / ,<  / /|_/ / /_/ /
//  / /___/ /| |/ /  / / _, _/ 
//  \____/_/ |_/_/  /_/_/ |_|  
//   
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
  Mosquito(const int& age_, const std::string& myID_,
           const std::string& mom_, const std::string& dad_);
  virtual ~Mosquito();

  // delete copy constructor and assignment operator
//  Mosquito(const Mosquito& h) = delete;
//  Mosquito& operator=(const Mosquito& h) = delete;

  // compiler generated move constructor & assignment operator
//  Mosquito(Mosquito&& h);
// Mosquito& operator=(Mosquito&& h);

  // getters 
  int                   get_age(){return age;};
  std::string&          get_myID(){return myID;};
  std::string&          get_mate(){return mate;};
  
  
  // setter
  void                  age_one_day(){age++;};
  void                  set_mate(const std::string& newMate){mate = newMate;};
  
  // print
  std::string           print_male();
  std::string           print_female();
  std::string           print_aquatic();
  
protected:

  int                   age;
  std::string           myID;
  std::string           mate;
  std::string           momID;
  std::string           dadID;

};

#endif