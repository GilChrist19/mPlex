//      __  _____________  ____  ______
//     /  |/  / ____/ __ \/ __ \/ ____/___  ____
//    / /|_/ / / __/ / / / / / / /   / __ \/ __ \
//   / /  / / /_/ / /_/ / /_/ / /___/ /_/ / /_/ /
//  /_/  /_/\____/_____/_____/\____/ .___/ .___/
//                                /_/   /_/

#include "1_Mosquito.hpp"
#include "2_Patch.hpp"
#include <vector>






/******************************************************************************
 * Daisy Class
******************************************************************************/
#ifndef DAISY_MPLEX
#define DAISY_MPLEX

class Daisy: public Patch{
public:
  
  // destructor
  virtual ~Daisy();
  
  // delete copy constructor and assignment operator
  Daisy(const Daisy&) = delete;
  Daisy& operator=(const Daisy&) = delete;
  
  // default move semantics
  Daisy(Daisy&&);
  Daisy& operator=(Daisy&&);
  
  // single function this class exists for
  void oneDay_layEggs();

};

#endif

/******************************************************************************
 * Multiplex oneLocus Class
******************************************************************************/
#ifndef ONELOCUS_MPLEX
#define ONELOCUS_MPLEX

class oneLocus: public Patch{
public:
  
  // destructor
  virtual ~oneLocus();
  
  // delete copy constructor and assignment operator
  oneLocus(const oneLocus&) = delete;
  oneLocus& operator=(const oneLocus&) = delete;
  
  // default move semantics
  oneLocus(oneLocus&&);
  oneLocus& operator=(oneLocus&&);
  
  // single function this class exists for
  void oneDay_layEggs();
  
};

#endif

/******************************************************************************
 * Multiplex multiLocus Class
******************************************************************************/
#ifndef MULTILOCUS_MPLEX
#define MULTILOCUS_MPLEX

class multiLocus: public Patch{
public:
  
  // destructor
  virtual ~multiLocus();
  
  // delete copy constructor and assignment operator
  multiLocus(const multiLocus&) = delete;
  multiLocus& operator=(const multiLocus&) = delete;
  
  // default move semantics
  multiLocus(multiLocus&&);
  multiLocus& operator=(multiLocus&&);
  
  // single function this class exists for
  void oneDay_layEggs();
  
};

#endif







