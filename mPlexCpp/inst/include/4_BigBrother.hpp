///////////////////////////////////////////////////////////////////////////////
//                           ____  __          ______          
//                ____ ___  / __ \/ /__  _  __/ ____/___  ____ 
//               / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
//              / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
//             /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/ 
//                                              /_/   /_/      
///////////////////////////////////////////////////////////////////////////////

#include <string>

#ifndef BIGBROTHER_MPLEX
#define BIGBROTHER_MPLEX

class BigBrother final{
public:
  /* utility functions */
  static BigBrother&          instance(); // get instance
  
  // getters/setters
  std::string                 get_ID();
  void                        reset(){idMem = 0;};
  
private:
  // constructor & destructor
  BigBrother();
  ~BigBrother();
  
  // delete all copy & move semantics
  BigBrother(const BigBrother&) = delete;
  BigBrother& operator=(const BigBrother&) = delete;
  BigBrother(BigBrother&&) = delete;
  BigBrother& operator=(BigBrother&&) = delete;
  
  // only member?
  unsigned long long int idMem;
  
};

#endif