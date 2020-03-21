///////////////////////////////////////////////////////////////////////////////
//     ________ __ __  _______ 
//    / ____/ //_//  |/  / __ \
//   / /   / ,<  / /|_/ / /_/ /
//  / /___/ /| |/ /  / / _, _/ 
//  \____/_/ |_/_/  /_/_/ |_|  
//   
///////////////////////////////////////////////////////////////////////////////

#include <string>

#ifndef bigBrother_CKMR
#define bigBrother_CKMR

///////////////////////////////////////////////////////////////////////////////
// class declaration
///////////////////////////////////////////////////////////////////////////////
/* threadsafe prng if one is created for each thread!! */


class bigBrother {
public:
  
  /* constructor & destructor */
  bigBrother(const int& init_, const int& step_);
  virtual ~bigBrother();
  //    ~prng() = default;
  
  // getters/setters
  std::string         get_ID();
  void                reset(){idMem = init;};
  
private:
  
  const int step;
  const int init;
  unsigned long long int idMem;

};

#endif