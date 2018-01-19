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

#include <string>

class Mosquito {
public:

  /* constructor & destructor */
  Mosquito(const int& _age);
  ~Mosquito();

  /* delete copy constructor and assignment operator */
  Mosquito(const Mosquito& h) = delete;
  Mosquito& operator=(const Mosquito& h) = delete;

  /* compiler generated move constructor & assignment operator */
  Mosquito(Mosquito&& h);
  Mosquito& operator=(Mosquito&& h);

  /* accessors */
  int                   get_age(){return age;};
  std::string&          get_mate(){return mate;};
  std::string&          get_genotype(){return genotype;};
  void                  set_mate(){}

private:

  int                   age;
  std::string           mate;
  std::string           genotype;

};
