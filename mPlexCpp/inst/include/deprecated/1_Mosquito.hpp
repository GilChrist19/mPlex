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

  // constructor & destructor
  Mosquito(const int& _age, const std::string& _genotype);
  virtual ~Mosquito();

  // delete copy constructor and assignment operator
  Mosquito(const Mosquito& h) = delete;
  Mosquito& operator=(const Mosquito& h) = delete;

  // compiler generated move constructor & assignment operator
  Mosquito(Mosquito&& h);
  Mosquito& operator=(Mosquito&& h);

  // getters 
  int                   get_age(){return age;};
  std::string&          get_genotype(){return genotype;};
  
  // setter
  void                  age_one_day(){age++;};
  
  // print
  std::string           print();
  
protected:

  int                   age;
  std::string           genotype;

};


class FemaleMos: public Mosquito{
public:
  
  // destructor
  virtual ~FemaleMos();
  
  // delete copy constructor and assignment operator
  FemaleMos(const FemaleMos& h) = delete;
  FemaleMos& operator=(const FemaleMos& h) = delete;
  
  // compiler generated move constructor & assignment operator
  FemaleMos(FemaleMos&& h);
  FemaleMos& operator=(FemaleMos&& h);
  
  // setters
  void                  set_mate(std::string newMate){mate = newMate;};
  
  // getters
  std::string&          get_mate(){return mate;};
  
  // print
  std::string           print();
  
private:
  
  std::string           mate;

};

