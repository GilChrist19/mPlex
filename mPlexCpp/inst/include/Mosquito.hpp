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


private:

  int                   age;


};
