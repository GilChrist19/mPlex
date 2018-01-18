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

#include "Mosquito.hpp"

/* constructor & destructor */
Mosquito::Mosquito(const int& _age) : age(_age) {};
Mosquito::~Mosquito(){};

/* compiler generated move constructor & assignment operator */
Mosquito::Mosquito(Mosquito&& h) = default;
Mosquito& Mosquito::operator=(Mosquito&& h) = default;
