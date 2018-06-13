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

#include "3_DriveDefinitions.hpp"



/******************************************************************************
 * Constructor & Destructor
 ******************************************************************************/


oneLocus::~oneLocus(){};

/******************************************************************************
 * Default (compiler-generated) move semantics
 ******************************************************************************/
oneLocus::oneLocus(oneLocus&& d) = default;
oneLocus& oneLocus::operator=(oneLocus&& d) = default;

/******************************************************************************
 * Define mating function here
 ******************************************************************************/

/******************************************************************************
 * Lay eggs
 ******************************************************************************/
void oneLocus::oneDay_layEggs(){
  
  
  Rcpp::Rcout<<"You chose the oneLocus function!"<<std::endl;
  
  
  
}