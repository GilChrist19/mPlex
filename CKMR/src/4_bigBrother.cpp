///////////////////////////////////////////////////////////////////////////////
//     ________ __ __  _______ 
//    / ____/ //_//  |/  / __ \
//   / /   / ,<  / /|_/ / /_/ /
//  / /___/ /| |/ /  / / _, _/ 
//  \____/_/ |_/_/  /_/_/ |_|  
//   
///////////////////////////////////////////////////////////////////////////////

#include "4_bigBrother.hpp"

///////////////////////////////////////////////////////////////////////////////
// constructor & destructor
///////////////////////////////////////////////////////////////////////////////
// constructor & destructor
bigBrother::bigBrother(const int& init_, const int& step_) : 
  init(init_), idMem(init_), step(step_){
};
bigBrother::~bigBrother(){};

// getters/setters
std::string bigBrother::get_ID(){
  
  idMem+=step; // use 0:step as a default for setting stuff
  
  return(std::to_string(idMem));
};









// auxiliary function?

// // try to write something like dat.atable
// static inline void reverse(char *upp, char *low)
// {
//   upp--;
//   while (upp>low) {
//     char tmp = *upp;
//     *upp = *low;
//     *low = tmp;
//     upp--;
//     low++;
//   }
// }
// 
// 
// void writeInt64(unsigned long long int x, char **pch)
// {
//   
//   char *ch = *pch;
//   
//   
//   char *low = ch;
//   
//   
//   do { *ch++ = '0'+x%10; x/=10;} while (x>0);
//   
//   reverse(ch, low);
//     
//   *pch = ch;
//   
// }
// 
// // random from google
// //https://stackoverflow.com/questions/8257714/how-to-convert-an-int-to-string-in-c
// 
// char * itoa (int value, char *result)
// {
// 
//     char* ptr = result, *ptr1 = result, tmp_char;
// 
//     do {*ptr++ = '0' + value%10; value/=10;} while ( value>0 );
// 
//     // add null terminator
//     *ptr-- = '\0';
//     
//     // reverse string
//     while (ptr1 < ptr) {
//         tmp_char = *ptr;
//         *ptr--= *ptr1;
//         *ptr1++ = tmp_char;
//     }
//     
//     // return pointer
//     return result;
// }
// 
// 
// void jared (int value, char *result)
// {
//   
//   char* ptr = result, *ptr1 = result, tmp_char;
//   
//   do {*ptr++ = '0' + value%10; value/=10;} while ( value>0 );
//   
//   // add null terminator
//   *ptr-- = '\0';
//   
//   // reverse string
//   while (ptr1 < ptr) {
//     tmp_char = *ptr;
//     *ptr--= *ptr1;
//     *ptr1++ = tmp_char;
//   }
//   
// }













