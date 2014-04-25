/*! 
 @file mt19937ar.h

 @brief Mersenne Twister‚É‚æ‚é—”¶¬
 
 @author Makoto Fujisawa
 @date 2010
*/


#ifndef _MT19937AR_H_
#define _MT19937AR_H_

#include "rx_sph_commons.h"

void init_genrand(unsigned long s);

unsigned long genrand_int32(void);
long genrand_int31(void);

RXREAL genrand_real1(void);
RXREAL genrand_real2(void);
RXREAL genrand_real3(void);

RXREAL genrand_res53(void);



#endif // #ifndef _MT19937AR_H_