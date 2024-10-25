#include <stdint.h>
#include "randombytes.h"

//TODO Maybe we do not want to use the hardware RNG for all randomness, but instead only read a seed and then expand that using fips202.

int randombytes(uint8_t *buf, uint64_t n)
{
  for(uint16_t cnt_i  = 0;  cnt_i < n ; cnt_i ++) buf[cnt_i] = cnt_i;
  
  return n;
}

