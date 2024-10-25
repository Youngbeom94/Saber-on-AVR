#ifndef RANDOMBYTES_H
#define RANDOMBYTES_H
#include <stdint.h>

void init_randombytes();
int randombytes(uint8_t *buf, uint64_t n);

#endif
