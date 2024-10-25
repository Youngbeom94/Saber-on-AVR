#ifndef HOM_H
#define HOM_H

#include <stdint.h>

void negacyclic_4x4(int32_t *des, int32_t *src1, int32_t *src2);
void negacyclic_acc_4x4(int32_t *des, int32_t *src1, int32_t *src2);

void TC(int32_t *des, int32_t *src);
void iTC(int32_t *des, int32_t *src);

void TC_Hom(int32_t *des, int32_t *src);
void iTC_Hom(int32_t *des, int32_t *src);
void TC_striding_mul(int32_t *des, int32_t *src1, int32_t *src2);

void TMVP(int32_t *des, int32_t *srcM, int32_t *srcV);
// TMVP_TC_Hom_V transforming a vector of 4 elements to a vector of 7 elements.
void TMVP_TC_Hom_V(int32_t *des, int32_t *src, size_t jump);
// TMVP_TC_Hom_I transforming a vector of 7 elements to a vector of 7 elements.
void TMVP_TC_Hom_I(int32_t *des, int32_t *src, size_t jump);
void TMVP_TC4_negacyclic16_mul(int32_t *des, int32_t *src1, int32_t *src2);
void TMVP_Hom_M(int32_t *des, int32_t *src);

void Nussbaumer(int32_t *des, int32_t *src);
void Nussbaumer_uint16(int32_t *des, uint16_t *src);
void iNussbaumer(int32_t *des, int32_t *src);
void iNussbaumer_uint16(uint16_t *des, int32_t *src);
void iNussbaumer_acc_uint16(uint16_t *des, int32_t *src);
void BiHom(int32_t *des, int32_t *src_M, int32_t *src_V);

#endif

