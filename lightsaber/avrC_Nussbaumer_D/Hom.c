
#include <stddef.h>
#include <stdbool.h>
#include <string.h>
#include "Hom.h"

// ================
// Nussbaumer

static
void shift_addsub(uint16_t *des, uint16_t *src, size_t jump, size_t shift){

    uint16_t buff[64];

    for(size_t i = 0; i < shift; i++){
        buff[i] = src[i] + src[jump * 64 + ( (i + 64 - shift) % 64)];
        des[i]  = src[i] - src[jump * 64 + ( (i + 64 - shift) % 64)];
    }

    for(size_t i = shift; i < 64; i++){
        buff[i] = src[i] - src[jump * 64 + (i - shift)];
        des[i]  = src[i] + src[jump * 64 + (i - shift)];
    }

    for(size_t i = 0; i < 64; i++){
        des[jump * 64 + i] = buff[i];
    }

}

// des and src should not overlap.
static
void CT_1(uint16_t *des, uint16_t *src){

    shift_addsub(des + 0 * 64, src + 0 * 64, 2, 0);
    shift_addsub(des + 1 * 64, src + 1 * 64, 2, 0);
    shift_addsub(des + 4 * 64, src + 4 * 64, 2, 32);
    shift_addsub(des + 5 * 64, src + 5 * 64, 2, 32);

}

// des and src should not overlap.
static
void CT_2(uint16_t *des, uint16_t *src){

    shift_addsub(des + 0 * 64, src + 0 * 64, 1, 0);
    shift_addsub(des + 2 * 64, src + 2 * 64, 1, 32);
    shift_addsub(des + 4 * 64, src + 4 * 64, 1, 16);
    shift_addsub(des + 6 * 64, src + 6 * 64, 1, 48);

}

static
void addsub_shift_neg(uint16_t *des, uint16_t *src, size_t jump, size_t shift, bool neg){

    uint16_t buff[64];

    if(neg){

        for(size_t i = 0; i < 64 - shift; i++){
            buff[i] = src[jump * 64 + i] - src[i];
            des[i]  = src[jump * 64 + i] + src[i];
        }

        for(size_t i = 64 - shift; i < 64; i++){
            buff[i] = src[i] - src[jump * 64 + i];
            des[i]  = src[i] + src[jump * 64 + i];
        }

    }else{

        for(size_t i = 0; i < 64 - shift; i++){
            buff[i] = src[i] - src[jump * 64 + i];
            des[i]  = src[i] + src[jump * 64 + i];
        }

        for(size_t i = 64 - shift; i < 64; i++){
            buff[i] = src[jump * 64 + i] - src[i];
            des[i]  = src[jump * 64 + i] + src[i];
        }

    }

    for(size_t i = 0; i < 64 - shift; i++){
        des[jump * 64 + i + shift] = buff[i];
    }

    for(size_t i = 64 - shift; i < 64; i++){
        des[jump * 64 + i + shift - 64] = buff[i];
    }

}

// des and src should not overlap.
static
void GS_1(uint16_t *des, uint16_t *src){

    addsub_shift_neg(des + 0 * 64, src + 0 * 64, 2, 0, 0);
    addsub_shift_neg(des + 1 * 64, src + 1 * 64, 2, 0, 0);
    addsub_shift_neg(des + 4 * 64, src + 4 * 64, 2, 32, 1);
    addsub_shift_neg(des + 5 * 64, src + 5 * 64, 2, 32, 1);

}

// des and src should not overlap.
static
void GS_2(uint16_t *des, uint16_t *src){

    addsub_shift_neg(des + 0 * 64, src + 0 * 64, 1, 0, 0);
    addsub_shift_neg(des + 2 * 64, src + 2 * 64, 1, 32, 1);
    addsub_shift_neg(des + 4 * 64, src + 4 * 64, 1, 48, 1);
    addsub_shift_neg(des + 6 * 64, src + 6 * 64, 1, 16, 1);

}

void Nussbaumer(uint16_t *des, uint16_t *src){

    for(size_t i = 0; i < 64; i++){
        for(size_t j = 0; j < 4; j++){
            des[(j + 4) * 64 + i] = des[j * 64 + i] = src[i * 4 + j];
        }
    }

    CT_1(des, des);
    CT_2(des, des);

}

void iNussbaumer(uint16_t *des, uint16_t *src){

    uint16_t buff[4];

    GS_2(src, src);
    GS_1(src, src);

    for(size_t j = 0; j < 4; j++){
        des[j]  = src[j * 64] + src[(j + 4) * 64];
        buff[j] = src[j * 64] - src[(j + 4) * 64];
    }

    for(size_t i = 1; i < 64; i++){
        for(size_t j = 0; j < 4; j++){
            des[i * 4 + j] = src[j * 64 + i] + src[(j + 4) * 64 + i];
            des[i * 4 + j] = des[i * 4 + j] + buff[j];
            des[i * 4 + j] = des[i * 4 + j] >> 3;
            buff[j]        = src[j * 64 + i] - src[(j + 4) * 64 + i];
        }
    }

    for(size_t j = 0; j < 4; j++){
        des[j] -= buff[j];
        des[j] >>= 3;
    }

}

#define OVERFLOWING_MUL(X, Y) (X * Y)

#define KARATSUBA_N 64
static void Karatsuba(uint16_t *result_final, const uint16_t *a_1, const uint16_t *b_1) {
    uint16_t d01[KARATSUBA_N / 2 - 1];
    uint16_t d0123[KARATSUBA_N / 2 - 1];
    uint16_t d23[KARATSUBA_N / 2 - 1];
    uint16_t result_d01[KARATSUBA_N - 1];

    int32_t i, j;

    memset(result_d01, 0, (KARATSUBA_N - 1)*sizeof(uint16_t));
    memset(d01, 0, (KARATSUBA_N / 2 - 1)*sizeof(uint16_t));
    memset(d0123, 0, (KARATSUBA_N / 2 - 1)*sizeof(uint16_t));
    memset(d23, 0, (KARATSUBA_N / 2 - 1)*sizeof(uint16_t));
    memset(result_final, 0, (2 * KARATSUBA_N - 1)*sizeof(uint16_t));

    uint16_t acc1, acc2, acc3, acc4, acc5, acc6, acc7, acc8, acc9, acc10;


    for (i = 0; i < KARATSUBA_N / 4; i++) {
        acc1 = a_1[i]; //a0
        acc2 = a_1[i + KARATSUBA_N / 4]; //a1
        acc3 = a_1[i + 2 * KARATSUBA_N / 4]; //a2
        acc4 = a_1[i + 3 * KARATSUBA_N / 4]; //a3
        for (j = 0; j < KARATSUBA_N / 4; j++) {

            acc5 = b_1[j]; //b0
            acc6 = b_1[j + KARATSUBA_N / 4]; //b1

            result_final[i + j + 0 * KARATSUBA_N / 4] =
                result_final[i + j + 0 * KARATSUBA_N / 4] +
                OVERFLOWING_MUL(acc1, acc5);
            result_final[i + j + 2 * KARATSUBA_N / 4] =
                result_final[i + j + 2 * KARATSUBA_N / 4] +
                OVERFLOWING_MUL(acc2, acc6);

            acc7 = acc5 + acc6; //b01
            acc8 = acc1 + acc2; //a01
            d01[i + j] = d01[i + j] + (uint16_t)(acc7 * (uint64_t)acc8);
            //--------------------------------------------------------

            acc7 = b_1[j + 2 * KARATSUBA_N / 4]; //b2
            acc8 = b_1[j + 3 * KARATSUBA_N / 4]; //b3
            result_final[i + j + 4 * KARATSUBA_N / 4] =
                result_final[i + j + 4 * KARATSUBA_N / 4] +
                OVERFLOWING_MUL(acc7, acc3);

            result_final[i + j + 6 * KARATSUBA_N / 4] =
                result_final[i + j + 6 * KARATSUBA_N / 4] +
                OVERFLOWING_MUL(acc8, acc4);

            acc9 = acc3 + acc4;
            acc10 = acc7 + acc8;
            d23[i + j] = d23[i + j] + OVERFLOWING_MUL(acc9, acc10);
            //--------------------------------------------------------

            acc5 = acc5 + acc7; //b02
            acc7 = acc1 + acc3; //a02
            result_d01[i + j + 0 * KARATSUBA_N / 4] =
                result_d01[i + j + 0 * KARATSUBA_N / 4] +
                OVERFLOWING_MUL(acc5, acc7);

            acc6 = acc6 + acc8; //b13
            acc8 = acc2 + acc4;
            result_d01[i + j + 2 * KARATSUBA_N / 4] =
                result_d01[i + j + 2 * KARATSUBA_N / 4] +
                OVERFLOWING_MUL(acc6, acc8);

            acc5 = acc5 + acc6;
            acc7 = acc7 + acc8;
            d0123[i + j] = d0123[i + j] + OVERFLOWING_MUL(acc5, acc7);
        }
    }

    // 2nd last stage

    for (i = 0; i < KARATSUBA_N / 2 - 1; i++) {
        d0123[i] = d0123[i] - result_d01[i + 0 * KARATSUBA_N / 4] - result_d01[i + 2 * KARATSUBA_N / 4];
        d01[i] = d01[i] - result_final[i + 0 * KARATSUBA_N / 4] - result_final[i + 2 * KARATSUBA_N / 4];
        d23[i] = d23[i] - result_final[i + 4 * KARATSUBA_N / 4] - result_final[i + 6 * KARATSUBA_N / 4];
    }

    for (i = 0; i < KARATSUBA_N / 2 - 1; i++) {
        result_d01[i + 1 * KARATSUBA_N / 4] = result_d01[i + 1 * KARATSUBA_N / 4] + d0123[i];
        result_final[i + 1 * KARATSUBA_N / 4] = result_final[i + 1 * KARATSUBA_N / 4] + d01[i];
        result_final[i + 5 * KARATSUBA_N / 4] = result_final[i + 5 * KARATSUBA_N / 4] + d23[i];
    }

    // Last stage
    for (i = 0; i < KARATSUBA_N - 1; i++) {
        result_d01[i] = result_d01[i] - result_final[i] - result_final[i + KARATSUBA_N];
    }

    for (i = 0; i < KARATSUBA_N - 1; i++) {
        result_final[i + 1 * KARATSUBA_N / 2] = result_final[i + 1 * KARATSUBA_N / 2] + result_d01[i];
    }

}

static
void negacyclic_Karatsuba(uint16_t *des, const uint16_t *src1, const uint16_t *src2){

    uint16_t buff[2 * KARATSUBA_N];

    Karatsuba(buff, src1, src2);

    for(size_t i = 0; i < KARATSUBA_N - 1; i++){
        des[i] = buff[i] - buff[i + KARATSUBA_N];
    }
    des[KARATSUBA_N - 1] = buff[KARATSUBA_N - 1];

}

static
void negacyclic_Karatsuba_striding(uint16_t *des, const uint16_t *src1, const uint16_t *src2){

    uint16_t res_p00[16], res_p01[16], res_p02[16],
             res_p10[16], res_p11[16], res_p12[16],
             res_p20[16], res_p21[16], res_p22[16];

    uint16_t p00, p01, p02,
             p10, p11, p12,
             p20, p21, p22;

    uint16_t _p01, _p11, _p21;

    uint16_t buff[3];

    memset(res_p00, 0, sizeof(res_p00));
    memset(res_p01, 0, sizeof(res_p01));
    memset(res_p02, 0, sizeof(res_p02));
    memset(res_p10, 0, sizeof(res_p10));
    memset(res_p11, 0, sizeof(res_p11));
    memset(res_p12, 0, sizeof(res_p12));
    memset(res_p20, 0, sizeof(res_p20));
    memset(res_p21, 0, sizeof(res_p21));
    memset(res_p22, 0, sizeof(res_p22));

    for(size_t i = 0; i < 16; i++){

        // Load 4, cache 9.
        p00 = src1[4 * i + 0];
        p02 = src1[4 * i + 1];
        p20 = src1[4 * i + 2];
        p22 = src1[4 * i + 3];

        p01 = p00 + p02;
        p21 = p20 + p22;

        p10 = p00 + p20;
        p12 = p02 + p22;

        p11 = p10 + p12;

        for(size_t j = 0; j < 16 - i; j++){

            res_p00[i + j]      += p00 * src2[4 * j + 0];
            res_p02[i + j]      += p02 * src2[4 * j + 1];
            _p01                 = src2[4 * j + 0] + src2[4 * j + 1];
            res_p01[i + j]      += p01 * _p01;

            res_p20[i + j]      += p20 * src2[4 * j + 2];
            res_p22[i + j]      += p22 * src2[4 * j + 3];
            _p21                 = src2[4 * j + 2] + src2[4 * j + 3];
            res_p21[i + j]      += p21 * _p21;

            res_p10[i + j]      += p10 * (src2[4 * j + 0] + src2[4 * j + 2]);
            res_p12[i + j]      += p12 * (src2[4 * j + 1] + src2[4 * j + 3]);
            _p11                 = _p01 + _p21;
            res_p11[i + j]      += p11 * _p11;

        }
        for(size_t j = 16 - i; j < 16; j++){

            res_p00[i + j - 16] -= p00 * src2[4 * j + 0];
            res_p02[i + j - 16] -= p02 * src2[4 * j + 1];
            _p01                 = src2[4 * j + 0] + src2[4 * j + 1];
            res_p01[i + j - 16] -= p01 * _p01;

            res_p20[i + j - 16] -= p20 * src2[4 * j + 2];
            res_p22[i + j - 16] -= p22 * src2[4 * j + 3];
            _p21                 = src2[4 * j + 2] + src2[4 * j + 3];
            res_p21[i + j - 16] -= p21 * _p21;

            res_p10[i + j - 16] -= p10 * (src2[4 * j + 0] + src2[4 * j + 2]);
            res_p12[i + j - 16] -= p12 * (src2[4 * j + 1] + src2[4 * j + 3]);
            _p11                 = _p01 + _p21;
            res_p11[i + j - 16] -= p11 * _p11;

        }
    }

    // Load 9, store 6.
    for(size_t i = 0; i < 16; i++){
        res_p01[i] = res_p01[i] - res_p00[i] - res_p02[i];
        res_p11[i] = res_p11[i] - res_p10[i] - res_p12[i];
        res_p21[i] = res_p21[i] - res_p20[i] - res_p22[i];

        res_p10[i] = res_p10[i] - res_p00[i] - res_p20[i];
        res_p11[i] = res_p11[i] - res_p01[i] - res_p21[i];
        res_p12[i] = res_p12[i] - res_p02[i] - res_p22[i];

        res_p10[i] += res_p02[i];
        res_p20[i] += res_p12[i];
    }

    // Load 11, store 4, Cache 3.
    des[0]  = res_p00[0] - res_p20[15];
    des[1]  = res_p01[0] - res_p21[15];
    des[2]  = res_p10[0] - res_p22[15];
    des[3]  = res_p11[0];
    buff[0] = res_p20[0];
    buff[1] = res_p21[0];
    buff[2] = res_p22[0];
    // Load 7, store 4, cache 3.
    for(size_t i = 1; i < 15; i++){
        des[4 * i + 0] = buff[0] + res_p00[i];
        des[4 * i + 1] = buff[1] + res_p01[i];
        des[4 * i + 2] = buff[2] + res_p10[i];
        des[4 * i + 3] =           res_p11[i];
        buff[0]        =           res_p20[i];
        buff[1]        =           res_p21[i];
        buff[2]        =           res_p22[i];
    }
    // Load 4, store 4.
    des[60] = buff[0] + res_p00[15];
    des[61] = buff[1] + res_p01[15];
    des[62] = buff[2] + res_p10[15];
    des[63] =           res_p11[15];

}

void BiHom(uint16_t *des, uint16_t *src1, uint16_t *src2){
    for(size_t i = 0; i < 8; i++){
        //negacyclic_Karatsuba(des + i * 64, src1 + i * 64, src2 + i * 64);
        negacyclic_Karatsuba_striding(des + i * 64, src1 + i * 64, src2 + i * 64);
    }
}



