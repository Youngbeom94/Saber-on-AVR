
#include <stddef.h>
#include <stdbool.h>
#include <string.h>
#include "Hom.h"

// ================
// Toeplitz TC-4

__flash
static
int32_t TC4_trunc[7][7] = {
{ 1,  0,  0,  0, 0, 0, 0},
{ 1,  1,  1,  1, 0, 0, 0},
{ 1, -1,  1, -1, 0, 0, 0},
{ 1,  2,  4,  8, 0, 0, 0},
{ 1, -2,  4, -8, 0, 0, 0},
{ 8,  4,  2,  1, 0, 0, 0},
{ 0,  0,  0,  1, 0, 0, 0}
};

__flash
static
int32_t iTC4_modified[7][7] = {
{   1,    0,    0,    0,     0,    0,    0},
{-360, -120,  -40,    5,     3,    8, -360},
{ -30,   16,   16,   -1,    -1,    0,   96},
{  45,   27,   -7,   -1,     0,   -1,   45},
{   6,   -4,   -4,    1,     1,    0, -120},
{ -90,  -60,   20,    5,    -3,    2,  -90},
{   0,    0,    0,    0,     0,    0,    1}
};

__flash
static
int32_t iTC4_modified_scale[7] = {
1, -1527099483, -1431655765, 954437177, -1431655765, -1527099483, 1
};

__flash
static
int32_t iTC4_T_modified_scale[7] = {
1, -1431655765, 954437177, 954437177, -286331153, -1527099483, 1
};

// int32_t iTC4_T_modified[7][7] = {
// {4,   -8,   -5,   10,    1,    -2, 0},
// {0,   -4,    4,    9,   -1,    -2, 0},
// {0,   -4,   12,   -7,   -3,     2, 0},
// {0,    2,   -3,   -4,    3,     2, 0},
// {0,    2,   -5,    0,    5,    -2, 0},
// {0,    4,    0,   -5,    0,     1, 0},
// {0,   -4,    8,    5,  -10,    -1, 2}
// };

// 1/8, 1/4, 1/2, 1/2
__flash
static
int32_t TC4_trunc_T_modified[7][7] = {
{ 2,  4,  4,  1,  1, 32,  0},
{ 0,  2, -2,  1, -1,  8,  0},
{ 0,  1,  1,  1,  1,  2,  0},
{ 0,  1, -1,  2, -2,  1,  1},
{ 0,  0,  0,  0,  0,  0,  0},
{ 0,  0,  0,  0,  0,  0,  0},
{ 0,  0,  0,  0,  0,  0,  0}
};

void TC(int32_t *des, int32_t *src){

    int32_t buff[7];

    for(size_t i = 4; i < 7; i++){
        buff[i] = 0;
    }

    for(size_t i = 0; i < 4; i++){
        for(size_t j = 0; j < 7; j++){
            buff[j] = 0;
            for(size_t k = 0; k < 4; k++){
                buff[j] += TC4_trunc[j][k] * src[i * 4 + k];
            }
        }
        for(size_t j = 0; j < 7; j++){
            des[j * 4 + i] = buff[j];
        }
    }

}

void TC_Hom(int32_t *des, int32_t *src){
    for(size_t i = 0; i < 32; i++){
        TC(des + i * 28, src + i * 16);
    }
}

void iTC(int32_t *des, int32_t *src){

    int32_t TC4_buff[7];
    int32_t buff[7];

    // Apply Toom-4 inversion matrix.
    for(size_t j = 0; j < 7; j++){
        TC4_buff[j] = 0;
        for(size_t k = 0; k < 7; k++){
            TC4_buff[j] += iTC4_modified[j][k] * src[k * 4 + 0];
        }
        TC4_buff[j] = TC4_buff[j] * iTC4_modified_scale[j];
    }
    // Multiply by powers of two.
    TC4_buff[1] >>= 2;
    TC4_buff[2] >>= 3;
    TC4_buff[3] >>= 1;
    TC4_buff[4] >>= 3;
    TC4_buff[5] >>= 2;
    for(size_t j = 0; j < 7; j++){
        src[j * 4 + 0] = TC4_buff[j];
    }

    for(size_t j = 4; j < 7; j++){
        buff[j] = TC4_buff[j];
    }

    for(size_t i = 1; i < 4; i++){
        for(size_t j = 0; j < 7; j++){
            TC4_buff[j] = 0;
            for(size_t k = 0; k < 7; k++){
                TC4_buff[j] += iTC4_modified[j][k] * src[k * 4 + i];
            }
            TC4_buff[j] = TC4_buff[j] * iTC4_modified_scale[j];
        }
        // Multiply by powers of two.
        TC4_buff[1] >>= 2;
        TC4_buff[2] >>= 3;
        TC4_buff[3] >>= 1;
        TC4_buff[4] >>= 3;
        TC4_buff[5] >>= 2;
        for(size_t j = 0; j < 4; j++){
            src[j * 4 + i] = TC4_buff[j];
        }
        for(size_t j = 0; j < 3; j++){
            src[j * 4 + i] += buff[j + 4];
        }
        for(size_t j = 4; j < 7; j++){
            buff[j] = TC4_buff[j];
        }

    }

    for(size_t j = 0; j < 3; j++){
        src[j * 4 + 0] -= buff[j + 4];
    }

    // Export the result.
    for(size_t i = 0; i < 4; i++){
        for(size_t j = 0; j < 4; j++){
            des[i * 4 + j] = src[j * 4 + i];
        }
    }

}

void iTC_Hom(int32_t *des, int32_t *src){
    for(size_t i = 0; i < 32; i++){
        iTC(des + i * 16, src + i * 28);
    }
}

void negacyclic_4x4(int32_t *des, int32_t *src1, int32_t *src2){
    int32_t buff[4];
    buff[0] = src1[0] * src2[0] - src1[3] * src2[1]
            - src1[2] * src2[2] - src1[1] * src2[3];
    buff[1] = src1[1] * src2[0] + src1[0] * src2[1]
            - src1[3] * src2[2] - src1[2] * src2[3];
    buff[2] = src1[2] * src2[0] + src1[1] * src2[1]
            + src1[0] * src2[2] - src1[3] * src2[3];
    buff[3] = src1[3] * src2[0] + src1[2] * src2[1]
            + src1[1] * src2[2] + src1[0] * src2[3];
    for(size_t i = 0; i < 4; i++){
        des[i] = buff[i];
    }
}

void negacyclic_acc_4x4(int32_t *des, int32_t *src1, int32_t *src2){
    int32_t buff[4];
    buff[0] = src1[0] * src2[0] - src1[3] * src2[1]
            - src1[2] * src2[2] - src1[1] * src2[3];
    buff[1] = src1[1] * src2[0] + src1[0] * src2[1]
            - src1[3] * src2[2] - src1[2] * src2[3];
    buff[2] = src1[2] * src2[0] + src1[1] * src2[1]
            + src1[0] * src2[2] - src1[3] * src2[3];
    buff[3] = src1[3] * src2[0] + src1[2] * src2[1]
            + src1[1] * src2[2] + src1[0] * src2[3];
    for(size_t i = 0; i < 4; i++){
        des[i] += buff[i];
    }
}

// This function computes the product of two size-4 polynomials in Z_Q[x]
// using Toom-4 with the point set {0, 1, -1, 2, -2, 1/2, \infty}.
// Notice that matrices are modifed to ensure the well-defineness over Z_{32}.
void TC_striding_mul(int32_t *des, int32_t *src1, int32_t *src2){

    int32_t src1_extended[7][4], src2_extended[7][4];
    int32_t res[7][4];

    // Apply Toom-4 evaluation matrix.
    TC((int32_t*)&src1_extended[0][0], src1);
    TC((int32_t*)&src2_extended[0][0], src2);

    // Compute small-dimensional products.
    for(size_t i = 0; i < 7; i++){
        negacyclic_4x4((int32_t*)&res[i][0], (int32_t*)&src1_extended[i][0], (int32_t*)&src2_extended[i][0]);
    }

    iTC(des, (int32_t*)&res[0][0]);

}

void TMVP(int32_t *des, int32_t *srcM, int32_t *srcV)
{

    int32_t buff[4];

    for(size_t i = 0; i < 4; i++){
        buff[i] = 0;
        for(size_t j = 0; j < 4; j++){
            buff[i] += srcM[4 - 1 - i + j] * srcV[j];
        }
    }

    memmove(des, buff, sizeof(buff));

}

// TMVP_TC_Hom_V transforming a vector of 4 elements to a vector of 7 elements.
void TMVP_TC_Hom_V(int32_t *des, int32_t *src, size_t jump){

    int32_t buff[7];

    for(size_t i = 0; i < 7; i++){
        buff[i] = 0;
        for(size_t j = 0; j < 4; j++){
            buff[i] += TC4_trunc[i][j] * src[j * jump];
        }
    }

    for(size_t i = 0; i < 7; i++){
        des[i * jump] = buff[i];
    }

}

// TMVP_TC_Hom_I transforming a vector of 7 elements to a vector of 7 elements.
void TMVP_TC_Hom_I(int32_t *des, int32_t *src, size_t jump){

    int32_t buff[4];

    for(size_t i = 0; i < 4; i++){
        buff[i] = 0;
        for(size_t j = 0; j < 7; j++){
            buff[i] += TC4_trunc_T_modified[i][j] * src[j * jump];
        }
    }

    des[0 * jump] = buff[3] >> 1;
    des[1 * jump] = buff[2] >> 1;
    des[2 * jump] = buff[1] >> 2;
    des[3 * jump] = buff[0] >> 3;

}

// TMVP_TC_Hom_M transforming a size-16 polynomial into its Toeplitz form followed by
// the transformation mapping seven 4x4 Toeplitz matrices to seven 4x4 Toeplitz matrices.
// We apply on-the-fly extraction for the small Toeplitz matrices and merge the negation with
// the map iTC4_T_modified itself.
static
void TMVP_TC_Hom_M(int32_t *des, int32_t *src){

    int32_t src_buff[4];
    int32_t buff[7];

    for(size_t i = 0; i < 4; i++){

        src_buff[0] = src[15 - i];
        src_buff[1] = src[11 - i];
        src_buff[2] = src[7 - i];
        src_buff[3] = src[3 - i];

        buff[0] =    3 * src_buff[0] - 6 * src_buff[1] -  5 * src_buff[2] + 10 * src_buff[3];
        buff[1] =        src_buff[0] - 2 * src_buff[1] +  4 * src_buff[2] +  9 * src_buff[3];
        buff[2] =    3 * src_buff[0] - 6 * src_buff[1] + 12 * src_buff[2] -  7 * src_buff[3];
        buff[3] = (-3) * src_buff[0]                   -  3 * src_buff[2] -  4 * src_buff[3];
        buff[4] = (-5) * src_buff[0] + 4 * src_buff[1] -  5 * src_buff[2];
        buff[5] =                      3 * src_buff[1]                    -  5 * src_buff[3];
        buff[6] =   10 * src_buff[0] - 3 * src_buff[1] +  6 * src_buff[2] +  5 * src_buff[3];

        for(size_t j = 0; j < 7; j++){
            des[j * 8 + i] = buff[j] * iTC4_T_modified_scale[j];
        }

    }

    for(size_t i = 0; i < 3; i++){

        src_buff[0] = src[11 - i];
        src_buff[1] = src[7 - i];
        src_buff[2] = src[3 - i];
        src_buff[3] = src[15 - i];

        buff[0] =    3 * src_buff[0] - 6 * src_buff[1] -  5 * src_buff[2] - 10 * src_buff[3];
        buff[1] =        src_buff[0] - 2 * src_buff[1] +  4 * src_buff[2] -  9 * src_buff[3];
        buff[2] =    3 * src_buff[0] - 6 * src_buff[1] + 12 * src_buff[2] +  7 * src_buff[3];
        buff[3] = (-3) * src_buff[0]                   -  3 * src_buff[2] +  4 * src_buff[3];
        buff[4] = (-5) * src_buff[0] + 4 * src_buff[1] -  5 * src_buff[2];
        buff[5] =                      3 * src_buff[1]                    +  5 * src_buff[3];
        buff[6] =   10 * src_buff[0] - 3 * src_buff[1] +  6 * src_buff[2] -  5 * src_buff[3];

        for(size_t j = 0; j < 7; j++){
            des[j * 8 + i + 4] = buff[j] * iTC4_T_modified_scale[j];
        }

    }

}

// Compute the product of two size-16 polynomials in Z_Q[x] / (x^16 + 1) via Toeplitz-TC4.
// Notice that in our use case, the input of this function has coefficients
// in Z_{2^24} so we don't really need to load the most-significant byte of the inputs.
// (We still need to compute the 32-bit results.)
void TMVP_TC4_negacyclic16_mul(int32_t *des, int32_t *src1, int32_t *src2){

    int32_t src1_V_full[7][4];
    int32_t src2_Toeplitz_full[7][8];
    int32_t res_V_full[7][4];

    for(size_t i = 0; i < 4; i++){
        TMVP_TC_Hom_V((int32_t*)&src1_V_full[0][i], src1 + i, 4);
    }

    TMVP_TC_Hom_M((int32_t*)&src2_Toeplitz_full[0][0], src2);

    for(size_t i = 0; i < 7; i++){
        TMVP((int32_t*)&res_V_full[i][0], (int32_t*)&src2_Toeplitz_full[i][0], (int32_t*)&src1_V_full[i][0]);
    }

    for(size_t i = 0; i < 4; i++){
        TMVP_TC_Hom_I(des + i, (int32_t*)&res_V_full[0][i], 4);
    }

}

void TMVP_Hom_M(int32_t *des, int32_t *src){
    for(size_t i = 0; i < 32; i++){
        TMVP_TC_Hom_M(des + i * 56, src + i * 16);
    }
}

// ================
// Nussbaumer

static
void shift_addsub(int32_t *des, int32_t *src, size_t jump, size_t shift){

    int32_t buff[16];

    for(size_t i = 0; i < shift; i++){
        buff[i] = src[i] + src[jump * 16 + ( (i + 16 - shift) % 16)];
        des[i]  = src[i] - src[jump * 16 + ( (i + 16 - shift) % 16)];
    }

    for(size_t i = shift; i < 16; i++){
        buff[i] = src[i] - src[jump * 16 + (i - shift)];
        des[i]  = src[i] + src[jump * 16 + (i - shift)];
    }

    for(size_t i = 0; i < 16; i++){
        des[jump * 16 + i] = buff[i];
    }

}

// des and src should not overlap.
static
void CT_1(int32_t *des, int32_t *src){

    for(size_t i = 0; i < 8; i++){
        shift_addsub(des + i * 16, src + i * 16, 8, 0);
    }

    for(size_t i = 16; i < 24; i++){
        shift_addsub(des + i * 16, src + i * 16, 8, 8);
    }

}

// des and src should not overlap.
static
void CT_2(int32_t *des, int32_t *src){

    for(size_t i = 0; i < 4; i++){
        shift_addsub(des + i * 16, src + i * 16, 4,  0);
    }

    for(size_t i = 8; i < 12; i++){
        shift_addsub(des + i * 16, src + i * 16, 4,  8);
    }

    for(size_t i = 16; i < 20; i++){
        shift_addsub(des + i * 16, src + i * 16, 4,  4);
    }

    for(size_t i = 24; i < 28; i++){
        shift_addsub(des + i * 16, src + i * 16, 4, 12);
    }

}

// des and src should not overlap.
static
void CT_3(int32_t *des, int32_t *src){

    for(size_t i = 0; i < 2; i++){
        shift_addsub(des + i * 16, src + i * 16, 2,  0);
    }

    for(size_t i = 4; i < 6; i++){
        shift_addsub(des + i * 16, src + i * 16, 2,  8);
    }

    for(size_t i = 8; i < 10; i++){
        shift_addsub(des + i * 16, src + i * 16, 2,  4);
    }

    for(size_t i = 12; i < 14; i++){
        shift_addsub(des + i * 16, src + i * 16, 2, 12);
    }

    for(size_t i = 16; i < 18; i++){
        shift_addsub(des + i * 16, src + i * 16, 2,  2);
    }

    for(size_t i = 20; i < 22; i++){
        shift_addsub(des + i * 16, src + i * 16, 2, 10);
    }

    for(size_t i = 24; i < 26; i++){
        shift_addsub(des + i * 16, src + i * 16, 2,  6);
    }

    for(size_t i = 28; i < 30; i++){
        shift_addsub(des + i * 16, src + i * 16, 2, 14);
    }

}

// des and src should not overlap.
static
void CT_4(int32_t *des, int32_t *src){

    shift_addsub(des +  0 * 16, src +  0 * 16, 1,  0);
    shift_addsub(des +  2 * 16, src +  2 * 16, 1,  8);
    shift_addsub(des +  4 * 16, src +  4 * 16, 1,  4);
    shift_addsub(des +  6 * 16, src +  6 * 16, 1, 12);
    shift_addsub(des +  8 * 16, src +  8 * 16, 1,  2);
    shift_addsub(des + 10 * 16, src + 10 * 16, 1, 10);
    shift_addsub(des + 12 * 16, src + 12 * 16, 1,  6);
    shift_addsub(des + 14 * 16, src + 14 * 16, 1, 14);
    shift_addsub(des + 16 * 16, src + 16 * 16, 1,  1);
    shift_addsub(des + 18 * 16, src + 18 * 16, 1,  9);
    shift_addsub(des + 20 * 16, src + 20 * 16, 1,  5);
    shift_addsub(des + 22 * 16, src + 22 * 16, 1, 13);
    shift_addsub(des + 24 * 16, src + 24 * 16, 1,  3);
    shift_addsub(des + 26 * 16, src + 26 * 16, 1, 11);
    shift_addsub(des + 28 * 16, src + 28 * 16, 1,  7);
    shift_addsub(des + 30 * 16, src + 30 * 16, 1, 15);

}

static
void addsub_shift_neg(int32_t *des, int32_t *src, size_t jump, size_t shift, bool neg){

    int32_t buff[16];

    if(neg){

        for(size_t i = 0; i < 16 - shift; i++){
            buff[i] = src[jump * 16 + i] - src[i];
            des[i]  = src[jump * 16 + i] + src[i];
        }

        for(size_t i = 16 - shift; i < 16; i++){
            buff[i] = src[i] - src[jump * 16 + i];
            des[i]  = src[i] + src[jump * 16 + i];
        }

    }else{

        for(size_t i = 0; i < 16 - shift; i++){
            buff[i] = src[i] - src[jump * 16 + i];
            des[i]  = src[i] + src[jump * 16 + i];
        }

        for(size_t i = 16 - shift; i < 16; i++){
            buff[i] = src[jump * 16 + i] - src[i];
            des[i]  = src[jump * 16 + i] + src[i];
        }

    }

    for(size_t i = 0; i < 16 - shift; i++){
        des[jump * 16 + i + shift] = buff[i];
    }

    for(size_t i = 16 - shift; i < 16; i++){
        des[jump * 16 + i + shift - 16] = buff[i];
    }

}

// des and src should not overlap.
static
void GS_1(int32_t *des, int32_t *src){

    for(size_t i = 0; i < 8; i++){
        addsub_shift_neg(des + i * 16, src + i * 16, 8, 0, 0);
    }

    for(size_t i = 16; i < 24; i++){
        addsub_shift_neg(des + i * 16, src + i * 16, 8, 8, 1);
    }

}

// des and src should not overlap.
static
void GS_2(int32_t *des, int32_t *src){

    for(size_t i = 0; i < 4; i++){
        addsub_shift_neg(des + i * 16, src + i * 16, 4,  0, 0);
    }

    for(size_t i = 8; i < 12; i++){
        addsub_shift_neg(des + i * 16, src + i * 16, 4,  8, 1);
    }

    for(size_t i = 16; i < 20; i++){
        addsub_shift_neg(des + i * 16, src + i * 16, 4, 12, 1);
    }

    for(size_t i = 24; i < 28; i++){
        addsub_shift_neg(des + i * 16, src + i * 16, 4,  4, 1);
    }

}

// des and src should not overlap.
static
void GS_3(int32_t *des, int32_t *src){

    for(size_t i = 0; i < 2; i++){
        addsub_shift_neg(des + i * 16, src + i * 16, 2,  0, 0);
    }

    for(size_t i = 4; i < 6; i++){
        addsub_shift_neg(des + i * 16, src + i * 16, 2,  8, 1);
    }

    for(size_t i = 8; i < 10; i++){
        addsub_shift_neg(des + i * 16, src + i * 16, 2, 12, 1);
    }

    for(size_t i = 12; i < 14; i++){
        addsub_shift_neg(des + i * 16, src + i * 16, 2,  4, 1);
    }

    for(size_t i = 16; i < 18; i++){
        addsub_shift_neg(des + i * 16, src + i * 16, 2, 14, 1);
    }

    for(size_t i = 20; i < 22; i++){
        addsub_shift_neg(des + i * 16, src + i * 16, 2,  6, 1);
    }

    for(size_t i = 24; i < 26; i++){
        addsub_shift_neg(des + i * 16, src + i * 16, 2, 10, 1);
    }

    for(size_t i = 28; i < 30; i++){
        addsub_shift_neg(des + i * 16, src + i * 16, 2,  2, 1);
    }

}

// des and src should not overlap.
static
void GS_4(int32_t *des, int32_t *src){

    addsub_shift_neg(des +  0 * 16, src +  0 * 16, 1,  0, 0);
    addsub_shift_neg(des +  2 * 16, src +  2 * 16, 1,  8, 1);
    addsub_shift_neg(des +  4 * 16, src +  4 * 16, 1, 12, 1);
    addsub_shift_neg(des +  6 * 16, src +  6 * 16, 1,  4, 1);
    addsub_shift_neg(des +  8 * 16, src +  8 * 16, 1, 14, 1);
    addsub_shift_neg(des + 10 * 16, src + 10 * 16, 1,  6, 1);
    addsub_shift_neg(des + 12 * 16, src + 12 * 16, 1, 10, 1);
    addsub_shift_neg(des + 14 * 16, src + 14 * 16, 1,  2, 1);
    addsub_shift_neg(des + 16 * 16, src + 16 * 16, 1, 15, 1);
    addsub_shift_neg(des + 18 * 16, src + 18 * 16, 1,  7, 1);
    addsub_shift_neg(des + 20 * 16, src + 20 * 16, 1, 11, 1);
    addsub_shift_neg(des + 22 * 16, src + 22 * 16, 1,  3, 1);
    addsub_shift_neg(des + 24 * 16, src + 24 * 16, 1, 13, 1);
    addsub_shift_neg(des + 26 * 16, src + 26 * 16, 1,  5, 1);
    addsub_shift_neg(des + 28 * 16, src + 28 * 16, 1,  9, 1);
    addsub_shift_neg(des + 30 * 16, src + 30 * 16, 1,  1, 1);

}

void Nussbaumer(int32_t *des, int32_t *src){

    for(size_t i = 0; i < 16; i++){
        for(size_t j = 0; j < 16; j++){
            des[(j + 16) * 16 + i] = des[j * 16 + i] = src[i * 16 + j];
        }
    }

    CT_1(des, des);
    CT_2(des, des);
    CT_3(des, des);
    CT_4(des, des);

}

void Nussbaumer_uint16(int32_t *des, uint16_t *src){

    for(size_t i = 0; i < 16; i++){
        for(size_t j = 0; j < 16; j++){
            des[(j + 16) * 16 + i] = des[j * 16 + i] = (int32_t)(uint16_t)src[i * 16 + j];
        }
    }

    CT_1(des, des);
    CT_2(des, des);
    CT_3(des, des);
    CT_4(des, des);

}

void iNussbaumer(int32_t *des, int32_t *src){

    int32_t buff[16];

    // Apply the inverse of symbolic FFT.
    GS_4(src, src);
    GS_3(src, src);
    GS_2(src, src);
    GS_1(src, src);

// ========

    // Apply the 0th layer of GS butterfly buffly and
    // map
    // ( Z_{32}[y] / (y^16 + 1) ) [x] / (x^32 - 1)
    // to
    // ( Z_{32}[y] / (y^16 + 1) ) [x] / (x^16 - y)
    // on the fly.

    for(size_t i = 0; i < 16; i++){
        des[i]  = src[i * 16] + src[(i + 16) * 16];
        buff[i] = src[i * 16] - src[(i + 16) * 16];
    }

    for(size_t i = 1; i < 16; i++){
        for(size_t j = 0; j < 16; j++){
            des[i * 16 + j] = src[j * 16 + i] + src[(j + 16) * 16 + i];
            des[i * 16 + j] = des[i * 16 + j] + buff[j];
            des[i * 16 + j] = des[i * 16 + j] >> 5;
            buff[j]         = src[j * 16 + i] - src[(j + 16) * 16 + i];
        }
    }

    for(size_t i = 0; i < 16; i++){
        des[i] -= buff[i];
        des[i] >>= 5;
    }

}

void iNussbaumer_uint16(uint16_t *des, int32_t *src){

    int32_t buff[16];
    int32_t des_buff[16];
    int32_t t;

    // Apply the inverse of symbolic FFT.
    GS_4(src, src);
    GS_3(src, src);
    GS_2(src, src);
    GS_1(src, src);

// ========

    // Apply the 0th layer of GS butterfly buffly and
    // map
    // ( Z_{32}[y] / (y^16 + 1) ) [x] / (x^32 - 1)
    // to
    // ( Z_{32}[y] / (y^16 + 1) ) [x] / (x^16 - y)
    // on the fly.

    for(size_t i = 0; i < 16; i++){
        des_buff[i] = src[i * 16] + src[(i + 16) * 16];
        buff[i]     = src[i * 16] - src[(i + 16) * 16];
    }

    for(size_t i = 1; i < 16; i++){
        for(size_t j = 0; j < 16; j++){
            t               = src[j * 16 + i] + src[(j + 16) * 16 + i];
            t               = t + buff[j];
            t               = t >> 5;
            des[i * 16 + j] = (uint16_t)t;
            buff[j]         = src[j * 16 + i] - src[(j + 16) * 16 + i];
        }
    }

    for(size_t i = 0; i < 16; i++){
        t      = des_buff[i] - buff[i];
        t      = t >> 5;
        des[i] = (uint16_t)t;
    }

}

void iNussbaumer_acc_uint16(uint16_t *des, int32_t *src){

    int32_t buff[16];
    int32_t des_buff[16];
    int32_t t;

    // Apply the inverse of symbolic FFT.
    GS_4(src, src);
    GS_3(src, src);
    GS_2(src, src);
    GS_1(src, src);

// ========

    // Apply the 0th layer of GS butterfly buffly and
    // map
    // ( Z_{32}[y] / (y^16 + 1) ) [x] / (x^32 - 1)
    // to
    // ( Z_{32}[y] / (y^16 + 1) ) [x] / (x^16 - y)
    // on the fly.

    for(size_t i = 0; i < 16; i++){
        des_buff[i] = src[i * 16] + src[(i + 16) * 16];
        buff[i]     = src[i * 16] - src[(i + 16) * 16];
    }

    for(size_t i = 1; i < 16; i++){
        for(size_t j = 0; j < 16; j++){
            t               = src[j * 16 + i] + src[(j + 16) * 16 + i];
            t               = t + buff[j];
            t               = t >> 5;
            des[i * 16 + j] += (uint16_t)t;
            buff[j]         = src[j * 16 + i] - src[(j + 16) * 16 + i];
        }
    }

    for(size_t i = 0; i < 16; i++){
        t      = des_buff[i] - buff[i];
        t      = t >> 5;
        des[i] += (uint16_t)t;
    }

}

void BiHom(int32_t *des, int32_t *src_M, int32_t *src_V){

    for(size_t i = 0; i < 32; i++){
        TMVP_TC4_negacyclic16_mul(des + i * 16, src_M + i * 16, src_V + i * 16);
    }

}

