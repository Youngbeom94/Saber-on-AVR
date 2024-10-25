#include "poly.h"
#include "cbd.h"
#include "fips202.h"
#include <string.h>
#include "pack_unpack.h"

#define h1 (1 << (SABER_EQ - SABER_EP - 1))
#define h2 ((1 << (SABER_EP - 2)) - (1 << (SABER_EP - SABER_ET - 1)) + (1 << (SABER_EQ - SABER_EP - 1)))
#define MAX(a,b) (((a)>(b))?(a):(b))
// #define OVERFLOWING_MUL(X, Y) (X * Y)

// Scaling for the inversion.
uint16_t iTC4_modified_scale[7] = 
{
//1, 20389, -21845, -29127, -21845, 20389, 1
1, 0x4FA5, 0xAAAB, 0x8E39, 0xAAAB, 0x4FA5, 1
};

// Multiply two size-64 polynomials in R[x] / (x^64 + 1) via
// R[x] / (x^64 + 1)
// to
// (R[y] / (y^16 + 1)) / (x^4 - y)
// to
// (R[y] / (y^16 + 1)) / (x^7) (2 layers of Karatsuba)
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

static
void TC_striding_mul(uint16_t *des, uint16_t *src1, uint16_t *src2){

    uint16_t src1_extended[7][64], src2_extended[7][64];
    uint16_t res[7][64];
    uint16_t TC4_buff[7];
    // uint16_t twiddle;
    uint16_t t0, t1, t2, t3;

    // Apply Toom-4 evaluation matrix.
    for(size_t i = 0; i < 64; i++){
        src1_extended[0][i] = src1[i * 4 + 0];
        t0 = src1[i * 4 + 0] + src1[i * 4 + 2];
        t1 = src1[i * 4 + 1] + src1[i * 4 + 3];
        src1_extended[1][i] = t0 + t1;
        src1_extended[2][i] = t0 - t1;
        t2 = src1[i * 4 + 0] + 4 * src1[i * 4 + 2];
        t3 = src1[i * 4 + 1] + 4 * src1[i * 4 + 3];
        src1_extended[3][i] = t2 + 2 * t3;
        src1_extended[4][i] = t2 - 2 * t3;
        src1_extended[5][i] = ((src1[i * 4 + 0] * 2 + src1[i * 4 + 1]) * 2 + src1[i * 4 + 2] ) * 2 + src1[i * 4 + 3];
        src1_extended[6][i] = src1[i * 4 + 3];
    }

    for(size_t i = 0; i < 64; i++){
        src2_extended[0][i] = src2[i * 4 + 0];
        t0 = src2[i * 4 + 0] + src2[i * 4 + 2];
        t1 = src2[i * 4 + 1] + src2[i * 4 + 3];
        src2_extended[1][i] = t0 + t1;
        src2_extended[2][i] = t0 - t1;
        t2 = src2[i * 4 + 0] + 4 * src2[i * 4 + 2];
        t3 = src2[i * 4 + 1] + 4 * src2[i * 4 + 3];
        src2_extended[3][i] = t2 + 2 * t3;
        src2_extended[4][i] = t2 - 2 * t3;
        src2_extended[5][i] = ((src2[i * 4 + 0] * 2 + src2[i * 4 + 1]) * 2 + src2[i * 4 + 2] ) * 2 + src2[i * 4 + 3];
        src2_extended[6][i] = src2[i * 4 + 3];
    }

    // Compute small-dimensional products.
    for(size_t i = 0; i < 7; i++){
        negacyclic_Karatsuba_striding((uint16_t*)&res[i][0],
                                      (uint16_t*)&src1_extended[i][0],
                                      (uint16_t*)&src2_extended[i][0]);
    }

    // Apply Toom-4 inversion matrix.
    for(size_t i = 0; i < 64; i++){

        // {-360, -120,  -40,    5,     3,    8, -360},
        t0 = ( (res[0][i] + res[6][i]) * 3 + res[1][i]) * 3 + res[2][i];
        t0 = res[5][i] - t0 * 5;
        t0 *= 8;
        t0 += res[3][i] * 5;
        t0 += res[4][i] * 3;
        TC4_buff[1] = t0 * iTC4_modified_scale[1];

        // { -30,   16,   16,   -1,    -1,    0,   96},
        t0 = res[1][i] + res[2][i];
        t1 = res[3][i] + res[4][i];
        t0 = 16 * t0 - t1;
        t0 = t0 - 30 * res[0][i];
        t0 += 96 * res[6][i];
        TC4_buff[2] = t0 * iTC4_modified_scale[2];

        // {  45,   27,   -7,   -1,     0,   -1,   45},
        t0 = res[0][i] + res[6][i];
        t1 = res[3][i] + res[5][i];
        t0 = 45 * t0 - t1;
        t0 += 27 * res[1][i];
        t0 = t0 - 7 * res[2][i];
        TC4_buff[3] = t0 * iTC4_modified_scale[3];

        // {   6,   -4,   -4,    1,     1,    0, -120},
        t0 = res[1][i] + res[2][i];
        t1 = res[3][i] + res[4][i];
        t0 = t1 - 4 * t0;
        t0 += 6 * res[0][i];
        t0 = t0 - 120 * res[6][i];
        TC4_buff[4] = t0 * iTC4_modified_scale[4];

        // { -90,  -60,   20,    5,    -3,    2,  -90},
        t0 = res[0][i] + res[6][i];
        t0 *= -90;
        t0 = t0 - 60 * res[1][i];
        t0 += 20 * res[2][i];
        t0 += 5 * res[3][i];
        t0 = t0 - 3 * res[4][i];
        t0 += 2 * res[5][i];
        TC4_buff[5] = t0 * iTC4_modified_scale[5];

        // Multiply by powers of two.
        res[1][i] = TC4_buff[1] >> 2;
        res[2][i] = TC4_buff[2] >> 3;
        res[3][i] = TC4_buff[3] >> 1;
        res[4][i] = TC4_buff[4] >> 3;
        res[5][i] = TC4_buff[5] >> 2;
    }

    // Export the result.
    for(size_t i = 0; i < 64; i++){
        des[i * 4 + 3] = res[3][i];
    }

    for(size_t i = 0; i < 63; i++){
        des[(i + 1) * 4 + 0] = res[0][i + 1] + res[4][i];
        des[(i + 1) * 4 + 1] = res[1][i + 1] + res[5][i];
        des[(i + 1) * 4 + 2] = res[2][i + 1] + res[6][i];
    }

    des[0] = res[0][0] - res[4][63];
    des[1] = res[1][0] - res[5][63];
    des[2] = res[2][0] - res[6][63];
}

static
void polymul(uint16_t des[SABER_N], uint16_t poly1[SABER_N], uint16_t poly2[SABER_N])
{
    TC_striding_mul(des, poly1, poly2);
}

static
void polymla(uint16_t des[SABER_N], uint16_t poly1[SABER_N], uint16_t poly2[SABER_N])
{
    TC_striding_mul(poly2, poly1, poly2);
    for(int cnt_i = 0 ; cnt_i < SABER_N ; cnt_i ++) des[cnt_i] += poly2[cnt_i];
}

static inline shake128incctx shake128_absorb_seed(const uint8_t seed[SABER_SEEDBYTES]){

    shake128incctx ctx;
    shake128_inc_init(&ctx);
    shake128_inc_absorb(&ctx, seed, SABER_SEEDBYTES);
    shake128_inc_finalize(&ctx);

    return ctx;

}

void MatrixVectorMulKeyPairNTT_D( uint8_t pk[SABER_INDCPA_PUBLICKEYBYTES], uint8_t sk[SABER_INDCPA_SECRETKEYBYTES]){

    uint16_t poly[SABER_N];
    uint16_t s_poly[SABER_N];
    uint16_t acc[SABER_L * SABER_N];

    uint8_t shake_out[MAX(SABER_POLYBYTES, SABER_POLYCOINBYTES)];

    uint8_t *seed_A = pk + SABER_POLYVECCOMPRESSEDBYTES;
    uint8_t *seed_s = sk;

    size_t i, j;

    shake128incctx shake_s_ctx = shake128_absorb_seed(seed_s);
    shake128incctx shake_A_ctx = shake128_absorb_seed(seed_A);

    for (i = 0; i < SABER_L; i++) {

        shake128_inc_squeeze(shake_out, SABER_POLYCOINBYTES, &shake_s_ctx);
        cbd(s_poly, shake_out);
#ifdef SABER_COMPRESS_SECRETKEY
        POLmu2BS(sk + i * SABER_POLYSECRETBYTES, s_poly); // sk <- s
#else
        POLq2BS(sk + i * SABER_POLYSECRETBYTES, s_poly);
#endif

        for (j = 0; j < SABER_L; j++) {

            shake128_inc_squeeze(shake_out, SABER_POLYBYTES, &shake_A_ctx);
            BS2POLq(shake_out, poly);

            if (i == 0) {
                polymul(acc + j * SABER_N, s_poly, poly);
            } else {
                polymla(acc + j * SABER_N, s_poly, poly);
            }

        }
    }

    shake128_inc_ctx_release(&shake_A_ctx);
    shake128_inc_ctx_release(&shake_s_ctx);

    for (i = 0; i < SABER_L; i++) {

        for (j = 0; j < SABER_N; j++) {
            poly[j] = ((acc[i * SABER_N + j] + h1) >> (SABER_EQ - SABER_EP));
        }

        POLp2BS(pk + i * SABER_POLYCOMPRESSEDBYTES, poly);
    }

}

uint32_t MatrixVectorMulEncNTT_D(uint8_t ct0[SABER_POLYVECCOMPRESSEDBYTES], uint8_t ct1[SABER_SCALEBYTES_KEM], const uint8_t seed_s[SABER_NOISE_SEEDBYTES], const uint8_t seed_A[SABER_SEEDBYTES], const uint8_t pk[SABER_INDCPA_PUBLICKEYBYTES], const uint8_t m[SABER_KEYBYTES], int compare){

    uint16_t poly[SABER_N];
    uint16_t s_poly[SABER_L * SABER_N];
    uint16_t acc[SABER_N];

    uint8_t shake_out[MAX(SABER_POLYBYTES, SABER_POLYCOINBYTES)];

    uint16_t *mp = poly;

    size_t i, j;
    uint32_t fail = 0;

    shake128incctx shake_s_ctx = shake128_absorb_seed(seed_s);

    for(i = 0; i < SABER_L; i++){
        shake128_inc_squeeze(shake_out, SABER_POLYCOINBYTES, &shake_s_ctx);
        cbd(s_poly + i * SABER_N, shake_out);
    }

    shake128_inc_ctx_release(&shake_s_ctx);

    shake128incctx shake_A_ctx = shake128_absorb_seed(seed_A);

    for (i = 0; i < SABER_L; i++) {

        for (j = 0; j < SABER_L; j++) {

            shake128_inc_squeeze(shake_out, SABER_POLYBYTES, &shake_A_ctx);
            BS2POLq(shake_out, poly);

            if (j == 0) {
                polymul(acc, s_poly + j * SABER_N, poly);
            } else {
                polymla(acc, s_poly + j * SABER_N, poly);
            }

        }

        for (j = 0; j < SABER_N; j++) {
            acc[j] = ((acc[j] + h1) >> (SABER_EQ - SABER_EP));
        }

        if (compare) {
            fail |= POLp2BS_cmp(ct0 + i * SABER_POLYCOMPRESSEDBYTES, acc);
        } else {
            POLp2BS(ct0 + i * SABER_POLYCOMPRESSEDBYTES, acc);
        }
    }

    shake128_inc_ctx_release(&shake_A_ctx);

    for(j = 0; j < SABER_L; j++){

        BS2POLp(pk + j * SABER_POLYCOMPRESSEDBYTES, poly);

        if(j == 0){
            polymul(acc, s_poly + j * SABER_N, poly);
        }else{
            polymla(acc, s_poly + j * SABER_N, poly);
        }

    }

    BS2POLmsg(m, mp);

    for(j = 0; j < SABER_N; j++){
        acc[j] = (acc[j] - (mp[j] << (SABER_EP - 1)) + h1) >> (SABER_EP - SABER_ET);
    }

    if(compare){
        fail |= POLT2BS_cmp(ct1, acc);
    }else{
        POLT2BS(ct1, acc);
    }

    return fail;

}



void InnerProdDecNTT(uint8_t m[SABER_KEYBYTES], const uint8_t ciphertext[SABER_BYTES_CCA_DEC], const uint8_t sk[SABER_INDCPA_SECRETKEYBYTES]){

    uint16_t poly[SABER_N];
    uint16_t buff[SABER_N];
    uint16_t acc[SABER_N];

    size_t i;

    for (i = 0; i < SABER_L; i++) {

#ifdef SABER_COMPRESS_SECRETKEY
        BS2POLmu(sk + i * SABER_POLYSECRETBYTES, buff);
#else
        BS2POLq(sk + i * SABER_POLYSECRETBYTES, buff);
#endif
        BS2POLp(ciphertext + i * SABER_POLYCOMPRESSEDBYTES, poly);

        if(i == 0){
            polymul(acc, buff, poly);
        }else{
            polymla(acc, buff, poly);
        }

    }

    BS2POLT(ciphertext + SABER_POLYVECCOMPRESSEDBYTES, buff);

    for (i = 0; i < SABER_N; i++) {
        poly[i] = (acc[i] + h2 - (buff[i] << (SABER_EP - SABER_ET))) >> (SABER_EP - 1);
    }

    POLmsg2BS(m, poly);

}




