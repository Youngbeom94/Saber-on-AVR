#include "poly.h"
#include "cbd.h"
#include "fips202.h"
#include <string.h>
#include "pack_unpack.h"

#define h1 (1 << (SABER_EQ - SABER_EP - 1))
#define h2 ((1 << (SABER_EP - 2)) - (1 << (SABER_EP - SABER_ET - 1)) + (1 << (SABER_EQ - SABER_EP - 1)))
#define MAX(a,b) (((a)>(b))?(a):(b))


// Multiply two size-64 polynomials in R[x] / (x^64 + 1) via
// R[x] / (x^64 + 1)
// to
// (R[y] / (y^16 + 1)) / (x^4 - y)
// to
// (R[y] / (y^16 + 1)) / (x^7) (2 layers of Karatsuba)

uint16_t res_p00[16], res_p01[16], res_p02[16], res_p10[16], res_p11[16], res_p12[16], res_p20[16], res_p21[16], res_p22[16];


extern void karatsuba_striding_asm(uint16_t *src1_des, uint16_t *src2);
extern void TC_striding_asm(uint16_t src_ex[7][64], uint16_t src[256]);
extern void iTC_striding_asm(uint16_t src[7][64]);

static
void negacyclic_Karatsuba_striding(uint16_t *src1_des, uint16_t *src2)
{
    memset(res_p00, 0, sizeof(res_p00));
    memset(res_p01, 0, sizeof(res_p01));
    memset(res_p02, 0, sizeof(res_p02));
    memset(res_p10, 0, sizeof(res_p10));
    memset(res_p11, 0, sizeof(res_p11));
    memset(res_p12, 0, sizeof(res_p12));
    memset(res_p20, 0, sizeof(res_p20));
    memset(res_p21, 0, sizeof(res_p21));
    memset(res_p22, 0, sizeof(res_p22));
    
    karatsuba_striding_asm(src1_des, src2);

}

static
void TC_striding_mul(uint16_t *des, uint16_t *src1, uint16_t *src2)
{    
    uint16_t src1_ex[7][64], src2_ex[7][64];    

    // Apply Toom-4 evaluation matrix.    
    TC_striding_asm(src1_ex, src1);
    TC_striding_asm(src2_ex, src2);
    
    // Compute small-dimensional products.
    for(size_t i = 0; i < 7; i++)   negacyclic_Karatsuba_striding((uint16_t*)&src1_ex[i][0], (uint16_t*)&src2_ex[i][0]);        

    // Apply Toom-4 inversion matrix.
    iTC_striding_asm(src1_ex);

    // Export the result.
    for(size_t i = 0; i < 64; i++)
    {
        des[i * 4 + 3] = src1_ex[3][i];
    }

    for(size_t i = 0; i < 63; i++)
    {
        des[(i + 1) * 4 + 0] = src1_ex[0][i + 1] + src1_ex[4][i];
        des[(i + 1) * 4 + 1] = src1_ex[1][i + 1] + src1_ex[5][i];
        des[(i + 1) * 4 + 2] = src1_ex[2][i + 1] + src1_ex[6][i];
    }

    des[0] = src1_ex[0][0] - src1_ex[4][63];
    des[1] = src1_ex[1][0] - src1_ex[5][63];
    des[2] = src1_ex[2][0] - src1_ex[6][63];
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




