#include "poly.h"
#include "cbd.h"
#include "fips202.h"
#include "pack_unpack.h"
#include "Hom.h"

#define h1 (1 << (SABER_EQ - SABER_EP - 1))
#define h2 ((1 << (SABER_EP - 2)) - (1 << (SABER_EP - SABER_ET - 1)) + (1 << (SABER_EQ - SABER_EP - 1)))
#define MAX(a,b) (((a)>(b))?(a):(b))

static
void polymul(uint16_t des[SABER_N], uint16_t poly1[SABER_N], uint16_t poly2[SABER_N]){

    int32_t poly1_Hom[32 * 16];
    int32_t poly2_Hom[32 * 16];

    Nussbaumer_uint16(poly1_Hom, poly1);
    Nussbaumer_uint16(poly2_Hom, poly2);

    BiHom(poly1_Hom, poly2_Hom, poly1_Hom);

    iNussbaumer_uint16(des, poly1_Hom);

}

static
void polymla(uint16_t des[SABER_N], uint16_t poly1[SABER_N], uint16_t poly2[SABER_N]){

    int32_t poly1_Hom[32 * 16];
    int32_t poly2_Hom[32 * 16];

    Nussbaumer_uint16(poly1_Hom, poly1);
    Nussbaumer_uint16(poly2_Hom, poly2);

    BiHom(poly1_Hom, poly2_Hom, poly1_Hom);

    iNussbaumer_acc_uint16(des, poly1_Hom);

}

static inline shake128incctx shake128_absorb_seed(const uint8_t seed[SABER_SEEDBYTES]){

    shake128incctx ctx;
    shake128_inc_init(&ctx);
    shake128_inc_absorb(&ctx, seed, SABER_SEEDBYTES);
    shake128_inc_finalize(&ctx);

    return ctx;

}

void MatrixVectorMulKeyPairNTT_B( uint8_t pk[SABER_INDCPA_PUBLICKEYBYTES], uint8_t sk[SABER_INDCPA_SECRETKEYBYTES]){

    uint16_t poly[SABER_N];
    uint16_t acc[SABER_L * SABER_N];

    int32_t s_Hom[32 * 7 * 8];
    int32_t poly_Hom[32 * 16];
    int32_t buff_V[7 * 4];

    uint8_t shake_out[MAX(SABER_POLYBYTES, SABER_POLYCOINBYTES)];

    uint8_t *seed_A = pk + SABER_POLYVECCOMPRESSEDBYTES;
    uint8_t *seed_s = sk;

    size_t i, j;

    shake128incctx shake_s_ctx = shake128_absorb_seed(seed_s);
    shake128incctx shake_A_ctx = shake128_absorb_seed(seed_A);

    for (i = 0; i < SABER_L; i++) {

        shake128_inc_squeeze(shake_out, SABER_POLYCOINBYTES, &shake_s_ctx);
        cbd(poly, shake_out);
#ifdef SABER_COMPRESS_SECRETKEY
        POLmu2BS(sk + i * SABER_POLYSECRETBYTES, poly); // sk <- s
#else
        POLq2BS(sk + i * SABER_POLYSECRETBYTES, poly);
#endif

        Nussbaumer_uint16(poly_Hom, poly);
        TMVP_Hom_M(s_Hom, poly_Hom);

        for (j = 0; j < SABER_L; j++) {

            shake128_inc_squeeze(shake_out, SABER_POLYBYTES, &shake_A_ctx);
            BS2POLq(shake_out, poly);

            Nussbaumer_uint16(poly_Hom, poly);

            for(size_t k = 0; k < 32; k++){

                for(size_t h = 0; h < 4; h++){
                    TMVP_TC_Hom_V(buff_V + h, poly_Hom + k * 16 + h, 4);
                }

                for(size_t h = 0; h < 7; h++){
                    TMVP(buff_V + h * 4, s_Hom + k * 56 + h * 8, buff_V + h * 4);
                }

                for(size_t h = 0; h < 4; h++){
                    TMVP_TC_Hom_I(poly_Hom + k * 16 + h , buff_V + h, 4);
                }

            }

            if (i == 0) {
                iNussbaumer_uint16(acc + j * SABER_N, poly_Hom);
            } else {
                iNussbaumer_acc_uint16(acc + j * SABER_N, poly_Hom);
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

uint32_t MatrixVectorMulEncNTT_C(uint8_t ct0[SABER_POLYVECCOMPRESSEDBYTES], uint8_t ct1[SABER_SCALEBYTES_KEM], const uint8_t seed_s[SABER_NOISE_SEEDBYTES], const uint8_t seed_A[SABER_SEEDBYTES], const uint8_t pk[SABER_INDCPA_PUBLICKEYBYTES], const uint8_t m[SABER_KEYBYTES], int compare){

    uint16_t poly[SABER_N];
    uint16_t s_poly[SABER_L * SABER_N];
    uint16_t acc[SABER_N];

    int32_t s_Hom[32 * 16];
    int32_t poly_Hom[32 * 16];

    int32_t acc_Hom_full[32 * 7 * 4];

    int32_t s_TC_buff[7 * 4];
    int32_t poly_TC_buff[7 * 4];

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

            Nussbaumer_uint16(s_Hom, s_poly + j * SABER_N);
            Nussbaumer_uint16(poly_Hom, poly);

            if (j == 0) {

                for(size_t k = 0; k < 32; k++){
                    TC(s_TC_buff, s_Hom + k * 16);
                    TC(poly_TC_buff, poly_Hom + k * 16);
                    for(size_t h = 0; h < 7; h++){
                        negacyclic_4x4(acc_Hom_full + (k * 7 + h) * 4, s_TC_buff + h * 4, poly_TC_buff + h * 4);
                    }
                }

            } else {

                for(size_t k = 0; k < 32; k++){
                    TC(s_TC_buff, s_Hom + k * 16);
                    TC(poly_TC_buff, poly_Hom + k * 16);
                    for(size_t h = 0; h < 7; h++){
                        negacyclic_acc_4x4(acc_Hom_full + (k * 7 + h) * 4, s_TC_buff + h * 4, poly_TC_buff + h * 4);
                    }
                }

            }

        }

        iTC_Hom(poly_Hom, acc_Hom_full);
        iNussbaumer_uint16(acc, poly_Hom);

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

        Nussbaumer_uint16(s_Hom, s_poly + j * SABER_N);
        Nussbaumer_uint16(poly_Hom, poly);

        if (j == 0) {

            for(size_t k = 0; k < 32; k++){
                TC(s_TC_buff, s_Hom + k * 16);
                TC(poly_TC_buff, poly_Hom + k * 16);
                for(size_t h = 0; h < 7; h++){
                    negacyclic_4x4(acc_Hom_full + (k * 7 + h) * 4, s_TC_buff + h * 4, poly_TC_buff + h * 4);
                }
            }

        } else {

            for(size_t k = 0; k < 32; k++){
                TC(s_TC_buff, s_Hom + k * 16);
                TC(poly_TC_buff, poly_Hom + k * 16);
                for(size_t h = 0; h < 7; h++){
                    negacyclic_acc_4x4(acc_Hom_full + (k * 7 + h) * 4, s_TC_buff + h * 4, poly_TC_buff + h * 4);
                }
            }

        }

    }

    iTC_Hom(poly_Hom, acc_Hom_full);
    iNussbaumer_uint16(acc, poly_Hom);

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




