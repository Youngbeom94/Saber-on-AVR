#include "poly.h"
#include "cbd.h"
#include "fips202.h"
#include "pack_unpack.h"

#define h1 (1 << (SABER_EQ - SABER_EP - 1))
#define h2 ((1 << (SABER_EP - 2)) - (1 << (SABER_EP - SABER_ET - 1)) + (1 << (SABER_EQ - SABER_EP - 1)))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define SABER_NTT_Q 26017793
#define SABER_NTT_QINV (-26017791)  // q^(-1) mod 2^32,  ( 4268949505 )
#define MONT 2031451                // 2^32 mod q,  ( -23986342 )

__flash
static const int32_t zetas[SABER_N] = {
    0, -12801892, 11254964, 10461102, 10834930, 3884520, -1243096, 4242387,
    -11878620, -1411302, -7034539, -9865511, -1089702, 10312081, 7234003, 2783608,
    3379100, 2996851, -10660734, 5016187, 1964411, 10215805, -8852218, -5777903,
    3110445, 11217276, 9064857, -2328877, -12612507, 10666144, 3417718, 8537670,
    157609, 5761559, 4376943, 10294824, -7007758, 2435542, 11505619, -5146905,
    -7274932, 4189852, 2951046, 9722261, -8812312, 6489555, -9407632, 7418655,
    10138881, 541874, 3865161, 394839, 2842122, -3795799, -4954840, 10501964,
    12446017, 6659183, 10286831, 7349045, 12810998, 501589, -7840175, 6191705,
    1, 12468715, -9505664, -9735085, 4638815, 10625597, -9136174, 7110376,
    10929114, -7784561, -10444142, 8039860, -3390925, -2718795, -7746019, 9647706,
    -7650810, 3675137, 1446934, -6459008, -2883401, -8309560, -11094930, -741134,
    6870299, 8600941, -12020096, 7048551, -7177184, -7369413, -9270010, -605137,
    -6395434, 541145, 9845685, 9290957, -4399979, -3178844, 968457, 10217802,
    9379922, 779735, -8676447, -2084786, 8161711, 795579, 533389, 7178475,
    -7932738, -12195153, 3554747, -4380905, -4335783, 1047651, -3320665, 4089002,
    -7209882, -10596173, -1623716, -10493162, -8259776, 12885567, 6919167, -6946499,
    -65536, -8881489, -6839688, -9789386, 8531365, 7104653, 828955, -6928906,
    -6591607, 12104552, -8808132, 12078876, 9690787, 9102656, 9941961, 12345070,
    -11422536, -7068631, 8588861, -11943822, -662623, -3101123, -1928491, -4261707,
    12010394, 4215969, 8292795, 12076379, -11749023, -6441091, 3908810, 7141900,
    10535187, -2226861, -5545760, 2251627, 1823925, 4251833, -11400825, 12084362,
    -172981, -1767508, 764577, 9104253, -12103602, 591828, 11732288, 5195426,
    -7622158, 6981634, -580870, 644325, 9557335, 2099791, 10280788, 6432828,
    -2311921, -10119235, -921594, 3578049, 12496571, -9011511, 9585685, 12434343,
    397, 6699185, -1168623, 11822412, -5653748, 3479543, -10587851, 12897628,
    -6113173, 5646650, -9495287, -8364119, 6728011, -12632102, -5069969, 5523711,
    6710211, 2032981, 2041352, 11535331, 72695, 5364391, -7680193, -8034475,
    -4359562, 6242694, -10721993, -11646897, 12615182, -11664145, -11685157, -6079252,
    10756416, 6692221, 6067995, -6016677, -3599532, 12870789, -5789466, -2308314,
    3284635, -2658721, -10200783, 4909334, -12024858, 3631347, 3613089, -12102655,
    -1144033, -2166243, 6273737, 3972846, -4131513, -367241, 8603438, 10230628,
    -365924, 8201785, 5829573, -2938434, -889154, -9935122, -10976759, 125955
};

int32_t montgomery_reduce(int64_t a) {
    int32_t t;

    t = (int32_t)((uint64_t)a * (uint64_t)SABER_NTT_QINV);
    t = (a - (int64_t)t * SABER_NTT_Q) >> 32;
    return t;
}


void ntt(int32_t a[SABER_N]) {
    unsigned int len, start, j, k;
    int32_t zeta, t;

    k = 0;
    for (len = 128; len > 0; len >>= 1) {
        for (start = 0; start < SABER_N; start = j + len) {
            zeta = zetas[++k];
            for (j = start; j < start + len; ++j) {
                t = montgomery_reduce((int64_t)zeta * a[j + len]);
                a[j + len] = (a[j] - t) % SABER_NTT_Q;
                a[j] = (a[j] + t)% SABER_NTT_Q;
            }
        }
    }
}

void invntt_tomont(int32_t a[SABER_N]) 
{
    unsigned int start, len, j, k;
    int32_t t, zeta;
    // const int32_t f = 41978; // mont^2/256
    const int32_t f = 6226687; // mont^2/256 mod q : SABER
    k = 256;
    for (len = 1; len < SABER_N; len <<= 1) {
        for (start = 0; start < SABER_N; start = j + len) {
            zeta = -zetas[--k];
            for (j = start; j < start + len; ++j) 
            {
                t = a[j];
                a[j] = (t + a[j + len])% SABER_NTT_Q;
                a[j + len] = (t - a[j + len])% SABER_NTT_Q;
                a[j + len] = montgomery_reduce((int64_t)zeta * a[j + len]);
            }
        }
    }

    for (j = 0; j < SABER_N; ++j) {
        a[j] = montgomery_reduce((int64_t)f * a[j]);
    }
}

static inline
void polymul(uint16_t des[SABER_N], uint16_t poly1[SABER_N], uint16_t poly2[SABER_N])
{

    int32_t a_tmp[SABER_N] = { 0 };
    int32_t b_tmp[SABER_N] = { 0 };
    int32_t res_tmp[SABER_N] = { 0 };

    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) 
    {
        a_tmp[cnt_i] = (int32_t)(int16_t)poly1[cnt_i];
        b_tmp[cnt_i] = (int32_t)poly2[cnt_i];
    }

    ntt(a_tmp);
    ntt(b_tmp);
    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) res_tmp[cnt_i] = montgomery_reduce((int64_t)a_tmp[cnt_i] * b_tmp[cnt_i]);
    invntt_tomont(res_tmp);
    
    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++)  des[cnt_i] = (uint16_t)((res_tmp[cnt_i] & 8191));  
}

static inline
void polymul_dec(uint16_t des[SABER_N], uint16_t poly1[SABER_N], uint16_t poly2[SABER_N])
{

    int32_t a_tmp[SABER_N] = { 0 };
    int32_t b_tmp[SABER_N] = { 0 };
    int32_t res_tmp[SABER_N] = { 0 };

    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) 
    {
        poly1[cnt_i] = poly1[cnt_i] & 8191;
        
        a_tmp[cnt_i] = (int32_t)(int16_t)poly1[cnt_i];

        if(a_tmp[cnt_i] >= 8192/2) a_tmp[cnt_i] -= 8192;
        if(a_tmp[cnt_i] < -8192/2) a_tmp[cnt_i] += 8192;
        b_tmp[cnt_i] = (int32_t)(poly2[cnt_i]);
    }

    ntt(a_tmp);
    ntt(b_tmp);

    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) res_tmp[cnt_i] = montgomery_reduce((int64_t)a_tmp[cnt_i] * b_tmp[cnt_i]);
    invntt_tomont(res_tmp);

    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) des[cnt_i] = (uint16_t)((res_tmp[cnt_i] & 8191));        
}

static inline
void polymla(uint16_t des[SABER_N], uint16_t poly1[SABER_N], uint16_t poly2[SABER_N])
{

    int32_t a_tmp[SABER_N] = { 0 };
    int32_t b_tmp[SABER_N] = { 0 };
    int32_t res_tmp[SABER_N] = { 0 };

    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) 
    {
        a_tmp[cnt_i] = (int32_t)(int16_t)poly1[cnt_i];
        b_tmp[cnt_i] = (int32_t)poly2[cnt_i];
    }

    ntt(a_tmp);
    ntt(b_tmp);

    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) res_tmp[cnt_i] = montgomery_reduce((int64_t)a_tmp[cnt_i] * b_tmp[cnt_i]);
    invntt_tomont(res_tmp);

    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) des[cnt_i] = (des[cnt_i]+ (uint16_t)((res_tmp[cnt_i] & 8191))) & 8191;
}

static inline
void polymla_dec(uint16_t des[SABER_N], uint16_t poly1[SABER_N], uint16_t poly2[SABER_N])
{

    int32_t a_tmp[SABER_N] = { 0 };
    int32_t b_tmp[SABER_N] = { 0 };
    int32_t res_tmp[SABER_N] = { 0 };

    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) 
    {
        poly1[cnt_i] = poly1[cnt_i] & 8191;        
        a_tmp[cnt_i] = (int32_t)(int16_t)poly1[cnt_i];
        if(a_tmp[cnt_i] >= 8192/2) a_tmp[cnt_i] -= 8192;
        if(a_tmp[cnt_i] < -8192/2) a_tmp[cnt_i] += 8192;
        b_tmp[cnt_i] = (int32_t)(poly2[cnt_i]);
    }

    ntt(a_tmp);
    ntt(b_tmp);

    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) res_tmp[cnt_i] = montgomery_reduce((int64_t)a_tmp[cnt_i] * b_tmp[cnt_i]);

    invntt_tomont(res_tmp);

    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) des[cnt_i] += (uint16_t)((res_tmp[cnt_i] & 8191));
}

static inline shake128incctx shake128_absorb_seed(const uint8_t seed[SABER_SEEDBYTES])
{

    shake128incctx ctx;
    shake128_inc_init(&ctx);
    shake128_inc_absorb(&ctx, seed, SABER_SEEDBYTES);
    shake128_inc_finalize(&ctx);

    return ctx;
}

void MatrixVectorMulKeyPairNTT_D( uint8_t pk[SABER_INDCPA_PUBLICKEYBYTES], uint8_t sk[SABER_INDCPA_SECRETKEYBYTES])
{

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

            if (j == 0) 
            {
                polymul(acc, s_poly + j * SABER_N, poly);
            } 
            else 
            {
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

        if(j == 0)
        {
            polymul(acc, s_poly + j * SABER_N, poly);
        }
        else
        {
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

        if(i == 0)
        {
            polymul_dec(acc, buff, poly);
        }else{
            polymla_dec(acc, buff, poly);
        }

    }

    BS2POLT(ciphertext + SABER_POLYVECCOMPRESSEDBYTES, buff);

    for (i = 0; i < SABER_N; i++) {
        poly[i] = (acc[i] + h2 - (buff[i] << (SABER_EP - SABER_ET))) >> (SABER_EP - 1);
    }

    POLmsg2BS(m, poly);

}




