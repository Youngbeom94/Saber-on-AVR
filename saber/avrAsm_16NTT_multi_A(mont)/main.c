#include <stdio.h>
#include <stdint.h>
#include "api.h"


int main()
{
    int ret = 0;

    uint8_t pk[CRYPTO_PUBLICKEYBYTES] = {0};
    uint8_t sk[CRYPTO_SECRETKEYBYTES] = {0};
    uint8_t ss[CRYPTO_BYTES] = {0};
    uint8_t ct[CRYPTO_CIPHERTEXTBYTES] = {0};
    uint8_t re_ss[CRYPTO_BYTES] = {0};

    
    ret = crypto_kem_keypair(pk, sk);
    // for(int cnt_i = 0 ; cnt_i < CRYPTO_PUBLICKEYBYTES ; cnt_i ++) printf("%02X", pk[cnt_i]);
    // printf("\n\n");
    // for(int cnt_i = 0 ; cnt_i < CRYPTO_SECRETKEYBYTES ; cnt_i ++) printf("%02X", sk[cnt_i]);
    // printf("\n\n");
    ret = crypto_kem_enc(ct, ss, pk);
    // for(int cnt_i = 0 ; cnt_i < CRYPTO_BYTES ; cnt_i ++) printf("%02X", ss[cnt_i]);
    // printf("\n\n");
    // for(int cnt_i = 0 ; cnt_i < CRYPTO_CIPHERTEXTBYTES ; cnt_i ++) printf("%02X", ct[cnt_i]);
    ret = crypto_kem_dec(re_ss, ct, sk);

    if(ret != 0) return 1;

    for(int cnt_i = 0 ; cnt_i < CRYPTO_BYTES ; cnt_i ++)
    {
        if(ss[cnt_i] != re_ss[cnt_i]) return 1;
    }
    
    // SABER2 : ss== re_ss : A74E899A32B7A1461F7620ADDA4C0533F12C2FE6A9212B81766AF34D4A9366C3
    // SABER3 : ss== re_ss : CC968158CA4B5E20027A66A7BA0B721FFD5BC98A890854FBADE95CF57ABE2E7A
    // SABER4 : ss== re_ss : A2B013043B74B60A50FB568DE4449580ACAD9E96418AE7F69A38CF9DF8ED442D
    printf("\nss == re_ss : "); 
    for(int cnt_i = 0 ; cnt_i < CRYPTO_BYTES ; cnt_i ++) printf("%02X", ss[cnt_i]);
    printf("\ntest finish\n");
    return 0;
}