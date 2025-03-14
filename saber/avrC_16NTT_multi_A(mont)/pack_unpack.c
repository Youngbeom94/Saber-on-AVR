#include "pack_unpack.h"
#include <stddef.h>

////////////////////////////////////////////////////////////////////////////////////////////////
///                                 Polynomial to bit-string                                 ///
////////////////////////////////////////////////////////////////////////////////////////////////


/* This function reduces its input mod 2 */
void POLmsg2BS(uint8_t bytes[SABER_KEYBYTES], const uint16_t data[SABER_N])
{
    size_t i, j;
    uint8_t byte;
    for (j = 0; j < SABER_KEYBYTES; j++) {
        byte = 0;
        for (i = 0; i < 8; i++) {
            byte |= ((data[j * 8 + i] & 0x01) << i);
        }
        bytes[j] = byte;
    }
}

/* This function reduces its input mod p */
void POLp2BS(uint8_t bytes[SABER_POLYCOMPRESSEDBYTES], const uint16_t data[SABER_N])
{
    size_t i;
    const uint16_t *in = data;
    uint8_t *out = bytes;
    for (i = 0; i < SABER_N / 4; i++) {
        out[0] = (uint8_t) (in[0]);
        out[1] = (uint8_t) (((in[0] >> 8) & 0x03) | (in[1] << 2));
        out[2] = (uint8_t) (((in[1] >> 6) & 0x0f) | (in[2] << 4));
        out[3] = (uint8_t) (((in[2] >> 4) & 0x3f) | (in[3] << 6));
        out[4] = (uint8_t) (in[3] >> 2);
        in += 4;
        out += 5;
    }
}

/* This function reduces its input mod p */
uint32_t POLp2BS_cmp(const uint8_t bytes[SABER_POLYCOMPRESSEDBYTES], const uint16_t data[SABER_N])
{
    size_t j;
    const uint16_t *in = data;
    const uint8_t *out = bytes;
    uint32_t fail = 0;

    for (j = 0; j < SABER_N / 4; j++) {
        fail |= out[0] ^ (uint8_t) (in[0]);
        fail |= out[1] ^ (uint8_t) (((in[0] >> 8) & 0x03) | (in[1] << 2));
        fail |= out[2] ^ (uint8_t) (((in[1] >> 6) & 0x0f) | (in[2] << 4));
        fail |= out[3] ^ (uint8_t) (((in[2] >> 4) & 0x3f) | (in[3] << 6));
        fail |= out[4] ^ (uint8_t) (in[3] >> 2);
        in += 4;
        out += 5;
    }
    return fail;
}

/* This function reduces its input mod q */
void POLq2BS(uint8_t bytes[SABER_POLYBYTES], const uint16_t data[SABER_N])
{
    size_t i;
    const uint16_t *in = data;
    uint8_t *out = bytes;
    for (i = 0; i < SABER_N / 8; i++) {
        out[0] = (uint8_t) (in[0]);
        out[1] = (uint8_t) (((in[0] >> 8) & 0x1f) | (in[1] << 5));
        out[2] = (uint8_t) (in[1] >> 3);
        out[3] = (uint8_t) (((in[1] >> 11) & 0x03) | (in[2] << 2));
        out[4] = (uint8_t) (((in[2] >> 6) & 0x7f) | (in[3] << 7));
        out[5] = (uint8_t) (in[3] >> 1);
        out[6] = (uint8_t) (((in[3] >> 9) & 0x0f) | (in[4] << 4));
        out[7] = (uint8_t) (in[4] >> 4);
        out[8] = (uint8_t) (((in[4] >> 12) & 0x01) | (in[5] << 1));
        out[9] = (uint8_t) (((in[5] >> 7) & 0x3f) | (in[6] << 6));
        out[10] = (uint8_t) (in[6] >> 2);
        out[11] = (uint8_t) (((in[6] >> 10) & 0x07) | (in[7] << 3));
        out[12] = (uint8_t) (in[7] >> 5);
        in += 8;
        out += 13;
    }
}


/* This function reduces its input mod T */
void POLT2BS(uint8_t bytes[SABER_SCALEBYTES_KEM], const uint16_t data[SABER_N])
{
    size_t j;
    const uint16_t *in = data;
    uint8_t *out = bytes;
#if SABER_ET == 3 // LightSaber
    for (j = 0; j < SABER_N / 8; j++) {
        out[0] = (uint8_t) ((in[0] & 0x7) | ((in[1] & 0x7) << 3) | (in[2] << 6));
        out[1] = (uint8_t) (((in[2] >> 2) & 0x01) | ((in[3] & 0x7) << 1) | ((in[4] & 0x7) << 4) | (in[5] << 7));
        out[2] = (uint8_t) (((in[5] >> 1) & 0x03) | ((in[6] & 0x7) << 2) | (in[7] << 5));
        in += 8;
        out += 3;
    }
#elif SABER_ET == 4 // Saber
    for (j = 0; j < SABER_N / 2; j++) {
        out[0] = (uint8_t) ((in[0] & 0x0f) | (in[1] << 4));
        in += 2;
        out += 1;
    }
#elif SABER_ET == 6 // FireSaber
    for (j = 0; j < SABER_N / 4; j++) {
        out[0] = (uint8_t) ((in[0] & 0x3f) | (in[1] << 6));
        out[1] = (uint8_t) (((in[1] >> 2) & 0x0f) | (in[2] << 4));
        out[2] = (uint8_t) (((in[2] >> 4) & 0x03) | (in[3] << 2));
        in += 4;
        out += 3;
    }
#else
#error "Unsupported SABER parameter."
#endif
}

/* This function reduces its input mod T */
uint32_t POLT2BS_cmp(const uint8_t bytes[SABER_SCALEBYTES_KEM], const uint16_t data[SABER_N])
{
    size_t j;
    const uint16_t *in = data;
    const uint8_t *out = bytes;
    uint32_t fail = 0;
#if SABER_ET == 3 // LightSaber
    for (j = 0; j < SABER_N / 8; j++) {
        fail |= out[0] ^ (uint8_t) ((in[0] & 0x7) | ((in[1] & 0x7) << 3) | (in[2] << 6));
        fail |= out[1] ^ (uint8_t) (((in[2] >> 2) & 0x01) | ((in[3] & 0x7) << 1) | ((in[4] & 0x7) << 4) | (in[5] << 7));
        fail |= out[2] ^ (uint8_t) (((in[5] >> 1) & 0x03) | ((in[6] & 0x7) << 2) | (in[7] << 5));
        in += 8;
        out += 3;
    }
#elif SABER_ET == 4 // Saber
    for (j = 0; j < SABER_N / 2; j++) {
        fail |= out[0] ^ (uint8_t) ((in[0] & 0x0f) | (in[1] << 4));
        in += 2;
        out += 1;
    }
#elif SABER_ET == 6 // FireSaber
    for (j = 0; j < SABER_N / 4; j++) {
        fail |= out[0] ^ (uint8_t) ((in[0] & 0x3f) | (in[1] << 6));
        fail |= out[1] ^ (uint8_t) (((in[1] >> 2) & 0x0f) | (in[2] << 4));
        fail |= out[2] ^ (uint8_t) (((in[2] >> 4) & 0x03) | (in[3] << 2));
        in += 4;
        out += 3;
    }
#else
#error "Unsupported SABER parameter."
#endif
    return fail;
}

////////////////////////////////////////////////////////////////////////////////////////////////
///                                 Bit-string to polynomial                                 ///
////////////////////////////////////////////////////////////////////////////////////////////////

/* This function does NOT reduce its output mod 2 */
void BS2POLmsg(const uint8_t bytes[SABER_KEYBYTES], uint16_t data[SABER_N])
{
    size_t i, j;
    uint8_t byte;
    for (j = 0; j < SABER_KEYBYTES; j++) {
        byte = bytes[j];
        for (i = 0; i < 8; i++) {
            data[j * 8 + i] = byte >> i;
        }
    }
}

/* This function does NOT reduce its output mod p */
void BS2POLp(const uint8_t bytes[SABER_POLYCOMPRESSEDBYTES], uint16_t data[SABER_N])
{
    size_t j;
    const uint8_t *in = bytes;
    int16_t *out = (int16_t *)data;

    for (j = 0; j < SABER_N / 4; j++) {
        out[0] = (in[0]) | (in[1] << 8);
        out[1] = (in[1] >> 2) | (in[2] << 6);
        out[2] = (in[2] >> 4) | (in[3] << 4);
        out[3] = (in[3] >> 6) | (in[4] << 2);
        in += 5;
        out += 4;
    }
}

/* This function does NOT reduce its output mod q */
void BS2POLq(const uint8_t bytes[SABER_POLYBYTES], uint16_t data[SABER_N])
{
    size_t i;
    const uint8_t *in = bytes;
    int16_t *out = (int16_t *)data;

    for (i = 0; i < SABER_N / 8; i++) {
        out[0] = (in[0]) | (in[1] << 8);
        out[1] = (in[1] >> 5) | (in[2] << 3) | (in[3] << 11);
        out[2] = (in[3] >> 2) | (in[4] << 6);
        out[3] = (in[4] >> 7) | (in[5] << 1) | (in[6] << 9);
        out[4] = (in[6] >> 4) | (in[7] << 4) | (in[8] << 12);
        out[5] = (in[8] >> 1) | (in[9] << 7);
        out[6] = (in[9] >> 6) | (in[10] << 2) | (in[11] << 10);
        out[7] = (in[11] >> 3) | (in[12] << 5);
        in += 13;
        out += 8;
    }
}

/* This function does NOT reduce its output mod T */
void BS2POLT(const uint8_t bytes[SABER_SCALEBYTES_KEM], uint16_t data[SABER_N])
{
    size_t j;
    const uint8_t *in = bytes;
    uint16_t *out = data;
#if SABER_ET == 3 // LightSaber
    for (j = 0; j < SABER_N / 8; j++) {
        out[0] = in[0];
        out[1] = in[0] >> 3;
        out[2] = (in[0] >> 6) | (in[1] << 2);
        out[3] = in[1] >> 1;
        out[4] = in[1] >> 4;
        out[5] = (in[1] >> 7) | (in[2] << 1);
        out[6] = in[2] >> 2;
        out[7] = in[2] >> 5;
        in += 3;
        out += 8;
    }
#elif SABER_ET == 4 // Saber
    for (j = 0; j < SABER_N / 2; j++) {
        out[0] = in[0];
        out[1] = in[0] >> 4;
        in += 1;
        out += 2;
    }
#elif SABER_ET == 6 // FireSaber
    for (j = 0; j < SABER_N / 4; j++) {
        out[0] = in[0];
        out[1] = (in[0] >> 6) | (in[1] << 2);
        out[2] = (in[1] >> 4) | (in[2] << 4);
        out[3] = (in[2] >> 2);
        in += 3;
        out += 4;
    }
#else
#error "Unsupported SABER parameter."
#endif
}


////////////////////////////////////////////////////////////////////////////////////////////////
///                             Polynomial vector to bit-string                              ///
////////////////////////////////////////////////////////////////////////////////////////////////


void POLVECp2BS(uint8_t bytes[SABER_POLYVECCOMPRESSEDBYTES], const uint16_t data[SABER_L][SABER_N])
{
    size_t i;

    for (i = 0; i < SABER_L; i++) {
        /* This function reduces its input mod p */
        POLp2BS(&bytes[i * SABER_POLYCOMPRESSEDBYTES], data[i]);
    }
}

uint32_t POLVECp2BS_cmp(const uint8_t bytes[SABER_POLYVECCOMPRESSEDBYTES], const uint16_t data[SABER_L][SABER_N])
{
    size_t i;
    uint32_t fail = 0;

    for (i = 0; i < SABER_L; i++) {
        /* This function reduces its input mod p */
        fail |= POLp2BS_cmp(&bytes[i * SABER_POLYCOMPRESSEDBYTES], data[i]);
    }

    return fail;
}


void POLVECq2BS(uint8_t bytes[SABER_POLYVECBYTES], const uint16_t data[SABER_L][SABER_N])
{
    size_t i;

    for (i = 0; i < SABER_L; i++) {
        /* This function reduces its input mod q */
        POLq2BS(&bytes[i * SABER_POLYBYTES], data[i]);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////
///                             Bit-string to polynomial vector                              ///
////////////////////////////////////////////////////////////////////////////////////////////////

void BS2POLVECp(const uint8_t bytes[SABER_POLYVECCOMPRESSEDBYTES], uint16_t data[SABER_L][SABER_N])
{
    size_t i;
    for (i = 0; i < SABER_L; i++) {
        BS2POLp(bytes + i * (SABER_EP * SABER_N / 8), data[i]);
    }
}

void BS2POLVECq(const uint8_t bytes[SABER_POLYVECBYTES], uint16_t data[SABER_L][SABER_N])
{
    size_t i;
    for (i = 0; i < SABER_L; i++) {
        BS2POLq(&bytes[i * SABER_POLYBYTES], data[i]);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// The following functions are for compressed secret. Secrets are stored with their 4-bit value in [-SABER_MU/2, SABER_MU/2]. ///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/* This function reduces its input mod 2**4 */
void POLmu2BS(uint8_t bytes[SABER_N / 2], const uint16_t data[SABER_N])
{
    size_t j;
    const uint16_t *in = data;
    uint8_t *out = bytes;
    for (j = 0; j < SABER_N / 2; j++) {
        out[0] = (uint8_t) ((in[0] & 0x0f) | (in[1] << 4));
        in += 2;
        out += 1;
    }
}

/* This function does NOT reduce its output mod q */
void BS2POLmu(const uint8_t bytes[SABER_N / 2], uint16_t data[SABER_N])
{
    size_t j;
    const uint8_t *in = bytes;
    int16_t *out = (int16_t *)&data[0];

    struct int4_t { // bitfield struct to sign-extend 4-bit to 16-bit.
        signed int bits: 4;
    } s0, s1;

    for (j = 0; j < SABER_N / 2; j++) {
        s0.bits = (in[0] & 0xf);
        s1.bits = (in[0] >> 4) & 0xf;
        out[0] = (int16_t) s0.bits;
        out[1] = (int16_t) s1.bits;
        in += 1;
        out += 2;
    }
}

void POLVECmu2BS(uint8_t bytes[SABER_L * SABER_N / 2], const uint16_t data[SABER_L][SABER_N])
{
    size_t i;
    for (i = 0; i < SABER_L; i++) {
        POLmu2BS(bytes + i * SABER_POLYSECRETBYTES, data[i]);
    }
}

void BS2POLVECmu(const uint8_t bytes[SABER_L * SABER_N / 2], uint16_t data[SABER_L][SABER_N])
{
    size_t i;
    for (i = 0; i < SABER_L; i++) {
        BS2POLmu(bytes + i * SABER_POLYSECRETBYTES, data[i]);
    }
}



