#include "poly.h"
#include "cbd.h"
#include "fips202.h"
#include <string.h>
#include "pack_unpack.h"

#define h1 (1 << (SABER_EQ - SABER_EP - 1))
#define h2 ((1 << (SABER_EP - 2)) - (1 << (SABER_EP - SABER_ET - 1)) + (1 << (SABER_EQ - SABER_EP - 1)))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

// todo params of SABER for NTT multi moduli
#define NTT_Q1 7681
#define NTT_Q2 10753
#define CRT_U_PINV 32747 // CRT_U * P1^-1(-10751) mod 2^16
#define CRT_U 3563       //
#define MULTI_MONT_PINV (-9)
#define MULTI_MONT (-3593)
#define q1_MONT (4088)     // 2^16 mod q, (4088)
#define NTT_Q1_INV (-7679) // q^(-1) mod 2^16, ( 57857 )

#define q2_MONT 1018        // 2^16 mod q,
#define NTT_Q2_INV (-10751) // q^(-1) mod 2^16, ( 54785 )

__flash
static int16_t q1_streamlined_CT_table[SABER_N] = {0, 3777, -3182, 3625, -3696, -1100, 2456, 2194, 121, -2250, 834, -2495, -2319, 2876, -1701, 1414, 2816, -2088, -2237, 1986, -1599, 1993, 3706, -2006, -1525, -2557, 1296, 1483, -2830, 3364, 617, 1921, -3689, -1738, 3266, -3600, 810, 1887, -638, -7, -438, -679, -1305, -1760, 396, -3174, -3555, -1881, 3772, -2535, -2440, -2555, 1535, -549, 3153, 2310, -1399, 1321, 514, -2956, -103, 2804, -2043, -1431, -1054, 1698, -3456, 1166, 2426, 3831, 915, -2, -3417, -194, 2919, 2789, 3405, 2385, -2113, -2732, 2175, 373, 3692, -730, -1756, 3135, -2391, 660, -1497, 2572, -3145, 1350, -2224, -3588, -1681, 2883, -1390, 1598, 3750, 2762, 2835, 2764, -2233, 3816, -1533, 1464, -727, 1521, 1386, -3428, -921, -2743, -2160, 2649, -859, 2579, 1532, 1919, -486, 404, -1056, 783, 1799, -2665, 3480, 2133, -3310, -1168, -17, 3744, 2422, 2001, 1278, 929, -1348, -2230, -179, -1242, -2059, -1070, 2161, 1649, 2072, 3177, -2071, 1121, -436, 236, 715, 670, -658, -1476, -2378, 2767, 3542, -226, 1203, 1181, -151, -3794, 1712, -222, 2786, -451, -3547, 1779, -1151, -434, 3568, -3693, 3581, -1586, 1509, 2918, 2339, -1407, 3434, -3550, 2340, 2891, 2998, -3314, 3461, -2719, -2247, -2589, 1144, 1072, 1295, -2815, -3770, 3450, 3781, -2258, 796, 3163, -3208, -589, 2963, -124, 3214, 3334, -3366, -3745, 3723, 1931, -429, -402, -3408, 83, -1526, 826, -1338, 2345, -2303, 2515, -642, -1837, -2965, -791, 370, 293, 3312, 2083, -1689, -777, 2070, 2262, -893, 2386, -188, -1519, -2874, -1404, 1012, 2130, 1441, 2532, -3335, -1084, -3343, 2937, 509, -1403, 2812, 3763, 592, 2005, 3657, 2460, -3677, 3752, 692, 1669, 2167, -3287};
__flash
static int16_t q1_streamlined_CT_inv_table[SABER_N] = {-3593, -3593, -3777, -3593, -3625, -3777, 3182, -3593, -2194, -3625, 1100, -3777, -2456, 3182, 3696, -3593, -1414, -2194, 2495, -3625, -2876, 1100, 2250, -3777, 1701, -2456, -834, 3182, 2319, 3696, -121, -3593, -1921, -1414, 2006, -2194, -1483, 2495, -1986, -3625, -3364, -2876, -1993, 1100, 2557, 2250, 2088, -3777, -617, 1701, -3706, -2456, -1296, -834, 2237, 3182, 2830, 2319, 1599, 3696, 1525, -121, -2816, -3593, 1431, -1921, 1881, -1414, -2310, 2006, 7, -2194, 2956, -1483, 1760, 2495, 2555, -1986, 3600, -3625, -2804, -3364, 3174, -2876, 549, -1993, -1887, 1100, -1321, 2557, 679, 2250, 2535, 2088, 1738, -3777, 2043, -617, 3555, 1701, -3153, -3706, 638, -2456, -514, -1296, 1305, -834, 2440, 2237, -3266, 3182, 103, 2830, -396, 2319, -1535, 1599, -810, 3696, 1399, 1525, 438, -121, -3772, -2816, 3689, -3593, 1168, 1431, -2883, -1921, 2743, 1881, 2732, -1414, -404, -2310, -660, 2006, -3816, 7, 2, -2194, 2665, 2956, -1350, -1483, -1521, 1760, -2789, 2495, -2579, 2555, 730, -1986, -2762, 3600, -1166, -3625, -2133, -2804, 3588, -3364, 3428, 3174, -2385, -2876, -1919, 549, -3135, -1993, -2764, -1887, -3831, 1100, -783, -1321, -2572, 2557, -1464, 679, 194, 2250, -2649, 2535, -373, 2088, -1598, 1738, -1698, -3777, 3310, 2043, 1681, -617, 921, 3555, 2113, 1701, 486, -3153, 2391, -3706, 2233, 638, -915, -2456, -1799, -514, 3145, -1296, 727, 1305, -2919, -834, 859, 2440, -3692, 2237, -3750, -3266, 3456, 3182, -3480, 103, 2224, 2830, -1386, -396, -3405, 2319, -1532, -1535, 1756, 1599, -2835, -810, -2426, 3696, 1056, 1399, 1497, 1525, 1533, 438, 3417, -121, 2160, -3772, -2175, -2816, 1390, 3689, 1054, 0};
__flash
static int16_t q1_last_twist_mont_table[SABER_N] = {1912, -3438, -551, 115, -2352, -3259, -672, 3458, -192, 988, 3237, -815, -2367, 3059, 421, 874, -977, 1347, -3571, -2907, 77, 1364, 22, 1487, -1091, -2867, -1409, 3570, 1792, 1020, 512, 2486, -951, -387, -1369, 2084, -3683, 2790, 45, -3592, -3279, 71, 2355, -1077, -2619, -1405, 349, -2596, 1197, -1839, 342, -2720, 1195, 3612, 2536, 1032, -1470, -2997, -420, 241, -120, -3223, 1063, 2371, 1401, 2872, -697, -1374, -3491, 1802, -3192, -2777, -912, -2988, 1934, -1951, -1642, -2752, -3761, 311, 1120, -3203, 320, 3474, 2286, -1202, -3736, -2538, -3262, 3664, -932, -2245, 831, -2836, 2432, 287, -2597, 82, -742, 2218, -212, 1731, 2134, -1700, 1707, -1583, 1585, 645, -2839, -913, 3578, 3031, -75, 866, -2216, 2442, 3756, 1795, -3316, -2779, -3142, -794, -1995, 3065, -570, 1973, 3129, 1661, 894, -1720, 2450, -2686, 700, -2962, 200, 251, 3349, 1169, -2335, 334, 3722, 2290, 3258, -443, -2361, 2068, 1520, -2701, -663, -1869, -2384, -534, 3708, 2042, 3254, 2778, 2027, 1891, -3810, -557, 1106, -3451, 316, -986, -1007, -1379, -1385, -394, -1493, 2082, 1768, -2697, 3797, 1424, -2207, -2885, 1564, 273, -2845, 78, 2479, -1075, -389, -3599, -3403, 69, 125, 1117, 1133, 3611, 1421, 2129, 406, -489, 116, -1237, 3325, -2548, 950, -728, 2466, -208, -1490, -2254, -1523, -644, -3727, -184, 2227, 2142, -461, 612, -1229, -3117, -3643, 1304, 2251, -1822, -3746, 1674, 27, -619, 1105, 3115, 1413, 890, 1501, -843, -2863, 3051, -818, 1969, -1331, -1632, 717, 631, -3087, -917, -882, -262, -252, 3217, -72, -3470, 2174, -3186, -3768, 187, 1118, 2248, 2514, -455, -379, -130, 989, -3329};
__flash
static int16_t q2_streamlined_CT_table[SABER_N] = {0, 223, 4188, -3688, 2413, -3686, 357, -376, 2695, -730, 4855, 2236, -425, 4544, 3364, -3784, 4875, -1520, -5063, -4035, 2503, 918, -3012, 4347, 1931, -1341, -3823, -341, -4095, -5175, -2629, -5213, -3091, 4129, -2935, 2790, 268, 1284, 4, 3550, 2982, 1287, 205, 4513, -2565, -2178, 4616, -193, -4102, 4742, -4876, -4744, -2984, -3062, -847, -4379, -2388, -1009, -3085, -1299, -2576, 4189, 1085, 544, 5023, 794, -567, -3198, 4734, -2998, 3441, -5341, 675, 2271, 1615, -2213, 512, 2774, 3057, -2045, 3615, -1458, -909, 5114, 2981, -4977, -116, 4580, -454, -5064, 4808, -1841, -886, -1356, -4828, -5156, 2737, 4286, -3169, -578, 5294, -636, 400, 151, -2884, -336, -1006, -326, 1572, -2740, -779, 2206, -1586, 1068, -3715, -1268, 2684, -5116, 1324, 2973, -2234, -4123, 3337, -864, 472, -467, 970, 635, -573, 2230, -1132, -4621, 2624, -4601, 3570, -3760, -5309, 3453, -5215, 854, -4250, 2428, 1381, 5172, -5015, -4447, 3135, 2662, 3524, -1573, 2139, 458, -2196, -2657, 4782, -3410, 2062, 2015, -4784, 1635, 1349, -1722, 2909, -4359, 2680, 2087, 40, 3241, -2439, 2117, 2050, 2118, -4144, -274, 3148, -1930, 1992, 4408, 5005, -4428, 2419, 1639, 2283, -778, -2374, 663, 1409, -2237, -4254, -1122, 97, -5313, -3535, -2813, 5083, 279, 4328, 2279, 2151, 355, -4003, 1204, -5356, -624, 5120, -4519, -1689, 1056, 3891, -3827, 1663, -2625, -2449, 3995, -1160, 2788, -4540, 3125, 5068, 3096, 1893, -2807, -5268, 2205, -4889, -152, 569, 4973, -825, 4393, 4000, 1510, 3419, -3360, 693, -3260, 4967, 4859, 2963, 554, -5107, -73, -4891, -1927, 5334, 2605, 2487, -2529, -834, 1782, 1111, 2113, 4720, -4670, -1053, -4403};
__flash
static int16_t q2_streamlined_CT_inv_table[SABER_N] = {1018, 1018, -223, 1018, 3688, -223, -4188, 1018, 376, 3688, 3686, -223, -357, -4188, -2413, 1018, 3784, 376, -2236, 3688, -4544, 3686, 730, -223, -3364, -357, -4855, -4188, 425, -2413, -2695, 1018, 5213, 3784, -4347, 376, 341, -2236, 4035, 3688, 5175, -4544, -918, 3686, 1341, 730, 1520, -223, 2629, -3364, 3012, -357, 3823, -4855, 5063, -4188, 4095, 425, -2503, -2413, -1931, -2695, -4875, 1018, -544, 5213, 193, 3784, 4379, -4347, -3550, 376, 1299, 341, -4513, -2236, 4744, 4035, -2790, 3688, -4189, 5175, 2178, -4544, 3062, -918, -1284, 3686, 1009, 1341, -1287, 730, -4742, 1520, -4129, -223, -1085, 2629, -4616, -3364, 847, 3012, -4, -357, 3085, 3823, -205, -4855, 4876, 5063, 2935, -4188, 2576, 4095, 2565, 425, 2984, -2503, -268, -2413, 2388, -1931, -2982, -2695, 4102, -4875, 3091, 1018, -635, -544, 5156, 5213, -2206, 193, 2045, 3784, -2973, 4379, -4580, -4347, -151, -3550, 5341, 376, 864, 1299, 1841, 341, 326, -4513, 2213, -2236, 1268, 4744, -5114, 4035, 578, -2790, 3198, 3688, 467, -4189, 1356, 5175, 2740, 2178, -2774, -4544, 5116, 3062, 4977, -918, 636, -1284, 2998, 3686, 4123, 1009, 5064, 1341, 336, -1287, -2271, 730, -1068, -4742, 1458, 1520, -4286, -4129, -794, -223, -970, -1085, 4828, 2629, 779, -4616, -3057, -3364, -1324, 847, 116, 3012, -400, -4, -3441, -357, -3337, 3085, -4808, 3823, 1006, -205, -1615, -4855, 3715, 4876, 909, 5063, 3169, 2935, 567, -4188, -472, 2576, 886, 4095, -1572, 2565, -512, 425, -2684, 2984, -2981, -2503, -5294, -268, -4734, -2413, 2234, 2388, 454, -1931, 2884, -2982, -675, -2695, 1586, 4102, -3615, -4875, -2737, 3091, -5023, 0};
__flash
static int16_t q2_last_twist_mont_table[SABER_N] = {2536, -1897, -1265, 5250, 525, -5324, -2683, 807, 1156, -2035, 5173, -558, -4357, -1511, -3377, -1413, 934, 2244, 2375, -5139, 2712, -4030, -403, 1035, -5273, 548, 4356, -1715, 5205, -4856, 1665, -5210, -521, -3278, -4629, 2763, -799, 3146, -1836, 1967, 1272, -4174, -2568, -4558, -4757, -1551, -3381, -3564, -2507, -1326, 2018, 4503, -625, 5314, 2682, -4033, 672, -4234, -2574, -2408, -4542, 3847, 1460, 146, -2136, 1937, 1269, -3099, 2916, -1859, 3040, 304, 2181, 3444, 2495, -5127, -1588, -4460, -446, 2106, -1940, -194, -2170, -217, -1097, -1185, 5258, 4827, 1558, 4457, 1521, 3378, 4639, -2762, 4025, -4974, -2648, -4566, 1694, 2320, 232, -4278, -4729, 2753, -800, -80, -8, -4302, 3871, 3613, -714, -2222, 4079, -2818, -4583, 617, 1137, 1189, -3107, -1386, 2012, -4100, -410, -41, -3230, -323, 1043, -971, -3323, 743, -1001, -3326, 1818, 4483, -627, -1138, -4415, 4935, -4883, 587, 1134, 2264, 2377, 1313, -944, -2245, 5152, -3786, 1772, -4124, -2563, 819, -3144, -2465, 5130, 513, -1024, -2253, 850, 85, -5368, -4838, -4785, 4898, 4791, 3705, -5006, 1650, 165, -5360, -536, 2097, 1285, -5248, -4826, 1668, 4468, 4748, 4776, -1673, 908, 4392, -3862, 3915, -4985, 4878, 4789, -2747, -1350, -135, 5363, -539, 3172, -3984, -2549, 2971, 3523, -723, 1003, -975, 5279, -2698, -4571, -3683, 707, 1146, -2036, 1947, 1270, 127, 1088, 4410, 441, 3270, 327, 1108, 4412, -3860, -386, 2112, -4090, -409, 3185, -5058, -4807, -1556, 1995, -5177, -1593, 916, -2059, 3020, 302, -4271, -3653, 710, 71, 3233, -752, 4226, -1728, -4474, -2598, -4561, -3682, 3933, -682, 4233, -652, 4236, -1727, -1248, -4426, 1708, 4472, -3854};


int16_t q1_montgomery_reduce(int32_t a)
{
    int16_t t;

    t = (int16_t)a * NTT_Q1_INV;
    t = (a - (int32_t)t * NTT_Q1) >> 16;
    return t;
}

int16_t q2_montgomery_reduce(int32_t a)
{
    int16_t t;

    t = (int16_t)a * NTT_Q2_INV;
    t = (a - (int32_t)t * NTT_Q2) >> 16;
    return t;
}

void q1_ntt(int16_t a[SABER_N])
{
    unsigned int len, start, j, k;
    int16_t zeta, t;

    k = 0;
    for (len = 128; len > 0; len >>= 1)
    { // 128, 64, 32, 16, 8, 4, 2, 1
        for (start = 0; start < SABER_N; start = j + len)
        {
            zeta = q1_streamlined_CT_table[++k];
            for (j = start; j < start + len; ++j)
            {
                t = q1_montgomery_reduce((int32_t)zeta * a[j + len]);
                if (len == 8)
                {
                    a[j + len] = (a[j] - t) % NTT_Q1;
                    a[j] = (a[j] + t) % NTT_Q1;
                }
                else
                {
                    a[j + len] = (a[j] - t);
                    a[j] = (a[j] + t);
                }
            }
        }
    }
}

void q2_ntt(int16_t a[SABER_N])
{
    unsigned int len, start, j, k;
    int16_t zeta, t;

    k = 0;
    for (len = 128; len > 0; len >>= 1)
    {
        for (start = 0; start < SABER_N; start = j + len)
        {
            zeta = q2_streamlined_CT_table[++k];
            for (j = start; j < start + len; ++j)
            {
                t = q2_montgomery_reduce((int32_t)zeta * a[j + len]);
                if (len == 32 || len == 2)
                {
                    // if (len == 8) {
                    a[j + len] = (a[j] - t) % NTT_Q2;
                    a[j] = (a[j] + t) % NTT_Q2;
                }
                else
                {
                    a[j + len] = (a[j] - t);
                    a[j] = (a[j] + t);
                }
            }
        }
    }
}

void q1_ct_invntt(int16_t a[SABER_N])
{
    unsigned int len, i, j, k;
    int16_t t, zeta;

    k = 0;
    for (len = 1; len < SABER_N; len <<= 1)
    {
        for (i = 0; i < len; i++)
        {
            zeta = q1_streamlined_CT_inv_table[k++];

            for (j = 0; j < SABER_N / (len << 1); j++)
            {
                t = q1_montgomery_reduce((int32_t)a[j * (len << 1) + i + len] * zeta);            
                if (len == 32)
                {
                    a[j * (len << 1) + i + len] = (a[j * (len << 1) + i] - t) % NTT_Q1;
                    a[j * (len << 1) + i] = (a[j * (len << 1) + i] + t) % NTT_Q1;
                }
                else
                {
                    a[j * (len << 1) + i + len] = (a[j * (len << 1) + i] - t);
                    a[j * (len << 1) + i] = (a[j * (len << 1) + i] + t);
                }
            }
        }
    }

    for (j = 0; j < SABER_N; ++j)
    {
        a[j] = q1_montgomery_reduce((int32_t)a[j] * q1_last_twist_mont_table[j]);        
    }
}

void q2_ct_invntt(int16_t a[SABER_N])
{
    unsigned int len, i, j, k;
    int16_t t, zeta;

    k = 0;
    for (len = 1; len < SABER_N; len <<= 1)
    {
        for (i = 0; i < len; i++)
        {
            zeta = q2_streamlined_CT_inv_table[k++];

            for (j = 0; j < SABER_N / (len << 1); j++)
            {
                t = q2_montgomery_reduce((int32_t)a[j * (len << 1) + i + len] * zeta);
                if (len == 4 || len == 64)
                {
                    a[j * (len << 1) + i + len] = (a[j * (len << 1) + i] - t) % NTT_Q2;
                    a[j * (len << 1) + i] = (a[j * (len << 1) + i] + t) % NTT_Q2;
                }
                else
                {
                    a[j * (len << 1) + i + len] = (a[j * (len << 1) + i] - t);
                    a[j * (len << 1) + i] = (a[j * (len << 1) + i] + t);
                }
            }
        }
    }

    for (j = 0; j < SABER_N; ++j)
    {
        a[j] = q2_montgomery_reduce((int32_t)a[j] * q2_last_twist_mont_table[j]);
    }
}

int16_t crt_mode(const int16_t a, int16_t b_pinv, int16_t b, int16_t p)
{
    int16_t t, u;
    t = ((int32_t)a * b_pinv) & 0xffff;
    u = (((int32_t)a * b) >> 16) & 0xffff;
    t = (((int32_t)t * p) >> 16) & 0xffff;
    t = u - t;
    return t;
};

void ntt_crt(int16_t *a, int16_t *b, uint16_t *c)
{
    int16_t tmp0[SABER_N] = {0x00};
    int16_t tmp1[SABER_N] = {0x00};

    for (int32_t i = 0; i < SABER_N; i++) tmp0[i] = crt_mode(a[i], MULTI_MONT_PINV, MULTI_MONT, NTT_Q1);

    for (int32_t i = 0; i < SABER_N; i++) tmp1[i] = b[i] - tmp0[i];

    for (int32_t i = 0; i < SABER_N; i++) tmp1[i] = crt_mode(tmp1[i], CRT_U_PINV, CRT_U, NTT_Q2);

    for (int32_t i = 0; i < SABER_N; i++) tmp1[i] = ((int32_t)tmp1[i] * NTT_Q1) & 0xffff;

    for (int32_t i = 0; i < SABER_N; i++) c[i] = (uint16_t)(tmp0[i] + tmp1[i]) & 8191;
}

void ntt_crt_acc(int16_t *a, int16_t *b, uint16_t *c)
{
    int16_t tmp0[SABER_N] = {0x00};
    int16_t tmp1[SABER_N] = {0x00};

    for (int32_t i = 0; i < SABER_N; i++) tmp0[i] = crt_mode(a[i], MULTI_MONT_PINV, MULTI_MONT, NTT_Q1);

    for (int32_t i = 0; i < SABER_N; i++) tmp1[i] = b[i] - tmp0[i];

    for (int32_t i = 0; i < SABER_N; i++) tmp1[i] = crt_mode(tmp1[i], CRT_U_PINV, CRT_U, NTT_Q2);

    for (int32_t i = 0; i < SABER_N; i++) tmp1[i] = ((int32_t)tmp1[i] * NTT_Q1) & 0xffff;

    for (int32_t i = 0; i < SABER_N; i++)
    {
        c[i] += (uint16_t)((uint16_t)(tmp0[i] + tmp1[i]) & 8191);
        c[i] &= 8191;
    } 
}

static inline void polymul(uint16_t des[SABER_N], uint16_t poly1[SABER_N], uint16_t poly2[SABER_N])
{
    int16_t a_tmp[SABER_N] = {0};
    int16_t b_tmp[SABER_N] = {0};
    int16_t q1_res[SABER_N] = {0};
    int16_t q2_res[SABER_N] = {0};
   
    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) 
    {
        poly2[cnt_i] = poly2[cnt_i] & 8191;        
        a_tmp[cnt_i] = (int16_t)poly1[cnt_i];
        b_tmp[cnt_i] = (int16_t)(poly2[cnt_i]);
        
        if (b_tmp[cnt_i] >= (8192 / 2)) b_tmp[cnt_i] -= 8192;
        if (b_tmp[cnt_i] < -(8192 / 2)) b_tmp[cnt_i] += 8192;
    }

    q1_ntt(a_tmp);
    q1_ntt(b_tmp);

    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) q1_res[cnt_i] = q1_montgomery_reduce((int32_t)a_tmp[cnt_i] * b_tmp[cnt_i]);
    
    q1_ct_invntt(q1_res);

   for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) 
    {
        poly2[cnt_i] = poly2[cnt_i] & 8191;        
        a_tmp[cnt_i] = (int16_t)poly1[cnt_i];
        b_tmp[cnt_i] = (int16_t)(poly2[cnt_i]);
        
        if (b_tmp[cnt_i] >= (8192 / 2)) b_tmp[cnt_i] -= 8192;
        if (b_tmp[cnt_i] < -(8192 / 2)) b_tmp[cnt_i] += 8192;
    }

    q2_ntt(a_tmp);
    q2_ntt(b_tmp);

    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) q2_res[cnt_i] = q2_montgomery_reduce((int32_t)a_tmp[cnt_i] * b_tmp[cnt_i]);

    q2_ct_invntt(q2_res);

    ntt_crt(q1_res, q2_res, des);
}

static inline void polymul_dec(uint16_t des[SABER_N], uint16_t poly1[SABER_N], uint16_t poly2[SABER_N])
{

    int16_t a_tmp[SABER_N] = {0};
    int16_t b_tmp[SABER_N] = {0};
    int16_t q1_res[SABER_N] = {0};
    int16_t q2_res[SABER_N] = {0};
   
    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) 
    {
        poly1[cnt_i] = poly1[cnt_i] & 8191;
        poly2[cnt_i] = poly2[cnt_i] & 8191;

        a_tmp[cnt_i] = (int16_t)poly1[cnt_i];
        if (a_tmp[cnt_i] >= 8192 / 2) a_tmp[cnt_i] -= 8192;
        if (a_tmp[cnt_i] < -8192 / 2) a_tmp[cnt_i] += 8192;

        b_tmp[cnt_i] = (int16_t)(poly2[cnt_i]);
        if (b_tmp[cnt_i] >= 8192 / 2) b_tmp[cnt_i] -= 8192;
        if (b_tmp[cnt_i] < -8192 / 2)  b_tmp[cnt_i] += 8192;
    }

    q1_ntt(a_tmp);
    q1_ntt(b_tmp);

    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) q1_res[cnt_i] = q1_montgomery_reduce((int32_t)a_tmp[cnt_i] * b_tmp[cnt_i]);
    
    q1_ct_invntt(q1_res);

   for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) 
    {
        poly1[cnt_i] = poly1[cnt_i] & 8191;
        poly2[cnt_i] = poly2[cnt_i] & 8191;

        a_tmp[cnt_i] = (int16_t)poly1[cnt_i];
        if (a_tmp[cnt_i] >= 8192 / 2) a_tmp[cnt_i] -= 8192;
        if (a_tmp[cnt_i] < -8192 / 2) a_tmp[cnt_i] += 8192;

        b_tmp[cnt_i] = (int16_t)(poly2[cnt_i]);
        if (b_tmp[cnt_i] >= 8192 / 2) b_tmp[cnt_i] -= 8192;
        if (b_tmp[cnt_i] < -8192 / 2)  b_tmp[cnt_i] += 8192;
    }

    q2_ntt(a_tmp);
    q2_ntt(b_tmp);

    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) q2_res[cnt_i] = q2_montgomery_reduce((int32_t)a_tmp[cnt_i] * b_tmp[cnt_i]);

    q2_ct_invntt(q2_res);

    ntt_crt(q1_res, q2_res, des);
}

static inline void polymla(uint16_t des[SABER_N], uint16_t poly1[SABER_N], uint16_t poly2[SABER_N])
{
    int16_t a_tmp[SABER_N] = {0};
    int16_t b_tmp[SABER_N] = {0};
    int16_t q1_res[SABER_N] = {0};
    int16_t q2_res[SABER_N] = {0};
   
    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) 
    {
        poly2[cnt_i] = poly2[cnt_i] & 8191;
        a_tmp[cnt_i] = (int16_t)(poly1[cnt_i]);
        b_tmp[cnt_i] = (int16_t)(poly2[cnt_i]);
        if (b_tmp[cnt_i] >= 8192 / 2) b_tmp[cnt_i] -= 8192;                        
        if (b_tmp[cnt_i] < -8192 / 2) b_tmp[cnt_i] += 8192;
    }

    q1_ntt(a_tmp);
    q1_ntt(b_tmp);
    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) q1_res[cnt_i] = q1_montgomery_reduce((int32_t)a_tmp[cnt_i] * b_tmp[cnt_i]);
    q1_ct_invntt(q1_res);

    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) 
    {
        poly2[cnt_i] = poly2[cnt_i] & 8191;
        a_tmp[cnt_i] = (int16_t)(poly1[cnt_i]);
        b_tmp[cnt_i] = (int16_t)(poly2[cnt_i]);
        if (b_tmp[cnt_i] >= 8192 / 2) b_tmp[cnt_i] -= 8192;                        
        if (b_tmp[cnt_i] < -8192 / 2) b_tmp[cnt_i] += 8192;
    }

    q2_ntt(a_tmp);
    q2_ntt(b_tmp);
    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) q2_res[cnt_i] = q2_montgomery_reduce((int32_t)a_tmp[cnt_i] * b_tmp[cnt_i]);
    q2_ct_invntt(q2_res);

    ntt_crt_acc(q1_res, q2_res, des);
}


static inline void polymla_dec(uint16_t des[SABER_N], uint16_t poly1[SABER_N], uint16_t poly2[SABER_N])
{
    int16_t a_tmp[SABER_N] = {0};
    int16_t b_tmp[SABER_N] = {0};
    int16_t q1_res[SABER_N] = {0};
    int16_t q2_res[SABER_N] = {0};
   
    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) 
    {
        poly1[cnt_i] = poly1[cnt_i] & 8191;
        if ((int16_t)poly1[cnt_i] >= 8192 / 2) poly1[cnt_i] -= 8192;
        if ((int16_t)poly1[cnt_i] < -8192 / 2) poly1[cnt_i] += 8192;
        poly2[cnt_i] = poly2[cnt_i] & 8191;
        if ((int16_t)poly2[cnt_i] >= 8192 / 2) poly2[cnt_i] -= 8192;
        if ((int16_t)poly2[cnt_i] < -8192 / 2) poly2[cnt_i] += 8192;

        a_tmp[cnt_i] = (int16_t)(poly1[cnt_i]);
        b_tmp[cnt_i] = (int16_t)(poly2[cnt_i]);
    }

    q1_ntt(a_tmp);
    q1_ntt(b_tmp);
    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) q1_res[cnt_i] = q1_montgomery_reduce((int32_t)a_tmp[cnt_i] * b_tmp[cnt_i]);
    q1_ct_invntt(q1_res);

    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) 
    {
        a_tmp[cnt_i] = (int16_t)(poly1[cnt_i]);
        b_tmp[cnt_i] = (int16_t)(poly2[cnt_i]);
    }

    q2_ntt(a_tmp);
    q2_ntt(b_tmp);
    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) q2_res[cnt_i] = q2_montgomery_reduce((int32_t)a_tmp[cnt_i] * b_tmp[cnt_i]);
    q2_ct_invntt(q2_res);

    ntt_crt_acc(q1_res, q2_res, des);
}

static inline shake128incctx shake128_absorb_seed(const uint8_t seed[SABER_SEEDBYTES])
{
    shake128incctx ctx;
    shake128_inc_init(&ctx);
    shake128_inc_absorb(&ctx, seed, SABER_SEEDBYTES);
    shake128_inc_finalize(&ctx);

    return ctx;
}

void MatrixVectorMulKeyPairNTT_D(uint8_t pk[SABER_INDCPA_PUBLICKEYBYTES], uint8_t sk[SABER_INDCPA_SECRETKEYBYTES])
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

    for (i = 0; i < SABER_L; i++)
    {

        shake128_inc_squeeze(shake_out, SABER_POLYCOINBYTES, &shake_s_ctx);
        cbd(s_poly, shake_out);
#ifdef SABER_COMPRESS_SECRETKEY
        POLmu2BS(sk + i * SABER_POLYSECRETBYTES, s_poly); // sk <- s
#else
        POLq2BS(sk + i * SABER_POLYSECRETBYTES, s_poly);
#endif

        for (j = 0; j < SABER_L; j++)
        {

            shake128_inc_squeeze(shake_out, SABER_POLYBYTES, &shake_A_ctx);
            BS2POLq(shake_out, poly);

            if (i == 0)
            {
                polymul(acc + j * SABER_N, s_poly, poly);
            }
            else
            {
                polymla(acc + j * SABER_N, s_poly, poly);
            }
        }
    }

    shake128_inc_ctx_release(&shake_A_ctx);
    shake128_inc_ctx_release(&shake_s_ctx);

    for (i = 0; i < SABER_L; i++)
    {

        for (j = 0; j < SABER_N; j++)
        {
            poly[j] = ((acc[i * SABER_N + j] + h1) >> (SABER_EQ - SABER_EP));
        }

        POLp2BS(pk + i * SABER_POLYCOMPRESSEDBYTES, poly);
    }
}

uint32_t MatrixVectorMulEncNTT_D(uint8_t ct0[SABER_POLYVECCOMPRESSEDBYTES], uint8_t ct1[SABER_SCALEBYTES_KEM], const uint8_t seed_s[SABER_NOISE_SEEDBYTES], const uint8_t seed_A[SABER_SEEDBYTES], const uint8_t pk[SABER_INDCPA_PUBLICKEYBYTES], const uint8_t m[SABER_KEYBYTES], int compare)
{

    uint16_t poly[SABER_N];
    uint16_t s_poly[SABER_L * SABER_N];
    uint16_t acc[SABER_N];

    uint8_t shake_out[MAX(SABER_POLYBYTES, SABER_POLYCOINBYTES)];

    uint16_t *mp = poly;

    size_t i, j;
    uint32_t fail = 0;

    shake128incctx shake_s_ctx = shake128_absorb_seed(seed_s);

    for (i = 0; i < SABER_L; i++)
    {
        shake128_inc_squeeze(shake_out, SABER_POLYCOINBYTES, &shake_s_ctx);
        cbd(s_poly + i * SABER_N, shake_out);
    }

    shake128_inc_ctx_release(&shake_s_ctx);

    shake128incctx shake_A_ctx = shake128_absorb_seed(seed_A);

    for (i = 0; i < SABER_L; i++)
    {

        for (j = 0; j < SABER_L; j++)
        {

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

        for (j = 0; j < SABER_N; j++)
        {
            acc[j] = ((acc[j] + h1) >> (SABER_EQ - SABER_EP));
        }

        if (compare)
        {
            fail |= POLp2BS_cmp(ct0 + i * SABER_POLYCOMPRESSEDBYTES, acc);
        }
        else
        {
            POLp2BS(ct0 + i * SABER_POLYCOMPRESSEDBYTES, acc);
        }
    }

    shake128_inc_ctx_release(&shake_A_ctx);

    for (j = 0; j < SABER_L; j++)
    {

        BS2POLp(pk + j * SABER_POLYCOMPRESSEDBYTES, poly);

        if (j == 0)
        {
            polymul(acc, s_poly + j * SABER_N, poly);
        }
        else
        {
            polymla(acc, s_poly + j * SABER_N, poly);
        }
    }

    BS2POLmsg(m, mp);

    for (j = 0; j < SABER_N; j++)
    {
        acc[j] = (acc[j] - (mp[j] << (SABER_EP - 1)) + h1) >> (SABER_EP - SABER_ET);
    }

    if (compare)
    {
        fail |= POLT2BS_cmp(ct1, acc);
    }
    else
    {
        POLT2BS(ct1, acc);
    }

    return fail;
}

void InnerProdDecNTT(uint8_t m[SABER_KEYBYTES], const uint8_t ciphertext[SABER_BYTES_CCA_DEC], const uint8_t sk[SABER_INDCPA_SECRETKEYBYTES])
{

    uint16_t poly[SABER_N];
    uint16_t buff[SABER_N];
    uint16_t acc[SABER_N];

    size_t i;

    for (i = 0; i < SABER_L; i++)
    {

#ifdef SABER_COMPRESS_SECRETKEY
        BS2POLmu(sk + i * SABER_POLYSECRETBYTES, buff);
#else
        BS2POLq(sk + i * SABER_POLYSECRETBYTES, buff);
#endif
        BS2POLp(ciphertext + i * SABER_POLYCOMPRESSEDBYTES, poly);

        if (i == 0)
        {
            polymul_dec(acc, buff, poly);
        }
        else
        {
            polymla_dec(acc, buff, poly);
        }
    }

    BS2POLT(ciphertext + SABER_POLYVECCOMPRESSEDBYTES, buff);

    for (i = 0; i < SABER_N; i++)
    {
        poly[i] = (acc[i] + h2 - (buff[i] << (SABER_EP - SABER_ET))) >> (SABER_EP - 1);
    }

    POLmsg2BS(m, poly);
}
