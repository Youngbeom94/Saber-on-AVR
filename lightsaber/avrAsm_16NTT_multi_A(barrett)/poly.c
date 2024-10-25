#include "poly.h"
#include "cbd.h"
#include "fips202.h"
#include <string.h>
#include "pack_unpack.h"

#define h1 (1 << (SABER_EQ - SABER_EP - 1))
#define h2 ((1 << (SABER_EP - 2)) - (1 << (SABER_EP - SABER_ET - 1)) + (1 << (SABER_EQ - SABER_EP - 1)))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

// 7 layer
#define Q1 3329
#define Q1pr 3
// omegaQ1 = Q1pr^((Q1 - 1) / (NTT_N << 1)) mod Q1
#define omegaQ1 (-268)
// invomegaQ1 = omegaQ1^{-1} mod Q1
#define invomegaQ1 (-1031)
// RmodQ1 = 2^16 mod^{+-} Q1
#define RmodQ1 (-1044)
// Q1prime = -Q1^{-1} mod^{+-} 2^16
#define Q1prime 3327
// RmodQ1hi = (R^2 mod^+- Q1) * Q1prime mod^+- R
#define RmodQ1hi (-20553)
// invNQ1 = NTT_N^{-1} mod Q1
#define invNQ1 (-26)
#define Q1_HF (Q1/2)
#define Q1_HFNEG -(Q1/2)

// 8 Layer
#define Q2 7681
#define Q2pr 17
// omegaQ2 = Q2pr^((Q2 - 1) / (NTT_N << 1)) mod Q2
#define omegaQ2 (-535)
// invomegaQ2 = omegaQ2^{-1} mod Q2
#define invomegaQ2 (-201)
// RmodQ2 = 2^16 mod^{+-} Q2
#define RmodQ2 (-3593)
// Q2prime = -Q2^{-1} mod^{+-} 2^16
#define Q2prime 7679
// RmodQ2hi = (R^2 mod^+- Q2) * Q2prime mod^+- R
#define RmodQ2hi (-30657)
// invNQ2 = NTT_N^{-1} mod Q2
#define invNQ2 (-30)
#define Q2_HF (Q2/2)
#define Q2_HFNEG -(Q2/2)

// CRT
#define CRT_U_PINV (26970)
#define CRT_U (-2726)
#define MONT_PINV (-20)
#define MONT (-1044)

__flash
int16_t q1_streamlined_CT_table[(SABER_N / 2) << 1] = {-1600, -31498, 749, 14745, 40, 787, -630, -12402, -687, -13525, 848, 16694, 1432, 28191, -797, -15690, 193, 3799, -1062, -20907, 1410, 27758, -569, -11202, 1583, 31164, -69, -1358, 543, 10690, -1333, -26242, -1089, -21438, 283, 5571, -56, -1102, 882, 17363, 296, 5827, 1476, 29057, -1339, -26360, 450, 8859, -936, -18426, 821, 16163, 1355, 26675, -535, -10532, 447, 8800, -1235, -24313, -1426, -28073, -1414, -27837, -1320, -25986, -464, -9134, 33, 650, -1352, -26616, -650, -12796, -632, -12442, -816, -16064, -1573, -30967, 76, 1496, 289, 5689, 331, 6516, -1052, -20710, -1274, -25080, 1025, 20179, 1197, 23565, 648, 12757, -1481, -29156, -682, -13426, -712, -14017, 1227, 24155, 910, 17915, 219, 4311, -855, -16832, -461, -9075, -1438, -28309, 927, 18249, 1534, 30199, 807, 15887, 452, 8898, -1435, -28250, -1010, -19883, -268, -5276, -641, -12619, -992, -19529, -733, -14430, -939, -18486, 1021, 20100, -892, -17560, -941, -18525, 540, 10631, 1540, 30317, 1651, 32502, 1626, 32010, -642, -12639, -1461, -28762, -1482, -29175, 952, 18741, 1041, 20494, -1100, -21655, 723, 14233, -1637, -32227, -17, -335, 568, 11182, 583, 11477, -680, -13387, -756, -14883, 1173, 23092, -314, -6182, -279, -5493, 233, 4587, 48, 945, 1409, 27738, -667, -13131, -554, -10906, 886, 17442, 1179, 23210, 1143, 22502, -525, -10335, 1092, 21498, -403, -7934, -1026, -20198, -1219, -23998, -394, -7756, -885, -17422, 1175, 23132, -1029, -20257, -1455, -28644, 1607, 31636, -1212, -23860, -556, -10946, 757, 14903, -319, -6280, 1063, 20927, 735, 14469, -863, -16989, 1230, 24214, -561, -11044, 375, 7382, -780, -15355, 1239, 24391, -1645, -32384, 109, 2146, -1292, -25435, -1584, -31183, 1031, 20297};
__flash
int16_t q1_streamlined_CT_asm_inv_table[(SABER_N / 2) << 1] = {1600, 31498, -40, -787, 1600, 31498, -749, -14745, -1432, -28191, -40, -787, 687, 13525, 1600, 31498, -848, -16694, -749, -14745, 630, 12402, -543, -10690, -1432, -28191, -1410, -27758, -40, -787, -1583, -31164, 687, 13525, -193, -3799, 1600, 31498, 69, 1358, -848, -16694, 1062, 20907, -749, -14745, 569, 11202, 630, 12402, 797, 15690, 1426, 28073, -543, -10690, 1339, 26360, -1432, -28191, -1355, -26675, -1410, -27758, 56, 1102, -40, -787, -447, -8800, -1583, -31164, -296, -5827, 687, 13525, 936, 18426, -193, -3799, 1089, 21438, 1600, 31498, 1235, 24313, 69, 1358, -1476, -29057, -848, -16694, -821, -16163, 1062, 20907, -283, -5571, -749, -14745, 535, 10532, 569, 11202, -882, -17363, 630, 12402, -450, -8859, 797, 15690, 1333, 26242, 1010, 19883, 1426, 28073, -1197, -23565, -543, -10690, 855, 16832, 1339, 26360, 816, 16064, -1432, -28191, -1534, -30199, -1355, -26675, -331, -6516, -1410, -27758, 712, 14017, 56, 1102, -33, -650, -40, -787, -452, -8898, -447, -8800, 1274, 25080, -1583, -31164, -910, -17915, -296, -5827, 650, 12796, 687, 13525, 1438, 28309, 936, 18426, -76, -1496, -193, -3799, 1481, 29156, 1089, 21438, 1320, 25986, 1600, 31498, 1435, 28250, 1235, 24313, -1025, -20179, 69, 1358, -219, -4311, -1476, -29057, 632, 12442, -848, -16694, -927, -18249, -821, -16163, -289, -5689, 1062, 20907, 682, 13426, -283, -5571, 464, 9134, -749, -14745, -807, -15887, 535, 10532, 1052, 20710, 569, 11202, -1227, -24155, -882, -17363, 1352, 26616, 630, 12402, 461, 9075, -450, -8859, 1573, 30967, 797, 15690, -648, -12757, 1333, 26242, 1414, 27837};
__flash
int16_t q1_last_twist_table[(SABER_N / 2) << 1] = {-26, -512, 174, 3425, 372, 7323, -697, -13721, -457, -8997, -1551, -30534, 1161, 22856, 1449, 28526, 802, 15788, -1270, -25002, 1073, 21123, -1035, -20375, -1524, -30002, -44, -866, -1242, -24450, -1163, -22895, 613, 12068, 507, 9981, -64, -1260, -596, -11733, -1389, -27344, 589, 11595, -1381, -27187, -1001, -19706, 41, 807, 1006, 19805, 1462, 28782, 715, 14076, -1456, -28663, -243, -4784, 858, 16891, 916, 18033, 1040, 20474, -302, -5945, -1564, -30790, 1248, 24569, 1635, 32187, -1211, -23840, 166, 3268, -1367, -26911, 1210, 23821, 865, 17029, 357, 7028, 1452, 28585, 1038, 20434, -1569, -30888, -255, -5020, -86, -1693, -1217, -23958, -306, -6024, -769, -15139, 537, 10572, -1033, -20336, -257, -5059, -1353, -26636, 92, 1811, -1640, -32286, -292, -5748, 1442, 28388, 1361, 26793, 1647, 32423, -267, -5256, -1030, -20277, -21, -413, -1652, -32522, -1236, -24332, -691, -13603, 15, 295, 1180, 23230, -1495, -29431, 18, 354, 1416, 27876, 1535, 30219, -1310, -25789, -964, -18978, -1487, -29274, -1572, -30947, -491, -9666, 213, 4193, 111, 2185, -1255, -24706, -1076, -21183, 799, 15729, -1506, -29648, 1372, 27010, 293, 5768, 856, 16852, -351, -6910, -980, -19293, -1636, -32207, -1087, -21399, -1176, -23151, 700, 13780, 693, 13643, 1252, 24647, 840, 16537, -500, -9843, -495, -9745, 1008, 19844, -600, -11812, -594, -11694, -122, -2402, -720, -14174, -47, -925, -1478, -29096, -864, -17009, -1388, -27325, -442, -8701, -371, -7304, -334, -6575, 1467, 28880, -1111, -21872, 265, 5217, -237, -4666, 1330, 26183, 318, 6260, -1616, -31813, 1596, 31419, -950, -18702, 724, 14253, -748, -14725, -1140, -22442, 203, 3996, 434, 8544, -1368, -26931, -1088, -21419, -145, -2855, -310, -6103};

__flash
int16_t q2_streamlined_CT_table[(SABER_N) << 1] = {3383, 28865, 1925, 16425, -1213, -10350, -583, -4974, 1728, 14744, -849, -7244, 527, 4496, -2446, -20870, -2381, -20315, -97, -828, 2132, 18191, -2648, -22593, -2138, -18242, 2784, 23754, 1366, 11655, -2268, -19351, 675, 5759, -3092, -26382, 1286, 10972, 1112, 9488, -1794, -15307, -2399, -20469, 3000, 25597, 1846, 15750, 365, 3114, -2753, -23489, 3654, 31177, -878, -7491, 2273, 19394, -330, -2816, -2645, -22568, 3449, 29428, 528, 4505, 2941, 25093, 2508, 21399, 1655, 14121, -584, -4983, -1740, -14846, -2774, -23668, -2516, -21467, -1080, -9215, 3411, 29103, 2551, 21766, -243, -2073, -202, -1724, 766, 6536, 2881, 24581, -3074, -26228, 732, 6246, -3080, -26279, 3477, 29667, 2469, 21066, 3380, 28839, -1714, -14624, 693, 5913, -695, -5930, -799, -6817, -1381, -11783, -1875, -15998, -1908, -16279, -2724, -23242, -1382, -11792, 2423, 20674, 2028, 17303, 1591, 13575, 1952, 16655, -2044, -17440, 550, 4693, 1848, 15768, -1228, -10478, 1097, 9360, 1438, 12269, 2681, 22875, 2990, 25511, -707, -6032, -1125, -9599, -3780, -32252, 417, 3558, -2593, -22124, 1415, 12073, 1682, 14351, -2880, -24573, -3532, -30136, -3078, -26262, 2562, 21860, -3099, -26441, 648, 5529, 3041, 25946, 2844, 24266, 1003, 8558, -1853, -15810, 1408, 12013, 1044, 8908, -993, -8472, -2722, -23225, -2819, -24052, 3125, 26663, -3789, -32329, 1402, 11962, -257, -2193, -1478, -12611, -3141, -26800, -3180, -27132, -2264, -19317, -1155, -9855, -3073, -26220, -3566, -30426, -1220, -10409, -2563, -21868, 1886, 16092, -2573, -21953, 2900, 24743, 2063, 17602, -1587, -13541, 198, 1689, -880, -7508, 3188, 27201, 3501, 29871, -219, -1869, -3837, -32738, 319, 2722, 2897, 24718, -405, -3456, 1800, 15358, -1633, -13933, 869, 7415, -1996, -17030, -535, -4565, 2811, 23984, -621, -5299, 3751, 32004, -3016, -25733, -2760, -23549, 1036, 8839, 2252, 19215, 2840, 24232, -1211, -10333, -1872, -15972, -3832, -32695, 3376, 28805, -639, -5452, 674, 5751, -1115, -9513, -218, -1860, -118, -1007, 2805, 23933, 3280, 27986, -3483, -29718, -335, -2858, 738, 6297, 329, 2807, 3239, 27636, -3250, -27730, -1897, -16186, 3765, 32124, 1189, 10145, -2457, -20964, -113, -964, 1771, 15111, -1775, -15145, 1717, 14650, 1170, 9983, 2395, 20435, -2110, -18003, -2481, -21168, 1499, 12790, 1657, 14138, 1885, 16083, 1725, 14718, 3193, 27243, 2433, 20759, -572, -4880, 536, 4573, -2717, -23182, 2546, 21723, 856, 7304, 111, 947, -3615, -30844, -1393, -11885, 217, 1851, -3265, -27858, 2951, 25179, -2067, -17636, 3137, 26766, -2671, -22790, 1459, 12449, -3086, -26330, -793, -6766, -2050, -17491, 1994, 17013, 1784, 15221, -1959, -16715, 1406, 11996, 296, 2526, 2838, 24214, -2372, -20238, 2169, 18506, -3586, -30597, -3139, -26783, -1230, -10495, 2012, 17167, -2002, -17082, 1876, 16006, 2757, 23523, 2197, 18745, -346, -2952, -3006, -25648, 3394, 28958, -1193, -10179, -3081, -26288, 94, 802, 2996, 25563, -3452, -29453, -1131, -9650, -1035, -8831, 1437, 12261, -702, -5990, 1065, 9087, 506, 4317, -542, -4624, 2173, 18541, 1266, 10802, -3120, -26621, 2689, 22943, 2583, 22039, -669, -5708, 2668, 22764, -763, -6510, -413, -3524, -1704, -14539, 3799, 32414, -2358, -20119, 3445, 29394, 321, 2739, 2922, 24931, -185, -1578, -3694, -31518, -2799, -23882, 1656, 14129, 62, 529, 2359, 20128, -3546, -30255, 1604, 13686, 2259, 19274, -398, -3396, 1129, 9633, 1950, 16638, 1968, 16791, -1683, -14360, 1667, 14223, 1607, 13711, -2875, -24530, -1979, -16885, 3626, 30938, 201, 1715};
__flash
int16_t q2_streamlined_CT_asm_inv_table[(SABER_N) << 1] = {-3383, -28865, 1213, 10350, -3383, -28865, -1925, -16425, -527, -4496, 1213, 10350, -1728, -14744, -3383, -28865, 849, 7244, -1925, -16425, 583, 4974, -1366, -11655, -527, -4496, -2132, -18191, 1213, 10350, 2138, 18242, -1728, -14744, 2381, 20315, -3383, -28865, -2784, -23754, 849, 7244, 97, 828, -1925, -16425, 2648, 22593, 583, 4974, 2446, 20870, 2645, 22568, -1366, -11655, -3000, -25597, -527, -4496, -3654, -31177, -2132, -18191, -1286, -10972, 1213, 10350, -2273, -19394, 2138, 18242, 1794, 15307, -1728, -14744, -365, -3114, 2381, 20315, -675, -5759, -3383, -28865, 330, 2816, -2784, -23754, 2399, 20469, 849, 7244, 2753, 23489, 97, 828, 3092, 26382, -1925, -16425, 878, 7491, 2648, 22593, -1112, -9488, 583, 4974, -1846, -15750, 2446, 20870, 2268, 19351, -2423, -20674, 2645, 22568, -2881, -24581, -1366, -11655, -693, -5913, -3000, -25597, 2774, 23668, -527, -4496, 1875, 15998, -3654, -31177, -2551, -21766, -2132, -18191, -3477, -29667, -1286, -10972, -2508, -21399, 1213, 10350, 2724, 23242, -2273, -19394, 202, 1724, 2138, 18242, -3380, -28839, 1794, 15307, 584, 4983, -1728, -14744, 799, 6817, -365, -3114, 1080, 9215, 2381, 20315, -732, -6246, -675, -5759, -528, -4505, -3383, -28865, 1382, 11792, 330, 2816, -766, -6536, -2784, -23754, 1714, 14624, 2399, 20469, 1740, 14846, 849, 7244, 1381, 11783, 2753, 23489, -3411, -29103, 97, 828, 3080, 26279, 3092, 26382, -2941, -25093, -1925, -16425, 1908, 16279, 878, 7491, 243, 2073, 2648, 22593, -2469, -21066, -1112, -9488, -1655, -14121, 583, 4974, 695, 5930, -1846, -15750, 2516, 21467, 2446, 20870, 3074, 26228, 2268, 19351, -3449, -29428, 1996, 17030, -2423, -20674, 2722, 23225, 2645, 22568, 2573, 21953, -2881, -24581, 2593, 22124, -1366, -11655, 219, 1869, -693, -5913, -648, -5529, -3000, -25597, 3180, 27132, 2774, 23668, -1097, -9360, -527, -4496, 405, 3456, 1875, 15998, 1853, 15810, -3654, -31177, 3566, 30426, -2551, -21766, 707, 6032, -2132, -18191, -198, -1689, -3477, -29667, 3532, 30136, -1286, -10972, -1402, -11962, -2508, -21399, 2044, 17440, 1213, 10350, 1633, 13933, 2724, 23242, -1044, -8908, -2273, -19394, 2563, 21868, 202, 1724, 3780, 32252, 2138, 18242, -3188, -27201, -3380, -28839, -2562, -21860, 1794, 15307, 1478, 12611, 584, 4983, -1848, -15768, -1728, -14744, -319, -2722, 799, 6817, -2844, -24266, -365, -3114, 1155, 9855, 1080, 9215, -2681, -22875, 2381, 20315, -2063, -17602, -732, -6246, -1682, -14351, -675, -5759, -3125, -26663, -528, -4505, -1591, -13575, -3383, -28865, -869, -7415, 1382, 11792, 993, 8472, 330, 2816, -1886, -16092, -766, -6536, -417, -3558, -2784, -23754, -3501, -29871, 1714, 14624, 3099, 26441, 2399, 20469, 3141, 26800, 1740, 14846, 1228, 10478, 849, 7244, -2897, -24718, 1381, 11783, -1003, -8558, 2753, 23489, 3073, 26220, -3411, -29103, -2990, -25511, 97, 828, 1587, 13541, 3080, 26279, 2880, 24573, 3092, 26382, 3789, 32329, -2941, -25093, -1952, -16655, -1925, -16425, -1800, -15358, 1908, 16279, -1408, -12013, 878, 7491, 1220, 10409, 243, 2073, 1125, 9599, 2648, 22593, 880, 7508, -2469, -21066, 3078, 26262, -1112, -9488, 257, 2193, -1655, -14121, -550, -4693, 583, 4974, 3837, 32738, 695, 5930, -3041, -25946, -1846, -15750, 2264, 19317, 2516, 21467, -1438, -12269, 2446, 20870, -2900, -24743, 3074, 26228, -1415, -12073, 2268, 19351, 2819, 24052, -3449, -29428, -2028, -17303};
__flash
int16_t q2_last_twist_table[(SABER_N) << 1] = {-30, -256, -1651, -14087, 1568, 13379, -247, -2107, 3561, 30383, -1428, -12184, 2831, 24155, -637, -5435, -2540, -21672, 3594, 30665, -380, -3242, -430, -3669, 1939, 16544, 1992, 16996, -980, -8362, -2726, -23259, 2575, 21970, -2948, -25153, 1111, 9479, -562, -4795, -2253, -19223, -326, -2782, -3603, -30742, 2189, 18677, -2172, -18532, -1245, -10623, -3228, -27542, 3624, 30921, 1271, 10844, -1998, -17047, 2186, 18651, -1569, -13387, 448, 3822, 2124, 18122, 3212, 27405, -408, -3481, -2483, -21186, -182, -1553, -1823, -15554, -2265, -19325, 2086, 17798, 3169, 27039, 554, 4727, -3820, -32593, -280, -2389, 2513, 21441, 1833, 15640, 255, 2176, 2512, 21433, 2034, 17355, -1741, -14855, -3385, -28882, -3224, -27508, 2820, 24061, 1574, 13430, -1453, -12397, 175, 1493, 3230, 27559, 3655, 31185, 2721, 23216, -1570, -13396, 649, 5537, 128, 1092, -2685, -22909, 2015, 17192, 2078, 17730, -2904, -24778, -52, -444, 2771, 23643, 3742, 31928, 596, 5085, 3100, 26450, -939, -8012, -3286, -28037, -80, -683, 718, 6126, 1621, 13831, -3219, -27465, 1815, 15486, -3808, -32491, -2692, -22969, 3422, 29197, 3468, 29590, 1903, 16237, 1547, 13199, -3707, -31629, 50, 427, -2369, -20213, -53, -452, 2972, 25358, 1746, 14897, 2380, 20307, -2158, -18413, 3622, 30904, 1673, 14274, 1691, 14428, -1927, -16442, 3277, 27960, 1889, 16117, -3320, -28327, -927, -7909, 1983, 16919, 829, 7073, 2353, 20076, 3269, 27892, 3497, 29837, 3755, 32038, -2017, -17209, -1676, -14300, -1088, -9283, 3620, 30887, 2075, 17704, -2301, -19633, 1641, 14001, 442, 3771, 3330, 28412, -1083, -9240, 2615, 22312, -3307, -28216, -3540, -30204, -2793, -23830, 680, 5802, 1578, 13464, -2257, -19257, 478, 4078, 3775, 32209, 1644, 14027, -161, -1374, 1637, 13967, 1246, 10631, 3027, 25827, -1628, -13890, -3055, -26066, -425, -3626, 934, 7969, -3390, -28924, -2219, -18933, 521, 4445, 2813, 24001, 2981, 25435, -63, -538, -2699, -23028, -2852, -24334, -2823, -24086, -971, -8285, 3146, 26842, -2504, -21365, -3642, -31074, 2347, 20025, -3206, -27354, -798, -6809, -903, -7705, -2841, -24240, 2647, 22585, -2058, -17559, -1116, -9522, 1567, 13370, -46, -392, 1565, 13353, 356, 3037, -2427, -20708, -3757, -32056, 2419, 20639, -2316, -19761, -3025, -25810, 1226, 10461, -634, -5409, -3143, -26817, 1901, 16220, 1949, 16629, -18, -154, 3618, 30870, 2477, 21134, 1388, 11843, -2472, -21092, -2393, -20418, -2910, -24829, 1154, 9846, -1524, -13003, -916, -7816, -228, -1945, -258, -2201, -1909, -16288, -341, -2909, -588, -5017, 2973, 25366, 1545, 13182, -3305, -28199, 3739, 31902, 1199, 10230, -2888, -24641, -3268, -27883, -3698, -31552, -1759, -15008, 233, 1988, -747, -6374, -3473, -29632, -898, -7662, 3835, 32721, -2735, -23336, -3297, -28131, 2131, 18182, 1805, 15401, -1798, -15341, 391, 3336, -1781, -15196, -3026, -25819, 1427, 12175, -2630, -22440, -1359, -11595, -3357, -28643, -1171, -9991, -2740, -23378, -2292, -19556, -168, -1433, 3044, 25972, 2636, 22491, 153, 1305, -29, -247, -1852, -15802, 3564, 30409, -2031, -17329, 1138, 9710, 1692, 14437, -2128, -18157, -2408, -20546, 105, 896, 1938, 16535, 2193, 18711, -2976, -25392, -942, -8037, -2683, -22892, 1613, 13762, -1611, -13745, 1209, 10315, 2783, 23745, 1330, 11348, 1505, 12841, -2946, -25136, 709, 6049, 3430, 29266, 1860, 15870, 2509, 21407, 2637, 22499, -48, -410, 1967, 16783, -3636, -31023, 1141, 9735, 1089, 9292, -3821, -32602, -79, -674, 517, 4411, 3617, 30861, 2678, 22849, -608, -5188, -688, -5870};

extern void q1_barrett_ntt_asm(int16_t a[256]);
extern void q2_barrett_ntt_asm(int16_t a[256]);
extern void q1_barrett_invntt_asm(int16_t a[256]);
extern void q2_barrett_invntt_asm(int16_t a[256]);


static int16_t q1_barrett_reduce(int16_t a)
{
    int16_t t;
    int32_t rounding = 0x2000000;
    const int16_t v =  20159;

    t = ((int32_t)v * a + rounding) >> 26;
    t *= Q1;
    return a - t;
}

static int16_t q2_barrett_reduce(int16_t a)
{
    int16_t t;
    int32_t rounding = 0x2000000;
    const int16_t v = 8737;

    t = ((int32_t)v * a + rounding) >> 26;
    t *= Q2;
    return a - t;
}

static int ct_lt_s16(uint16_t x, uint16_t y) {return (x ^ ((x ^ (x - y)) & (y ^ (x - y)))) >> 15;}
static int ct_gt_s16(uint16_t x, uint16_t y) {return ct_lt_s16(y, x);}

static int16_t q1_freeze_16(int32_t a)
{    
    int16_t r;
    r = q1_barrett_reduce(a);

    uint16_t mask;
    mask = 0 -ct_gt_s16(r, Q1_HF); // if (a >= Q1 / 2) a -= Q1;
	r -= mask & Q1;
	mask = 0 - ct_lt_s16(r, (uint16_t) Q1_HFNEG); // if (a <= -Q1 / 2) a += Q1;
	r += mask & Q1;

    return r;
}

static int16_t q2_freeze_16(int32_t a)
{    
    int16_t r;
    r = q2_barrett_reduce(a);

    uint16_t mask;
    mask = 0 -ct_gt_s16(r, Q2_HF); // if (a >= Q2 / 2) a -= Q2;
	r -= mask & Q2;
	mask = 0 - ct_lt_s16(r, (uint16_t) Q2_HFNEG); // if (a <= -Q2 / 2) a += Q2;
	r += mask & Q2;
    
    return r;
}

// we assume -Q1 / 2 < a < Q1 / 2
static int16_t q1_gethi(int16_t a)
{

    int16_t lo, hi;

    int16_t bp_lo, bp_hi;
    int16_t a_lo, a_hi;

    // below we compute lo = a R mod^+- Q
    lo = a * RmodQ1;

    bp_lo = (int16_t)(uint8_t)RmodQ1hi;
    bp_hi = (int16_t)(int8_t)(RmodQ1hi >> 8);

    a_lo = (int16_t)(uint8_t)a;
    a_hi = (int16_t)(int8_t)(a >> 8);

    hi = a_hi * bp_hi;
    hi += ((a_lo * bp_hi) >> 8);
    hi += ((a_hi * bp_lo) >> 8);

    lo -= hi * Q1;

    lo = q1_freeze_16(lo);

    // return (a R mod^+- Q) * Qprime mod^+- R
    // we reduce mod^+- R by typing
    return lo * Q1prime;
}

// we assume -Q2 / 2 < a < Q2 / 2
static int16_t q2_gethi(int16_t a)
{

    int16_t lo, hi;

    int16_t bp_lo, bp_hi;
    int16_t a_lo, a_hi;

    // below we compute lo = a R mod^+- Q
    lo = a * RmodQ2;

    bp_lo = (int16_t)(uint8_t)RmodQ2hi;
    bp_hi = (int16_t)(int8_t)(RmodQ2hi >> 8);

    a_lo = (int16_t)(uint8_t)a;
    a_hi = (int16_t)(int8_t)(a >> 8);

    hi = a_hi * bp_hi;
    hi += ((a_lo * bp_hi) >> 8);
    hi += ((a_hi * bp_lo) >> 8);

    lo -= hi * Q2;

    lo = q2_freeze_16(lo);

    // return (a R mod^+- Q) * Qprime mod^+- R
    // we reduce mod^+- R by typing
    return lo * Q2prime;
}

static int16_t barrett_mul_approx(int16_t a, int16_t b, int16_t bp, int16_t q)
{
    int16_t lo, hi;
    int32_t full;

    lo = a * b;

    // we compute the full product and round it to the upper 16-bit instead.
    full = (int32_t)a * bp;
    full += 0x8000;
    hi = (int16_t)(full >> 16);

    return lo - hi * q;
}

void q1_barrett_ntt(int16_t a[SABER_N])
{
    q1_barrett_ntt_asm(a);
}

void q2_barrett_ntt(int16_t a[SABER_N])
{
    q2_barrett_ntt_asm(a);
}

void q1_basemul(int16_t r[2], const int16_t a[2], const int16_t b[2], int16_t zeta)
{
    int16_t tmp[2];
    tmp[0]  = barrett_mul_approx(a[1], b[1], q1_gethi(b[1]), Q1);
    tmp[0]  = barrett_mul_approx(tmp[0], zeta, q1_gethi(zeta), Q1);
    tmp[0] += barrett_mul_approx(a[0], b[0], q1_gethi(b[0]), Q1);
    tmp[1]  = barrett_mul_approx(a[0], b[1], q1_gethi(b[1]), Q1);
    tmp[1] += barrett_mul_approx(a[1], b[0], q1_gethi(b[0]), Q1);

    r[0] = tmp[0];
    r[1] = tmp[1];
}

void q1_point_mul_on_the_fly(int16_t res[SABER_N], int16_t src1[SABER_N], int16_t src2[SABER_N])
{
    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) src1[cnt_i] = q1_freeze_16(src1[cnt_i]);
    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++) src2[cnt_i] = q1_freeze_16(src2[cnt_i]);
        
    for (int i = 0; i < SABER_N / 4; i++)
    {
        q1_basemul(&res[4 * i],     &src1[4 * i],     &src2[4 * i],      q1_streamlined_CT_table[126 + (i << 1)]);
        q1_basemul(&res[4 * i + 2], &src1[4 * i + 2], &src2[4 * i + 2], -q1_streamlined_CT_table[126 + (i << 1)]);
    }
}

void q2_point_mul_on_the_fly(int16_t res[SABER_N], int16_t src1[SABER_N], int16_t src2[SABER_N])
{
    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++)
    {
        src2[cnt_i] = q2_freeze_16(src2[cnt_i]);
        res[cnt_i] = barrett_mul_approx(src1[cnt_i], src2[cnt_i], q2_gethi(src2[cnt_i]), Q2);
    }
}

void q1_point_mul_pre(int16_t res[SABER_N], int16_t src1[SABER_N])
{
    for (int cnt_i = 0; cnt_i < SABER_N / 2; cnt_i++)
    {
        res[2 * cnt_i + 0] = barrett_mul_approx(src1[2 * cnt_i + 0], q1_last_twist_table[2 * cnt_i + 0], q1_last_twist_table[2 * cnt_i + 1], Q1);
        res[2 * cnt_i + 1] = barrett_mul_approx(src1[2 * cnt_i + 1], q1_last_twist_table[2 * cnt_i + 0], q1_last_twist_table[2 * cnt_i + 1], Q1);
        res[2 * cnt_i + 0] = q1_freeze_16(res[2 * cnt_i + 0]);
        res[2 * cnt_i + 1] = q1_freeze_16(res[2 * cnt_i + 1]);
    }
}

void q2_point_mul_pre(int16_t res[SABER_N], int16_t src1[SABER_N])
{
    for (int cnt_i = 0; cnt_i < SABER_N; cnt_i++)
    {
        res[cnt_i] = barrett_mul_approx(src1[cnt_i], q2_last_twist_table[2 * cnt_i + 0], q2_last_twist_table[2 * cnt_i + 1], Q2);
        res[cnt_i] = q2_freeze_16(res[cnt_i]);
    }
}

void q1_barrett_invntt(int16_t a[SABER_N])
{
    q1_barrett_invntt_asm(a);
    q1_point_mul_pre(a, a);
}

void q2_barrett_invntt(int16_t a[SABER_N])
{
    q2_barrett_invntt_asm(a);
    q2_point_mul_pre(a, a);
}

int16_t crt_mode(const int16_t a, int16_t b_pinv, int16_t b, int16_t p){
  int16_t t,u;
  t = ((int32_t)a * b_pinv)& 0xffff;
  u = (((int32_t)a * b) >> 16) & 0xffff;
  t = (((int32_t)t * p) >> 16) & 0xffff;
  t = u-t;
  return t;
};
void ntt_crt(int16_t *a, int16_t *b, uint16_t *c) 
{
  int16_t tmp0[SABER_N] = {0x00};
  int16_t tmp1[SABER_N] = {0x00};
  
  for (int32_t i = 0; i < SABER_N; i++) {
    tmp0[i] = crt_mode(a[i],MONT_PINV,MONT,Q1 );
  }

  for (int32_t i = 0; i < SABER_N; i++) {
    tmp1[i] = b[i] - tmp0[i];
  }

  for (int32_t i = 0; i < SABER_N; i++) {
    tmp1[i] = crt_mode(tmp1[i],CRT_U_PINV,CRT_U,Q2);
  }

  for (int32_t i = 0; i < SABER_N; i++) {
    tmp1[i] = ((int32_t)tmp1[i] * Q1)& 0xffff;
  }

  for (int32_t i = 0; i < SABER_N; i++) {    
    c[i] = tmp0[i] + tmp1[i];
  }
}

static inline void double_NTT(int16_t NTT_q1[SABER_N], int16_t NTT_q2[SABER_N], uint16_t poly[SABER_N])
{

    for (size_t i = 0; i < SABER_N; i++)
    {
        NTT_q1[i] = (int16_t)(poly[i] & 8191);
        if (NTT_q1[i] >= 4096)
            NTT_q1[i] -= 8192;
        if (NTT_q1[i] < -4096)
            NTT_q1[i] += 8192;
        NTT_q2[i] = NTT_q1[i];
    }

    q1_barrett_ntt(NTT_q1);
    q2_barrett_ntt(NTT_q2);
}

static inline void double_basemul(int16_t des_q1[SABER_N], int16_t des_q2[SABER_N],
                                  int16_t src1_q1[SABER_N], int16_t src1_q2[SABER_N],
                                  int16_t src2_q1[SABER_N], int16_t src2_q2[SABER_N])
{
    q1_point_mul_on_the_fly(des_q1, src1_q1, src2_q1);
    q2_point_mul_on_the_fly(des_q2, src1_q2, src2_q2);
}

static inline void double_reduce(int16_t des_q1[SABER_N], int16_t des_q2[SABER_N],
                                 int16_t src_q1[SABER_N], int16_t src_q2[SABER_N])
{
    for (size_t i = 0; i < SABER_N; i++)
    {
        des_q1[i] = q1_barrett_reduce(src_q1[i]);
    }
    for (size_t i = 0; i < SABER_N; i++)
    {
        des_q2[i] = q2_barrett_reduce(src_q2[i]);
    }
}

static inline void double_add(int16_t des_q1[SABER_N], int16_t des_q2[SABER_N],
                              int16_t src1_q1[SABER_N], int16_t src1_q2[SABER_N],
                              int16_t src2_q1[SABER_N], int16_t src2_q2[SABER_N])
{
    for (size_t i = 0; i < SABER_N; i++)
    {
        des_q1[i] = src1_q1[i] + src2_q1[i];
    }
    for (size_t i = 0; i < SABER_N; i++)
    {
        des_q2[i] = src1_q2[i] + src2_q2[i];
    }
}

static inline void double_iNTT(int16_t src_q1[SABER_N], int16_t src_q2[SABER_N])
{
    q1_barrett_invntt(src_q1);
    q2_barrett_invntt(src_q2);
}

static inline shake128incctx shake128_absorb_seed(const uint8_t seed[SABER_SEEDBYTES])
{
    shake128incctx ctx;
    shake128_inc_init(&ctx);
    shake128_inc_absorb(&ctx, seed, SABER_SEEDBYTES);
    shake128_inc_finalize(&ctx);

    return ctx;
}

void MatrixVectorMulKeyPairNTT_A(uint8_t pk[SABER_INDCPA_PUBLICKEYBYTES], uint8_t sk[SABER_INDCPA_SECRETKEYBYTES])
{

    uint16_t poly[SABER_N];
    uint16_t acc[SABER_N];

    int16_t s_NTTq1[SABER_N];
    int16_t s_NTTq2[SABER_N];
    int16_t poly_NTTq1[SABER_N];
    int16_t poly_NTTq2[SABER_N];
    int16_t acc_q1[SABER_L * SABER_N];
    int16_t acc_q2[SABER_L * SABER_N];

    uint8_t shake_out[MAX(SABER_POLYBYTES, SABER_POLYCOINBYTES)];

    uint8_t *seed_A = pk + SABER_POLYVECCOMPRESSEDBYTES;
    uint8_t *seed_s = sk;

    size_t i, j;

    shake128incctx shake_s_ctx = shake128_absorb_seed(seed_s);
    shake128incctx shake_A_ctx = shake128_absorb_seed(seed_A);

    for (i = 0; i < SABER_L; i++)
    {

        shake128_inc_squeeze(shake_out, SABER_POLYCOINBYTES, &shake_s_ctx);
        cbd(poly, shake_out);
#ifdef SABER_COMPRESS_SECRETKEY
        POLmu2BS(sk + i * SABER_POLYSECRETBYTES, poly); // sk <- s
#else
        POLq2BS(sk + i * SABER_POLYSECRETBYTES, poly);
#endif

        double_NTT(s_NTTq1, s_NTTq2, poly);

        for (j = 0; j < SABER_L; j++)
        {

            shake128_inc_squeeze(shake_out, SABER_POLYBYTES, &shake_A_ctx);
            BS2POLq(shake_out, poly);

            double_NTT(poly_NTTq1, poly_NTTq2, poly);

            if (i == 0)
            {
                double_basemul(acc_q1 + j * SABER_N, acc_q2 + j * SABER_N, s_NTTq1, s_NTTq2, poly_NTTq1, poly_NTTq2);
            }
            else
            {
                double_basemul(poly_NTTq1, poly_NTTq2, s_NTTq1, s_NTTq2, poly_NTTq1, poly_NTTq2);
                double_add(acc_q1 + j * SABER_N, acc_q2 + j * SABER_N,
                           acc_q1 + j * SABER_N, acc_q2 + j * SABER_N,
                           poly_NTTq1, poly_NTTq2);
            }
        }
    }

    shake128_inc_ctx_release(&shake_A_ctx);
    shake128_inc_ctx_release(&shake_s_ctx);

    for (i = 0; i < SABER_L; i++)
    {

        double_reduce(acc_q1 + i * SABER_N, acc_q2 + i * SABER_N, acc_q1 + i * SABER_N, acc_q2 + i * SABER_N);
        double_iNTT(acc_q1 + i * SABER_N, acc_q2 + i * SABER_N);
        ntt_crt(acc_q1 + i * SABER_N, acc_q2 + i * SABER_N, acc);

        for (j = 0; j < SABER_N; j++)
        {
            poly[j] = ((acc[j] + h1) >> (SABER_EQ - SABER_EP));
        }

        POLp2BS(pk + i * SABER_POLYCOMPRESSEDBYTES, poly);
    }
}

uint32_t MatrixVectorMulEncNTT_A(uint8_t ct0[SABER_POLYVECCOMPRESSEDBYTES], uint8_t ct1[SABER_SCALEBYTES_KEM], const uint8_t seed_s[SABER_NOISE_SEEDBYTES], const uint8_t seed_A[SABER_SEEDBYTES], const uint8_t pk[SABER_INDCPA_PUBLICKEYBYTES], const uint8_t m[SABER_KEYBYTES], int compare)
{

    uint16_t poly[SABER_N];
    uint16_t acc[SABER_N];

    int16_t poly_NTTq1[SABER_N];
    int16_t poly_NTTq2[SABER_N];
    int16_t s_NTTq1[SABER_L * SABER_N];
    int16_t s_NTTq2[SABER_L * SABER_N];
    int16_t acc_q1[SABER_N];
    int16_t acc_q2[SABER_N];

    uint8_t shake_out[MAX(SABER_POLYBYTES, SABER_POLYCOINBYTES)];

    uint16_t *mp = poly;

    size_t i, j;
    uint32_t fail = 0;

    shake128incctx shake_s_ctx = shake128_absorb_seed(seed_s);

    for (i = 0; i < SABER_L; i++)
    {
        shake128_inc_squeeze(shake_out, SABER_POLYCOINBYTES, &shake_s_ctx);
        cbd(poly, shake_out);
        double_NTT(s_NTTq1 + i * SABER_N, s_NTTq2 + i * SABER_N, poly);
    }

    shake128_inc_ctx_release(&shake_s_ctx);

    shake128incctx shake_A_ctx = shake128_absorb_seed(seed_A);

    for (i = 0; i < SABER_L; i++)
    {

        for (j = 0; j < SABER_L; j++)
        {

            shake128_inc_squeeze(shake_out, SABER_POLYBYTES, &shake_A_ctx);
            BS2POLq(shake_out, poly);

            double_NTT(poly_NTTq1, poly_NTTq2, poly);

            if (j == 0)
            {
                double_basemul(acc_q1, acc_q2, s_NTTq1 + j * SABER_N, s_NTTq2 + j * SABER_N, poly_NTTq1, poly_NTTq2);
            }
            else
            {
                double_basemul(poly_NTTq1, poly_NTTq2, s_NTTq1 + j * SABER_N, s_NTTq2 + j * SABER_N, poly_NTTq1, poly_NTTq2);
                double_add(acc_q1, acc_q2, acc_q1, acc_q2, poly_NTTq1, poly_NTTq2);
            }
        }

        double_reduce(acc_q1, acc_q2, acc_q1, acc_q2);
        double_iNTT(acc_q1, acc_q2);
        ntt_crt(acc_q1, acc_q2, acc);

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
        double_NTT(poly_NTTq1, poly_NTTq2, poly);

        if (j == 0)
        {
            double_basemul(acc_q1, acc_q2, s_NTTq1 + j * SABER_N, s_NTTq2 + j * SABER_N, poly_NTTq1, poly_NTTq2);
        }
        else
        {
            double_basemul(poly_NTTq1, poly_NTTq2, s_NTTq1 + j * SABER_N, s_NTTq2 + j * SABER_N, poly_NTTq1, poly_NTTq2);
            double_add(acc_q1, acc_q2, acc_q1, acc_q2, poly_NTTq1, poly_NTTq2);
        }
    }

    BS2POLmsg(m, mp);

    double_reduce(acc_q1, acc_q2, acc_q1, acc_q2);
    double_iNTT(acc_q1, acc_q2);
    ntt_crt(acc_q1, acc_q2, acc);

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

    int16_t poly_NTTq1[SABER_N];
    int16_t poly_NTTq2[SABER_N];
    int16_t buff_NTTq1[SABER_N];
    int16_t buff_NTTq2[SABER_N];
    int16_t acc_q1[SABER_N];
    int16_t acc_q2[SABER_N];

    size_t i;

    for (i = 0; i < SABER_L; i++)
    {

#ifdef SABER_COMPRESS_SECRETKEY
        BS2POLmu(sk + i * SABER_POLYSECRETBYTES, buff);
#else
        BS2POLq(sk + i * SABER_POLYSECRETBYTES, buff);
#endif
        BS2POLp(ciphertext + i * SABER_POLYCOMPRESSEDBYTES, poly);

        double_NTT(buff_NTTq1, buff_NTTq2, buff);
        double_NTT(poly_NTTq1, poly_NTTq2, poly);

        if (i == 0)
        {
            double_basemul(acc_q1, acc_q2, buff_NTTq1, buff_NTTq2, poly_NTTq1, poly_NTTq2);
        }
        else
        {
            double_basemul(poly_NTTq1, poly_NTTq2, buff_NTTq1, buff_NTTq2, poly_NTTq1, poly_NTTq2);
            double_add(acc_q1, acc_q2, acc_q1, acc_q2, poly_NTTq1, poly_NTTq2);
        }
    }

    double_reduce(acc_q1, acc_q2, acc_q1, acc_q2);
    double_iNTT(acc_q1, acc_q2);
    ntt_crt(acc_q1, acc_q2, acc);

    BS2POLT(ciphertext + SABER_POLYVECCOMPRESSEDBYTES, buff);

    for (i = 0; i < SABER_N; i++)
    {
        poly[i] = (acc[i] + h2 - (buff[i] << (SABER_EP - SABER_ET))) >> (SABER_EP - 1);
    }

    POLmsg2BS(m, poly);
}