/* __gmpfr_invsqrt_limb_approx -- reciprocal approximate square root of a limb

Copyright 2017-2019 Free Software Foundation, Inc.
Contributed by the AriC and Caramba projects, INRIA.

This file is part of the GNU MPFR Library.

The GNU MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MPFR Library; see the file COPYING.LESSER.  If not, see
https://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

/* for now, we only provide __gmpfr_invert_limb for 64-bit limb */
#if GMP_NUMB_BITS == 64

/* For 257 <= d10 <= 1024, T[d10-257] = floor(sqrt(2^30/d10)).
   Sage code:
   T = [floor(sqrt(2^30/d10)) for d10 in [257..1024]] */
static const mp_limb_t T[768] =
  {       2044, 2040, 2036, 2032, 2028, 2024, 2020, 2016, 2012, 2009, 2005,
    2001, 1997, 1994, 1990, 1986, 1983, 1979, 1975, 1972, 1968, 1965, 1961,
    1958, 1954, 1951, 1947, 1944, 1941, 1937, 1934, 1930, 1927, 1924, 1920,
    1917, 1914, 1911, 1907, 1904, 1901, 1898, 1895, 1891, 1888, 1885, 1882,
    1879, 1876, 1873, 1870, 1867, 1864, 1861, 1858, 1855, 1852, 1849, 1846,
    1843, 1840, 1837, 1834, 1831, 1828, 1826, 1823, 1820, 1817, 1814, 1812,
    1809, 1806, 1803, 1801, 1798, 1795, 1792, 1790, 1787, 1784, 1782, 1779,
    1777, 1774, 1771, 1769, 1766, 1764, 1761, 1759, 1756, 1754, 1751, 1749,
    1746, 1744, 1741, 1739, 1736, 1734, 1731, 1729, 1727, 1724, 1722, 1719,
    1717, 1715, 1712, 1710, 1708, 1705, 1703, 1701, 1698, 1696, 1694, 1692,
    1689, 1687, 1685, 1683, 1680, 1678, 1676, 1674, 1672, 1670, 1667, 1665,
    1663, 1661, 1659, 1657, 1655, 1652, 1650, 1648, 1646, 1644, 1642, 1640,
    1638, 1636, 1634, 1632, 1630, 1628, 1626, 1624, 1622, 1620, 1618, 1616,
    1614, 1612, 1610, 1608, 1606, 1604, 1602, 1600, 1598, 1597, 1595, 1593,
    1591, 1589, 1587, 1585, 1583, 1582, 1580, 1578, 1576, 1574, 1572, 1571,
    1569, 1567, 1565, 1563, 1562, 1560, 1558, 1556, 1555, 1553, 1551, 1549,
    1548, 1546, 1544, 1542, 1541, 1539, 1537, 1536, 1534, 1532, 1531, 1529,
    1527, 1526, 1524, 1522, 1521, 1519, 1517, 1516, 1514, 1513, 1511, 1509,
    1508, 1506, 1505, 1503, 1501, 1500, 1498, 1497, 1495, 1494, 1492, 1490,
    1489, 1487, 1486, 1484, 1483, 1481, 1480, 1478, 1477, 1475, 1474, 1472,
    1471, 1469, 1468, 1466, 1465, 1463, 1462, 1461, 1459, 1458, 1456, 1455,
    1453, 1452, 1450, 1449, 1448, 1446, 1445, 1443, 1442, 1441, 1439, 1438,
    1436, 1435, 1434, 1432, 1431, 1430, 1428, 1427, 1426, 1424, 1423, 1422,
    1420, 1419, 1418, 1416, 1415, 1414, 1412, 1411, 1410, 1408, 1407, 1406,
    1404, 1403, 1402, 1401, 1399, 1398, 1397, 1395, 1394, 1393, 1392, 1390,
    1389, 1388, 1387, 1385, 1384, 1383, 1382, 1381, 1379, 1378, 1377, 1376,
    1374, 1373, 1372, 1371, 1370, 1368, 1367, 1366, 1365, 1364, 1362, 1361,
    1360, 1359, 1358, 1357, 1355, 1354, 1353, 1352, 1351, 1350, 1349, 1347,
    1346, 1345, 1344, 1343, 1342, 1341, 1339, 1338, 1337, 1336, 1335, 1334,
    1333, 1332, 1331, 1330, 1328, 1327, 1326, 1325, 1324, 1323, 1322, 1321,
    1320, 1319, 1318, 1317, 1315, 1314, 1313, 1312, 1311, 1310, 1309, 1308,
    1307, 1306, 1305, 1304, 1303, 1302, 1301, 1300, 1299, 1298, 1297, 1296,
    1295, 1294, 1293, 1292, 1291, 1290, 1289, 1288, 1287, 1286, 1285, 1284,
    1283, 1282, 1281, 1280, 1279, 1278, 1277, 1276, 1275, 1274, 1273, 1272,
    1271, 1270, 1269, 1268, 1267, 1266, 1265, 1264, 1264, 1263, 1262, 1261,
    1260, 1259, 1258, 1257, 1256, 1255, 1254, 1253, 1252, 1252, 1251, 1250,
    1249, 1248, 1247, 1246, 1245, 1244, 1243, 1242, 1242, 1241, 1240, 1239,
    1238, 1237, 1236, 1235, 1234, 1234, 1233, 1232, 1231, 1230, 1229, 1228,
    1228, 1227, 1226, 1225, 1224, 1223, 1222, 1222, 1221, 1220, 1219, 1218,
    1217, 1216, 1216, 1215, 1214, 1213, 1212, 1211, 1211, 1210, 1209, 1208,
    1207, 1207, 1206, 1205, 1204, 1203, 1202, 1202, 1201, 1200, 1199, 1198,
    1198, 1197, 1196, 1195, 1194, 1194, 1193, 1192, 1191, 1190, 1190, 1189,
    1188, 1187, 1187, 1186, 1185, 1184, 1183, 1183, 1182, 1181, 1180, 1180,
    1179, 1178, 1177, 1177, 1176, 1175, 1174, 1174, 1173, 1172, 1171, 1171,
    1170, 1169, 1168, 1168, 1167, 1166, 1165, 1165, 1164, 1163, 1162, 1162,
    1161, 1160, 1159, 1159, 1158, 1157, 1157, 1156, 1155, 1154, 1154, 1153,
    1152, 1152, 1151, 1150, 1149, 1149, 1148, 1147, 1147, 1146, 1145, 1145,
    1144, 1143, 1142, 1142, 1141, 1140, 1140, 1139, 1138, 1138, 1137, 1136,
    1136, 1135, 1134, 1133, 1133, 1132, 1131, 1131, 1130, 1129, 1129, 1128,
    1127, 1127, 1126, 1125, 1125, 1124, 1123, 1123, 1122, 1121, 1121, 1120,
    1119, 1119, 1118, 1118, 1117, 1116, 1116, 1115, 1114, 1114, 1113, 1112,
    1112, 1111, 1110, 1110, 1109, 1109, 1108, 1107, 1107, 1106, 1105, 1105,
    1104, 1103, 1103, 1102, 1102, 1101, 1100, 1100, 1099, 1099, 1098, 1097,
    1097, 1096, 1095, 1095, 1094, 1094, 1093, 1092, 1092, 1091, 1091, 1090,
    1089, 1089, 1088, 1088, 1087, 1086, 1086, 1085, 1085, 1084, 1083, 1083,
    1082, 1082, 1081, 1080, 1080, 1079, 1079, 1078, 1077, 1077, 1076, 1076,
    1075, 1075, 1074, 1073, 1073, 1072, 1072, 1071, 1071, 1070, 1069, 1069,
    1068, 1068, 1067, 1067, 1066, 1065, 1065, 1064, 1064, 1063, 1063, 1062,
    1062, 1061, 1060, 1060, 1059, 1059, 1058, 1058, 1057, 1057, 1056, 1055,
    1055, 1054, 1054, 1053, 1053, 1052, 1052, 1051, 1051, 1050, 1049, 1049,
    1048, 1048, 1047, 1047, 1046, 1046, 1045, 1045, 1044, 1044, 1043, 1043,
    1042, 1041, 1041, 1040, 1040, 1039, 1039, 1038, 1038, 1037, 1037, 1036,
    1036, 1035, 1035, 1034, 1034, 1033, 1033, 1032, 1032, 1031, 1031, 1030,
    1030, 1029, 1029, 1028, 1028, 1027, 1027, 1026, 1026, 1025, 1025, 1024,
    1024 };

/* table of v0^3 */
static const mp_limb_t T3[768] =
  {             8539701184, 8489664000, 8439822656, 8390176768, 8340725952,
    8291469824, 8242408000, 8193540096, 8144865728, 8108486729, 8060150125,
    8012006001, 7964053973, 7928215784, 7880599000, 7833173256, 7797729087,
    7750636739, 7703734375, 7668682048, 7622111232, 7587307125, 7541066681,
    7506509912, 7460598664, 7426288351, 7380705123, 7346640384, 7312680621,
    7267563953, 7233848504, 7189057000, 7155584983, 7122217024, 7077888000,
    7044762213, 7011739944, 6978821031, 6935089643, 6902411264, 6869835701,
    6837362792, 6804992375, 6761990971, 6729859072, 6697829125, 6665900968,
    6634074439, 6602349376, 6570725617, 6539203000, 6507781363, 6476460544,
    6445240381, 6414120712, 6383101375, 6352182208, 6321363049, 6290643736,
    6260024107, 6229504000, 6199083253, 6168761704, 6138539191, 6108415552,
    6088387976, 6058428767, 6028568000, 5998805513, 5969141144, 5949419328,
    5919918129, 5890514616, 5861208627, 5841725401, 5812581592, 5783534875,
    5754585088, 5735339000, 5706550403, 5677858304, 5658783768, 5630252139,
    5611284433, 5582912824, 5554637011, 5535839609, 5507723096, 5489031744,
    5461074081, 5442488479, 5414689216, 5396209064, 5368567751, 5350192749,
    5322708936, 5304438784, 5277112021, 5258946419, 5231776256, 5213714904,
    5186700891, 5168743489, 5150827583, 5124031424, 5106219048, 5079577959,
    5061868813, 5044200875, 5017776128, 5000211000, 4982686912, 4956477625,
    4939055927, 4921675101, 4895680392, 4878401536, 4861163384, 4843965888,
    4818245769, 4801149703, 4784094125, 4767078987, 4741632000, 4724717752,
    4707843776, 4691010024, 4674216448, 4657463000, 4632407963, 4615754625,
    4599141247, 4582567781, 4566034179, 4549540393, 4533086375, 4508479808,
    4492125000, 4475809792, 4459534136, 4443297984, 4427101288, 4410944000,
    4394826072, 4378747456, 4362708104, 4346707968, 4330747000, 4314825152,
    4298942376, 4283098624, 4267293848, 4251528000, 4235801032, 4220112896,
    4204463544, 4188852928, 4173281000, 4157747712, 4142253016, 4126796864,
    4111379208, 4096000000, 4080659192, 4073003173, 4057719875, 4042474857,
    4027268071, 4012099469, 3996969003, 3981876625, 3966822287, 3959309368,
    3944312000, 3929352552, 3914430976, 3899547224, 3884701248, 3877292411,
    3862503009, 3847751263, 3833037125, 3818360547, 3811036328, 3796416000,
    3781833112, 3767287616, 3760028875, 3745539377, 3731087151, 3716672149,
    3709478592, 3695119336, 3680797184, 3666512088, 3659383421, 3645153819,
    3630961153, 3623878656, 3609741304, 3595640768, 3588604291, 3574558889,
    3560550183, 3553559576, 3539605824, 3525688648, 3518743761, 3504881359,
    3491055413, 3484156096, 3470384744, 3463512697, 3449795831, 3436115229,
    3429288512, 3415662216, 3408862625, 3395290527, 3381754501, 3375000000,
    3361517992, 3354790473, 3341362375, 3334661784, 3321287488, 3307949000,
    3301293169, 3288008303, 3281379256, 3268147904, 3261545587, 3248367641,
    3241792000, 3228667352, 3222118333, 3209046875, 3202524424, 3189506048,
    3183010111, 3170044709, 3163575232, 3150662696, 3144219625, 3131359847,
    3124943128, 3118535181, 3105745579, 3099363912, 3086626816, 3080271375,
    3067586677, 3061257408, 3048625000, 3042321849, 3036027392, 3023464536,
    3017196125, 3004685307, 2998442888, 2992209121, 2979767519, 2973559672,
    2961169856, 2954987875, 2948814504, 2936493568, 2930345991, 2924207000,
    2911954752, 2905841483, 2899736776, 2887553024, 2881473967, 2875403448,
    2863288000, 2857243059, 2851206632, 2839159296, 2833148375, 2827145944,
    2815166528, 2809189531, 2803221000, 2791309312, 2785366143, 2779431416,
    2767587264, 2761677827, 2755776808, 2749884201, 2738124199, 2732256792,
    2726397773, 2714704875, 2708870984, 2703045457, 2697228288, 2685619000,
    2679826869, 2674043072, 2668267603, 2656741625, 2650991104, 2645248887,
    2639514968, 2633789341, 2622362939, 2616662152, 2610969633, 2605285376,
    2593941624, 2588282117, 2582630848, 2576987811, 2571353000, 2560108032,
    2554497863, 2548895896, 2543302125, 2537716544, 2526569928, 2521008881,
    2515456000, 2509911279, 2504374712, 2498846293, 2487813875, 2482309864,
    2476813977, 2471326208, 2465846551, 2460375000, 2454911549, 2444008923,
    2438569736, 2433138625, 2427715584, 2422300607, 2416893688, 2411494821,
    2400721219, 2395346472, 2389979753, 2384621056, 2379270375, 2373927704,
    2368593037, 2363266368, 2357947691, 2352637000, 2342039552, 2336752783,
    2331473976, 2326203125, 2320940224, 2315685267, 2310438248, 2305199161,
    2299968000, 2294744759, 2289529432, 2284322013, 2273930875, 2268747144,
    2263571297, 2258403328, 2253243231, 2248091000, 2242946629, 2237810112,
    2232681443, 2227560616, 2222447625, 2217342464, 2212245127, 2207155608,
    2202073901, 2197000000, 2191933899, 2186875592, 2181825073, 2176782336,
    2171747375, 2166720184, 2161700757, 2156689088, 2151685171, 2146689000,
    2141700569, 2136719872, 2131746903, 2126781656, 2121824125, 2116874304,
    2111932187, 2106997768, 2102071041, 2097152000, 2092240639, 2087336952,
    2082440933, 2077552576, 2072671875, 2067798824, 2062933417, 2058075648,
    2053225511, 2048383000, 2043548109, 2038720832, 2033901163, 2029089096,
    2024284625, 2019487744, 2019487744, 2014698447, 2009916728, 2005142581,
    2000376000, 1995616979, 1990865512, 1986121593, 1981385216, 1976656375,
    1971935064, 1967221277, 1962515008, 1962515008, 1957816251, 1953125000,
    1948441249, 1943764992, 1939096223, 1934434936, 1929781125, 1925134784,
    1920495907, 1915864488, 1915864488, 1911240521, 1906624000, 1902014919,
    1897413272, 1892819053, 1888232256, 1883652875, 1879080904, 1879080904,
    1874516337, 1869959168, 1865409391, 1860867000, 1856331989, 1851804352,
    1851804352, 1847284083, 1842771176, 1838265625, 1833767424, 1829276567,
    1824793048, 1824793048, 1820316861, 1815848000, 1811386459, 1806932232,
    1802485313, 1798045696, 1798045696, 1793613375, 1789188344, 1784770597,
    1780360128, 1775956931, 1775956931, 1771561000, 1767172329, 1762790912,
    1758416743, 1758416743, 1754049816, 1749690125, 1745337664, 1740992427,
    1736654408, 1736654408, 1732323601, 1728000000, 1723683599, 1719374392,
    1719374392, 1715072373, 1710777536, 1706489875, 1702209384, 1702209384,
    1697936057, 1693669888, 1689410871, 1685159000, 1685159000, 1680914269,
    1676676672, 1672446203, 1672446203, 1668222856, 1664006625, 1659797504,
    1655595487, 1655595487, 1651400568, 1647212741, 1643032000, 1643032000,
    1638858339, 1634691752, 1630532233, 1630532233, 1626379776, 1622234375,
    1618096024, 1618096024, 1613964717, 1609840448, 1605723211, 1605723211,
    1601613000, 1597509809, 1593413632, 1593413632, 1589324463, 1585242296,
    1581167125, 1581167125, 1577098944, 1573037747, 1568983528, 1568983528,
    1564936281, 1560896000, 1556862679, 1556862679, 1552836312, 1548816893,
    1548816893, 1544804416, 1540798875, 1536800264, 1536800264, 1532808577,
    1528823808, 1528823808, 1524845951, 1520875000, 1516910949, 1516910949,
    1512953792, 1509003523, 1509003523, 1505060136, 1501123625, 1501123625,
    1497193984, 1493271207, 1489355288, 1489355288, 1485446221, 1481544000,
    1481544000, 1477648619, 1473760072, 1473760072, 1469878353, 1466003456,
    1466003456, 1462135375, 1458274104, 1454419637, 1454419637, 1450571968,
    1446731091, 1446731091, 1442897000, 1439069689, 1439069689, 1435249152,
    1431435383, 1431435383, 1427628376, 1423828125, 1423828125, 1420034624,
    1416247867, 1416247867, 1412467848, 1408694561, 1408694561, 1404928000,
    1401168159, 1401168159, 1397415032, 1397415032, 1393668613, 1389928896,
    1389928896, 1386195875, 1382469544, 1382469544, 1378749897, 1375036928,
    1375036928, 1371330631, 1367631000, 1367631000, 1363938029, 1363938029,
    1360251712, 1356572043, 1356572043, 1352899016, 1349232625, 1349232625,
    1345572864, 1341919727, 1341919727, 1338273208, 1338273208, 1334633301,
    1331000000, 1331000000, 1327373299, 1327373299, 1323753192, 1320139673,
    1320139673, 1316532736, 1312932375, 1312932375, 1309338584, 1309338584,
    1305751357, 1302170688, 1302170688, 1298596571, 1298596571, 1295029000,
    1291467969, 1291467969, 1287913472, 1287913472, 1284365503, 1280824056,
    1280824056, 1277289125, 1277289125, 1273760704, 1270238787, 1270238787,
    1266723368, 1266723368, 1263214441, 1259712000, 1259712000, 1256216039,
    1256216039, 1252726552, 1249243533, 1249243533, 1245766976, 1245766976,
    1242296875, 1242296875, 1238833224, 1235376017, 1235376017, 1231925248,
    1231925248, 1228480911, 1228480911, 1225043000, 1221611509, 1221611509,
    1218186432, 1218186432, 1214767763, 1214767763, 1211355496, 1207949625,
    1207949625, 1204550144, 1204550144, 1201157047, 1201157047, 1197770328,
    1197770328, 1194389981, 1191016000, 1191016000, 1187648379, 1187648379,
    1184287112, 1184287112, 1180932193, 1180932193, 1177583616, 1174241375,
    1174241375, 1170905464, 1170905464, 1167575877, 1167575877, 1164252608,
    1164252608, 1160935651, 1160935651, 1157625000, 1154320649, 1154320649,
    1151022592, 1151022592, 1147730823, 1147730823, 1144445336, 1144445336,
    1141166125, 1141166125, 1137893184, 1137893184, 1134626507, 1134626507,
    1131366088, 1128111921, 1128111921, 1124864000, 1124864000, 1121622319,
    1121622319, 1118386872, 1118386872, 1115157653, 1115157653, 1111934656,
    1111934656, 1108717875, 1108717875, 1105507304, 1105507304, 1102302937,
    1102302937, 1099104768, 1099104768, 1095912791, 1095912791, 1092727000,
    1092727000, 1089547389, 1089547389, 1086373952, 1086373952, 1083206683,
    1083206683, 1080045576, 1080045576, 1076890625, 1076890625, 1073741824,
    1073741824 };

/* umul_hi(h, x, y) puts in h the high part of x*y */
#ifdef HAVE_MULX_U64
#include <immintrin.h>
#define umul_hi(h, x, y) _mulx_u64 (x, y, (unsigned long long *) &h)
#else
#define umul_hi(h, x, y)                        \
  do {                                          \
    mp_limb_t _l;                               \
    umul_ppmm (h, _l, x, y);                    \
  } while (0)
#endif

/* given 2^62 <= d < 2^64, put in r an approximation of
   s = floor(2^96/sqrt(r)) - 2^64, with r <= s <= r + 15 */
#define __gmpfr_invsqrt_limb_approx(r, d)                               \
  do {                                                                  \
    mp_limb_t _d, _i, _v0, _e0, _d37, _v1, _e1, _h, _v2, _e2;           \
    _d = (d);                                                           \
    _i = (_d >> 54) - 256;                                              \
    /* i = d10 - 256 */                                                 \
    _v0 = T[_i];                                                        \
    _d37 = 1 + (_d >> 27);                                              \
    _e0 = T3[_i] * _d37;                                                \
    /* the value (_v0 << 57) - _e0 is less than 2^61 */                 \
    _v1 = (_v0 << 11) + (((_v0 << 57) - _e0) >> 47);                    \
    _e1 = - _v1 * _v1 * _d37;                                           \
    umul_hi (_h, _v1, _e1);                                             \
    /* _h = floor(e_1*v_1/2^64) */                                      \
    _v2 = (_v1 << 10) + (_h >> 6);                                      \
    umul_hi (_h, _v2 * _v2, _d);                                        \
    /* in _h + 2, one +1 accounts for the lower neglected part of */    \
    /* v2^2*d. the other +1 is to compute ceil((h+1)/2) */              \
    _e2 = (MPFR_LIMB_ONE << 61) - ((_h + 2) >> 1);                      \
    _h = _v2 * _e2;                                                     \
    (r) = (_v2 << 33) + (_h >> 29);                                     \
  } while (0)

/* given 2^62 <= d < 2^64, return a 32-bit approximation r of
   sqrt(2^126/d) */
#define __gmpfr_invsqrt_halflimb_approx(r, d)                           \
  do {                                                                  \
    mp_limb_t _d, _i, _v0, _e0, _d37, _v1, _e1, _h;                     \
    _d = (d);                                                           \
    _i = (_d >> 54) - 256;                                              \
    /* i = d10 - 256 */                                                 \
    _v0 = T[_i];                                                        \
    _d37 = 1 + (_d >> 27);                                              \
    _e0 = T3[_i] * _d37;                                                \
    /* the value (_v0 << 57) - _e0 is less than 2^61 */                 \
    _v1 = (_v0 << 11) + (((_v0 << 57) - _e0) >> 47);                    \
    _e1 = - _v1 * _v1 * _d37;                                           \
    umul_hi (_h, _v1, _e1);                                             \
    /* _h = floor(e_1*v_1/2^64) */                                      \
    (r) = (_v1 << 10) + (_h >> 6);                                      \
  } while (0)

/* given 2^62 <= n < 2^64, put in s an approximation of sqrt(2^64*n),
   with: s <= floor(sqrt(2^64*u)) <= s + 7 */
#define __gmpfr_sqrt_limb_approx(s, n)                                  \
  do {                                                                  \
    mp_limb_t _n, _x, _y, _z, _t;                                       \
    _n = (n);                                                           \
    __gmpfr_invsqrt_halflimb_approx (_x, _n);                           \
    MPFR_ASSERTD(_x < MPFR_LIMB_ONE << 32);                             \
    /* x has 32 bits, and is near (by below) sqrt(2^126/n) */           \
    _y = (_x * (_n >> 31)) >> 32;                                       \
    MPFR_ASSERTD(_y < MPFR_LIMB_ONE << 32);                             \
    /* y is near (by below) sqrt(n) */                                  \
    _z = _n - _y * _y;                                                  \
    /* reduce _z so that _z <= 2*_y */                                  \
    /* the maximal value of _z is 2*(2^32-1) */                         \
    while (_z > 2 * ((MPFR_LIMB_ONE << 32) - 1))                        \
      {                                                                 \
        _z -= (_y + _y) + 1;                                            \
        _y ++;                                                          \
      }                                                                 \
    /* now z <= 2*(2^32-1): one reduction is enough */                  \
    if (_z > _y + _y)                                                   \
      {                                                                 \
        _z -= (_y + _y) + 1;                                            \
        _y ++;                                                          \
      }                                                                 \
    /* _x * _z should be < 2^64 */                                      \
    _t = (_x * _z) >> 32;                                               \
    (s) = (_y << 32) + _t;                                              \
  } while (0)

/* given 2^62 <= u < 2^64, put in s the value floor(sqrt(2^64*u)),
   and in rh in rl the remainder: 2^64*u - s^2 = 2^64*rh + rl, with
   2^64*rh + rl <= 2*s, and in invs the approximation of 2^96/sqrt(u) */
#define __gmpfr_sqrt_limb(s, rh, rl, invs, u)                   \
  do {                                                          \
    mp_limb_t _u, _invs, _r, _h, _l;                            \
    _u = (u);                                                   \
    __gmpfr_invsqrt_limb_approx (_invs, _u);                    \
    umul_hi (_h, _invs, _u);                                    \
    _r = _h + _u;                                               \
    /* make sure _r has its most significant bit set */         \
    if (MPFR_UNLIKELY(_r < MPFR_LIMB_HIGHBIT))                  \
      _r = MPFR_LIMB_HIGHBIT;                                   \
    /* we know _r <= sqrt(2^64*u) <= _r + 16 */                 \
    umul_ppmm (_h, _l, _r, _r);                                 \
    sub_ddmmss (_h, _l, _u, 0, _h, _l);                         \
    /* now h:l <= 30*_r */                                      \
    MPFR_ASSERTD(_h < 30);                                      \
    if (_h >= 16)                                               \
      { /* subtract 16r+64 to h:l, add 8 to _r */               \
        sub_ddmmss (_h, _l, _h, _l, _r >> 60, _r << 4);         \
        sub_ddmmss (_h, _l, _h, _l, 0, 64);                     \
        _r += 8;                                                \
      }                                                         \
    if (_h >= 8)                                                \
      { /* subtract 8r+16 to h:l, add 4 to _r */                \
        sub_ddmmss (_h, _l, _h, _l, _r >> 61, _r << 3);         \
        sub_ddmmss (_h, _l, _h, _l, 0, 16);                     \
        _r += 4;                                                \
      }                                                         \
    if (_h >= 4)                                                \
      { /* subtract 4r+4 to h:l, add 2 to _r */                 \
        sub_ddmmss (_h, _l, _h, _l, _r >> 62, _r << 2);         \
        sub_ddmmss (_h, _l, _h, _l, 0, 4);                      \
        _r += 2;                                                \
      }                                                         \
    while (_h > 1 || ((_h == 1) && (_l > 2 * _r)))              \
      { /* subtract 2r+1 to h:l, add 1 to _r */                 \
        sub_ddmmss (_h, _l, _h, _l, _r >> 63, (_r << 1) + 1);   \
        _r ++;                                                  \
      }                                                         \
    (s) = _r;                                                   \
    (rh) = _h;                                                   \
    (rl) = _l;                                                   \
    (invs) = _invs;                                             \
  } while (0)

#endif /* GMP_NUMB_BITS == 64 */
