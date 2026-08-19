
//RAIL
//height of rail:
H_R = 0.01;
//width of rail:
W_R = 0.01;
//length of rail:
L = 0.2;

//AIR
//width of air:
W = 0.08;
//height of air:
H = 0.08;

//air gap:
A = 0.002;

//CORE
//height of core:
H_I = 0.002;
//length of core:
L_I = 0.01;

//width coil:
W_C = 0.002;

//x-gap coil:
G_X = 0.002;

//factor for size:
fact = 0.6;

//skin depth:
SKIND = 0.12;

//Abstand von Mitte weg:
P = -0.03;

//Abstand von oberer Eisenkernkante weg:
DIST_Y = 0.001;
//Abstand von hinterer Eisenkernkante weg:
DIST_X = 0.001;


//fuer feine Bereiche
k = 12;
//fue grobe horizontale Bereiche
q_h = 8;
//fue grobe vertikale Bereiche
q_v = q_h;
//fuer air-gap Bereich
//a = INDEX3_;

//Mesh.Algorithm = 8;
//Mesh.Algorithm3D = 1;
//Mesh.RecombinationAlgorithm = 1;
//Mesh.RecombinationAlgorithm = 3;
//Mesh.RecombineAll =1;
//Mesh.SubdivisionAlgorithm = 2;

Mesh.MshFileVersion = 2.2;
//Mesh.PartitionOldStyleMsh2 = 1;
//Mesh.PartitionCreateGhostCells = 1;

//Mesh.MeshSizeMin = 0.0000001;
//Mesh.MeshSizeMax = 10;

//Mesh.MeshSizeFromPoints = 0;
//Mesh.MeshSizeFromCurvature = 0;
//Mesh.MeshSizeExtendFromBoundary = 0;

//Points in z-direction (should not be changed)
N_z = 11;
//Points in y-direction (should not be changed)
N_y = 6;

//If N_z and N_y are being changed, these two lists have to be changed accordingly.
//These lists have the y-, and z-coordinates of the symmetry plane of all nodes
// mylist_zdirection have N_z entries
// mylist_ydirection have N_y entries
mylist_zdirection[] = {-L/2, -L_I/2-G_X-W_C+P, -L_I/2-G_X+P, -L_I/2+P, -L_I/2+L_I/4+P, 0+P, L_I/4+P, L_I/2+P, L_I/2+G_X+P, L_I/2+G_X+W_C+P, L/2};
mylist_ydirection[] = {0,H_R/2,H_R/2+A,H_R/2+A+H_I,H_R/2+A+H_I+DIST_Y,H/2}; 

//Create all Points of the symmetry plane
//if N_z is changing, this loop is not to changed
//if N_y is changing, the number of rows in the loop has to change, number of rows = N_y
For i In {0:N_z-1}
  Point(1+i*N_y) = {0, mylist_ydirection(0), mylist_zdirection(i)};
  Point(2+i*N_y) = {0, mylist_ydirection(1), mylist_zdirection(i)};
  Point(3+i*N_y) = {0, mylist_ydirection(2), mylist_zdirection(i)};
  Point(4+i*N_y) = {0, mylist_ydirection(3), mylist_zdirection(i)};
  Point(5+i*N_y) = {0, mylist_ydirection(4), mylist_zdirection(i)};
  Point(6+i*N_y) = {0, mylist_ydirection(5), mylist_zdirection(i)};
EndFor

//Create all lines of the symmetry plane
For j In {0:(N_z-1)}
  For i In {0:(N_y-2)}
    Line(1+i+(N_y-1)*j) = {i+N_y*j+1, i+N_y*j+2};
  EndFor
EndFor

For j In {0:(N_y-1)}
  For i In {0:(N_z-2)}
    If (j == 0)
      Line((1+(N_y-2)+(N_y-1)*(N_z-1))+1+i) = {(N_y*i+1),(N_y*i+1+N_y)};
    Else
      Line((1+(N_y-2)+(N_y-1)*(N_z-1))+1+i+(N_z-1)*j) = {(N_y*i+1*(j+1)),(N_y*i+1*(j+1)+N_y)};
    EndIf
  EndFor
EndFor

//Create all surfaces of the symmetry plane
For j In {0:(N_y-2)}
  For i In {0:(N_z-2)}
    If (j == 0)
      Curve Loop(i+1) = {(N_y-1)*i+1, N_y*N_z+i, -( (N_y-1)*(i)+N_y ), -(N_y*N_z-(N_z-1)+i)};
      Plane Surface(i+1) = {i+1};
    Else
      Curve Loop(j*(N_z-1)+i+1) = {(N_y-1)*i+1+j, N_y*N_z+i+(N_z-1)*j, -( (N_y-1)*(i)+N_y+j ), -(N_y*N_z-(N_z-1)+i+(N_z-1)*j)};
      Plane Surface(j*(N_z-1)+i+1) = {j*(N_z-1)+i+1};
    EndIf
  EndFor
EndFor

Coherence;


//outer areas z-direction hinten
Transfinite Curve {106, 96, 86, 76, 66, 56} = ((L/2+P-(L_I/2+G_X+W_C))*(2*k/W_R))/1 Using Progression 1;
//outer areas z-direction vorne
Transfinite Curve {115, 105, 95, 85, 75, 65} = ((L/2-P-(L_I/2+G_X+W_C))*(2*k/W_R))/1 Using Progression 1;
//inner areas z-direction
Transfinite Curve {114, 104, 94, 84, 74, 64, 107, 97, 87, 77, 67, 57} = ((W_C)*(2*k/W_R))+2 Using Progression 1;
//post_processing areas z-direction
Transfinite Curve {63, 73, 83, 93, 103, 113, 62, 72, 82, 92, 102, 112, 61, 71, 81, 91, 101, 111, 60, 70, 80, 90, 100, 110, 59, 69, 79, 89, 99, 109, 58, 68, 78, 88, 98, 108} = ((L_I/4)*(2*k/W_R))+2 Using Progression 1;


//outer areas y-direction
Transfinite Curve {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55} = q_v+1 Using Progression 1;
//inner areas y-direction
Transfinite Curve {4, 9, 14, 19, 24, 29, 34, 39, 44, 49, 54} = ((DIST_Y)*(2*k/W_R))+1 Using Progression 1;
//coil-core areas y-direction
Transfinite Curve {53, 48, 43, 38, 33, 28, 23, 18, 13, 8, 3} = ((H_I)*(2*k/W_R))+2 Using Progression 1;
//Air-Gap areas y-direction
Transfinite Curve {52, 47, 42, 37, 32, 27, 22, 17, 12, 7, 2} = ((A)*(2*k/W_R))+1 Using Progression 1;
//Rail areas y-direction
Transfinite Curve {51, 46, 41, 36, 31, 26, 21, 16, 11, 6, 1} = k+1 Using Progression 1;


Extrude {W_R/2, 0, 0} {
  Surface{1:((N_y-2)*(N_z-1)+(N_z-2)+1)}; Layers {k}; Recombine;
}
Coherence;
eps = 1e-6;
Surfaces_1() = Surface In BoundingBox{W_R/2-eps, 0-eps, -L/2-eps, W_R/2+eps, H/2+eps, L/2+eps};

For i In {0:#Surfaces_1()-1}

Extrude {DIST_X, 0, 0} {
   Surface{Surfaces_1(i)}; Layers {((DIST_X)*(2*k/W_R))+1}; Recombine;
}
EndFor
Coherence;

Surfaces_2() = Surface In BoundingBox{W_R/2+DIST_X-eps, 0-eps, -L/2-eps, W_R/2+DIST_X+eps, H/2+eps, L/2+eps};

For i In {0:#Surfaces_1()-1}
Extrude {W/2-DIST_X-W_R/2, 0, 0} {
   Surface{Surfaces_2(i)}; Layers {q_h}; Recombine;
}
EndFor
Coherence;

Physical Volume("Rail") = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
Physical Volume("Core") = {24, 25, 26, 27};
Physical Volume("Coil1") = {122, 72, 22};
Physical Volume("Coil2") = {129, 79, 29};
Physical Volume("Air_Coarse") = {150, 110, 120, 130, 140, 119, 109, 139, 149, 148, 118, 128, 108, 138, 107, 117, 127, 137, 147, 106, 136, 126, 146, 116, 105, 135, 125, 145, 115, 124, 114, 104, 134, 144, 123, 113, 143, 103, 133, 142, 102, 132, 112, 101, 121, 111, 141, 131, 100, 50, 99, 49, 48, 98, 47, 97, 96, 46, 95, 45, 94, 44, 93, 43, 92, 42, 91, 41};
Physical Volume("Air_Fine") = {40, 90, 30, 80, 20, 70, 60, 39, 89, 19, 69, 59, 38, 88, 28, 78, 18, 68, 58, 37, 87, 77, 17, 67, 57, 36, 86, 76, 16, 66, 56, 35, 85, 75, 15, 65, 55, 34, 84, 74, 14, 64, 54, 33, 83, 23, 73, 13, 63, 53, 32, 82, 12, 62, 52, 31, 81, 21, 71, 11, 61, 51};


Physical Surface("schiene_vorne") = {330};
Physical Surface("schiene_hinten") = {124};
Physical Surface("luft_aussen") = {2310, 3410, 1210, 3415, 1215, 2315, 50, 2306, 2288, 3388, 3406, 1206, 1188, 3393, 2293, 1193, 49, 2284, 3384, 1166, 1184, 2266, 3366, 2271, 3371, 1171, 48, 1162, 2244, 2262, 3362, 1144, 3344, 3349, 1149, 2249, 47, 2240, 2222, 3322, 1140, 3340, 1122, 2227, 1127, 46, 3327, 1100, 3300, 3318, 1118, 2200, 2218, 1105, 3305, 45, 2205, 2196, 2178, 1096, 1078, 3296, 3278, 44, 1083, 3283, 2183, 2156, 2174, 3274, 1074, 3256, 1056, 3261, 43, 2161, 1061, 3252, 2152, 2134, 3234, 1034, 1052, 2139, 1039, 42, 3239, 3212, 3230, 1012, 1030, 2112, 2130, 2117, 3217, 1017, 41, 2108, 2104, 1008, 1004, 3204, 3208, 1004, 2104, 3204, 784, 1884, 2984, 1224, 344, 564, 1444, 1664, 2324, 2544, 2764, 1210, 2310, 3410, 990, 2090, 3190, 770, 550, 1650, 1430, 2530, 2750, 2970, 1870, 2750, 2530, 2970, 3410, 3190, 2746, 3186, 3406, 2966, 2526, 2534, 2728, 2535, 2508, 2975, 2948, 3168, 3195, 2755, 3388, 3415, 2504, 3164, 3384, 2944, 2724, 2512, 2926, 3146, 2733, 2706, 3173, 2513, 3366, 2486, 3393, 3362, 2490, 3142, 2702, 2482, 2922, 3151, 3371, 3344, 2711, 3124, 2931, 2904, 2684, 2491, 2464, 2900, 3340, 3120, 2680, 2460, 2468, 2909, 3129, 2469, 3349, 2689, 2882, 2442, 2662, 3102, 3322, 2658, 3318, 2878, 2438, 2446, 3098, 3300, 2887, 2860, 3327, 2667, 2640, 2420, 2447, 3107, 3080, 3296, 2636, 2424, 2416, 2856, 3076, 3085, 2645, 2425, 2618, 2398, 3058, 3278, 3305, 2865, 2838, 2614, 2402, 2394, 2834, 3054, 3274, 2596, 2623, 3256, 2816, 3063, 2843, 3283, 3036, 2376, 2403, 2380, 2592, 2372, 3032, 3252, 2812, 3041, 3261, 2821, 2601, 2381, 2794, 3014, 2574, 2354, 3234, 3230, 2790, 2350, 3010, 2358, 2570, 3239, 3212, 2772, 2552, 2579, 3019, 2359, 2332, 2992, 3208, 2768, 2548, 2328, 2988, 2336, 3217, 2557, 2544, 2997, 2984, 2337, 2777, 2764, 3204, 2324};
Physical Surface("luft_vertikal_sym") = {1210, 550, 770, 990, 30, 40, 20, 50, 326, 986, 968, 766, 748, 1188, 1206, 546, 528, 19, 39, 946, 506, 524, 964, 304, 726, 49, 1166, 744, 1184, 704, 502, 484, 722, 942, 924, 38, 28, 48, 18, 1162, 282, 1144, 37, 47, 17, 480, 260, 462, 920, 1140, 1122, 700, 902, 678, 46, 16, 458, 238, 1118, 1100, 36, 898, 880, 440, 216, 15, 418, 436, 1096, 1078, 35, 45, 876, 858, 656, 44, 194, 34, 396, 414, 616, 634, 836, 854, 1056, 1074, 14, 33, 13, 23, 43, 832, 814, 594, 612, 392, 172, 1052, 374, 1034, 792, 32, 810, 1012, 1030, 370, 150, 12, 42, 352, 572, 590, 128, 21, 41, 31, 784, 788, 568, 564, 1008, 11, 348, 344, 1004};
Physical Surface("luft_horizontal_sym") = {2530, 1430, 1434, 1408, 1435, 2535, 335, 2508, 2534, 2512, 2513, 1412, 1413, 313, 2486, 1386, 1364, 1391, 2491, 291, 2490, 1390, 2464, 1369, 269, 2469, 1368, 1342, 2442, 2468, 247, 1347, 2447, 2446, 1346, 1320, 2420, 225, 1325, 2398, 1324, 1298, 2424, 2425, 1276, 1302, 1303, 203, 2403, 2402, 2376, 2381, 1281, 181, 2354, 2380, 1280, 1254, 1259, 2359, 1258, 1232, 159, 2358, 2332, 137, 1237, 1236, 1224, 2337, 2324, 2336};

Physical Surface("luft_horizontal_sym") -= {2530, 1430, 335, 2535, 2508, 1435, 1408, 1413, 2486, 2513, 1386, 313, 2491, 291, 1391, 2464, 1364, 269, 1369, 2469, 1342, 2442, 1347, 247, 2447, 1320, 2420, 1298, 1325, 2398, 2425, 225, 2376, 1303, 1276, 2403, 203, 1281, 181, 2381, 2354, 1254, 2359, 159, 1232, 1259, 2332, 1224, 1237, 137, 2324, 2337};
Physical Surface("luft_aussen") -= {1210, 2310, 3410, 3415, 50, 1215, 2315, 1188, 2288, 3388, 3393, 49, 1193, 2293, 3366, 1166, 2266, 3371, 1171, 2271, 48, 2244, 1144, 3344, 3349, 47, 2249, 1149, 2222, 3322, 1122, 2227, 3327, 46, 1127, 1100, 3300, 2200, 45, 2205, 1105, 3305, 2178, 1078, 3278, 44, 2183, 3283, 1083, 2156, 3256, 1056, 2161, 3261, 43, 1061, 2134, 1034, 3234, 3239, 1039, 2139, 42, 3212, 2112, 1012, 3217, 41, 1017, 2117, 1004, 2104, 3204, 2799, 2953, 3410, 2530, 2750, 2970, 3190, 2948, 2534, 2728, 2526, 2746, 3388, 2508, 2966, 3406, 3186, 3168, 2486, 3164, 2504, 2512, 3384, 2706, 2724, 3366, 3146, 2926, 2944, 2702, 2490, 3142, 3124, 2464, 2482, 3344, 2922, 2904, 3362, 2684, 3102, 2468, 2460, 2442, 3340, 2900, 2882, 3120, 2662, 2680, 3322, 2878, 2860, 2640, 2658, 3318, 3300, 3098, 3080, 2420, 2438, 2446, 2636, 2618, 2856, 2838, 3076, 3058, 2398, 2416, 3296, 2424, 3278, 3036, 2816, 2834, 3274, 3054, 3256, 2402, 2614, 2596, 2376, 2394, 2372, 3032, 2380, 3014, 2354, 3252, 3234, 2794, 2812, 2592, 2574, 2332, 2992, 3230, 3212, 3010, 2772, 2350, 2790, 2358, 2570, 2552, 3208, 3204, 2988, 2984, 2768, 2764, 2548, 2544, 2336, 2328, 2324};
Physical Surface("luft_vertikal_sym") -= {1210, 990, 770, 550, 986, 968, 1206, 1188, 528, 546, 766, 748, 326, 726, 506, 946, 304, 964, 744, 524, 1184, 1166, 704, 722, 1162, 1144, 282, 484, 924, 502, 942, 462, 260, 700, 1140, 1122, 480, 920, 902, 678, 458, 440, 238, 880, 898, 1118, 1100, 858, 656, 876, 1096, 1078, 418, 436, 216, 836, 854, 1056, 1074, 194, 414, 396, 634, 616, 1034, 1052, 594, 832, 814, 172, 374, 392, 612, 1030, 1012, 572, 150, 810, 792, 590, 370, 352, 564, 568, 348, 344, 784, 788, 128, 1008, 1004};

Physical Surface("schiene_vertikal_sym") = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
Physical Surface("schiene_horizontal_sym") = {136, 158, 180, 202, 224, 246, 268, 312, 290, 334};
Physical Surface("Kern_vertikal_sym") = {24, 25, 26, 27};
//Physical Surface("Spule_vertikal_sym") = {43, 39};

Physical Surface("Spule1_ground") = {2799};
Physical Surface("Spule1_voltage") = {22};
Physical Surface("Spule2_ground") = {2953};
Physical Surface("Spule2_voltage") = {29};

Physical Surface("flux") = {414, 436, 458, 480};

Physical Surface("current") = {4, 5, 6, 7};

Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";