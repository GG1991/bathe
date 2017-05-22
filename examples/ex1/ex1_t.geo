lc=2.0;
lw=0.4;
lg=0.1;
wg=0.1;
t =0.03;
nx=20;
ny=2;
nz=5;

Point(1)  = {0.0         ,0.0, 0.0         ,0.1};
Point(2)  = {+lc*1.0/3.0 ,0.0, 0.0         ,0.1};
Point(3)  = {(lc-wg)/2.0 ,0.0, 0.0         ,0.1};
Point(4)  = {(lc+wg)/2.0 ,0.0, 0.0         ,0.1};
Point(5)  = {+lc*2.0/3.0 ,0.0, 0.0         ,0.1};
Point(6)  = {+lc         ,0.0, 0.0         ,0.1};

Point(7)  = {0.0         ,0.0, lw-lg       ,0.1};
Point(8)  = {+lc*1.0/3.0 ,0.0, lw-lg       ,0.1};
Point(9)  = {(lc-wg)/2.0 ,0.0, lw-lg       ,0.1};
Point(10) = {(lc+wg)/2.0 ,0.0, lw-lg       ,0.1};
Point(11) = {+lc*2.0/3.0 ,0.0, lw-lg       ,0.1};
Point(12) = {+lc         ,0.0, lw-lg       ,0.1};

Point(13) = {0.0         ,0.0, lw          ,0.1};
Point(14) = {+lc*1.0/3.0 ,0.0, lw          ,0.1};
Point(15) = {(lc-wg)/2.0 ,0.0, lw          ,0.1};
Point(16) = {(lc+wg)/2.0 ,0.0, lw          ,0.1};
Point(17) = {+lc*2.0/3.0 ,0.0, lw          ,0.1};
Point(18) = {+lc         ,0.0, lw          ,0.1};


// z= 0
Line(1)  = {1,2};
Line(2)  = {2,3};
Line(3)  = {3,4};
Line(4)  = {4,5};
Line(5)  = {5,6};

// z= lw-lg
Line(6)  = {7,8};
Line(7)  = {8,9};
Line(8)  = {9,10};
Line(9)  = {10,11};
Line(10) = {11,12};

// z= lw
Line(11)  = {13,14};
Line(12)  = {14,15};
//Line(13)  = {15,16};
Line(14)  = {16,17};
Line(15)  = {17,18};

Line(16) = {1,7};
Line(17) = {7,13};

Line(18) = {2,8};
Line(19) = {8,14};

Line(20) = {3,9};
Line(21) = {9,15};

Line(22) = {4,10};
Line(23) = {10,16};

Line(24) = {5,11};
Line(25) = {11,17};

Line(26) = {6,12};
Line(27) = {12,18};

Transfinite Line {16, 18, 20, 22, 24, 26, 2, 4, 7, 9, 12, 14} = 9 Using Progression 1;
Transfinite Line {1, 6, 11, 5, 10, 15} = 4 Using Progression 1;
Transfinite Line {17, 19, 21, 23, 25, 27} = 2 Using Progression 1;
Transfinite Line {3, 8} = 4 Using Progression 1;

Line Loop(28) = {16, 6, -18, -1};
Plane Surface(29) = {28};
Line Loop(30) = {24, 10, -26, -5};
Plane Surface(31) = {30};
Line Loop(32) = {6, 19, -11, -17};
Plane Surface(33) = {32};
Line Loop(34) = {25, 15, -27, -10};
Plane Surface(35) = {34};
Transfinite Surface {29};
Transfinite Surface {33};
Transfinite Surface {35};
Transfinite Surface {31};

Line Loop(36) = {19, 12, -21, -7};
Plane Surface(37) = {36};
Line Loop(38) = {23, 14, -25, -9};
Plane Surface(39) = {38};
Line Loop(40) = {20, 8, -22, -3};
Plane Surface(41) = {40};
Line Loop(42) = {22, 9, -24, -4};
Plane Surface(43) = {42};
Line Loop(44) = {18, 7, -20, -2};
Plane Surface(45) = {44};
Transfinite Surface {37};
Transfinite Surface {39};
Transfinite Surface {41};
Transfinite Surface {43};
Transfinite Surface {45};

Extrude {0, t, 0} {
  Surface{29, 33, 37, 45, 41, 39, 43, 31, 35};
  Layers{1};    
}

Physical Volume("MAT") = {1, 2, 3, 4, 6, 5, 7, 9, 8};
Physical Surface("x0") = {88, 54};
Physical Surface("x1") = {216, 238};
Physical Point("P100") = {6};
