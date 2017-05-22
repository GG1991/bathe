/*  CLWL-DCB SPECIMEN
 *  Autor: Guido Giuntoli
 */

le=0.457;
lc=0.003;
la=0.038;
lb=0.070;
t =0.0508;

Point(1)  = {0.0       ,0.0      ,0.0  ,0.1};
Point(2)  = {0.0       ,le       ,0.0  ,0.1};
Point(3)  = {le        ,le       ,0.0  ,0.1};
Point(4)  = {le        ,0.0      ,0.0  ,0.1};
Point(5)  = {le/2-lc/2 ,0.0      ,0.0  ,0.1};
Point(6)  = {le/2+lc/2 ,0.0      ,0.0  ,0.1};

Point(7)  = {le/2-lc/2 ,la       ,0.0  ,0.1};
Point(8)  = {le/2-lc/2 ,3*la     ,0.0  ,0.1};
Point(9)  = {le/2-lc/2 ,3*la+lb  ,0.0  ,0.1};

Point(10) = {le/2+lc/2 ,la       ,0.0  ,0.1};
Point(11) = {le/2+lc/2 ,3*la     ,0.0  ,0.1};
Point(12) = {le/2+lc/2 ,3*la+lb  ,0.0  ,0.1};

Point(13) = {le/2      ,2*la     ,0.0  ,0.1};
Point(14) = {le/2-la   ,2*la     ,0.0  ,0.1};
Point(15) = {le/2+la   ,2*la     ,0.0  ,0.1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 6};
Line(5) = {6, 10};
Line(6) = {11, 12};
Line(7) = {12, 9};
Line(8) = {9, 8};
Line(9) = {7, 5};
Line(10) = {5, 1};
Circle(11) = {10, 13, 15};
Circle(12) = {15, 13, 11};
Circle(13) = {8, 13, 14};
Circle(14) = {14, 13, 7};

Transfinite Line {7} = 2 Using Progression 1;
Transfinite Line {1, 2, 3} = 10 Using Progression 1;
Transfinite Line {10, 4} = 5 Using Progression 1;
Transfinite Line {8, 6} = 5 Using Progression 1;
Transfinite Line {9, 5} = 5 Using Progression 1;
Transfinite Line {14, 11, 13, 12} = 4 Using Progression 1;

Line Loop(15) = {1, 2, 3, 4, 5, 11, 12, 6, 7, 8, 13, 14, 9, 10};
Plane Surface(16) = {15};
//Recombine Surface{16};

vol[]=Extrude {0, 0, 0.05} {
  Surface{16};
  Layers{1};
  Recombine;
};

Physical Volume("MAT") = {1};

Physical Point("FIX_F") = {14};
Physical Point("MOV_F") = {15};
Physical Point("FIX_B") = {60};
Physical Point("MOV_B") = {38};
