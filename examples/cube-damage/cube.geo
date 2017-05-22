lc=0.5;

Point(1) = {-lc,-lc,-lc,0.1};
Point(2) = {lc,-lc,-lc,0.1};
Point(3) = {lc,lc,-lc,0.1};
Point(4) = {-lc, lc,-lc,0.1};

Line(1) = {2,3};
Line(2) = {3,4};
Line(3) = {4,1};
Line(4) = {1,2};

Transfinite Line {1, 3} = 1 Using Progression 1;
Transfinite Line {2, 4} = 2 Using Progression 1;
Line Loop(5) = {2,3,4,1};
Plane Surface(6) = {5};
Transfinite Surface {6};
Recombine Surface {6};

Physical Point("P000") = {1};
Physical Point("P100") = {2};
Physical Point("P110") = {3};
Physical Point("P010") = {4};
Physical Point("P001") = {10};
Physical Point("P101") = {14};
Physical Point("P111") = {5};
Physical Point("P011") = {6};

vol[]=Extrude {0, 0, 2*lc} {
  Surface{6};Layers{1};Recombine;
};
Physical Volume("MAT") = {vol[1]};
