lc = 0.05;
R = 0.254;
L = 1.27;

Point(1) = {L,0.0,0.0,lc};
Point(2) = {L,R,0,lc};
Point(3) = {L,0,R,lc};
Point(4) = {L,-R,0,lc};
Point(5) = {L,0,-R,lc};
Point(6) = {0.0,0.0,0.0,lc/4};
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
Line(5) = {6,2};
Line(6) = {6,3};
Line(7) = {6,4};
Line(8) = {6,5};
Line Loop(9) = {6,-1,-5};
Surface(10) = {9};
Line Loop(11) = {7,-2,-6};
Surface(12) = {11};
Line Loop(13) = {8,-3,-7};
Surface(14) = {13};
Line Loop(15) = {5,-4,-8};
Surface(16) = {15};
Line Loop(17) = {4,3,2,1};
Surface(18) = {17};
Surface Loop(19) = {10,12,14,16,18};

lc2 = 1;
w = 10;
Point(10) = {-w,-w,-w,lc2};
Point(11) = {w,-w,-w,lc2};
Point(12) = {w,w,-w,lc2};
Point(13) = {-w,w,-w,lc2};
Point(14) = {-w,-w,w,lc2};
Point(15) = {w,-w,w,lc2};
Point(16) = {w,w,w,lc2};
Point(17) = {-w,w,w,lc2};
Line(31) = {14,15};
Line(32) = {15,11};
Line(33) = {11,10};
Line(34) = {10,14};
Line(35) = {14,17};
Line(36) = {17,13};
Line(37) = {13,10};
Line(38) = {13,12};
Line(39) = {12,16};
Line(40) = {16,17};
Line(41) = {16,15};
Line(42) = {11,12};
Line Loop(43) = {38,39,40,36};
Plane Surface(44) = {43};
Line Loop(45) = {36,37,34,35};
Plane Surface(46) = {45};
Line Loop(47) = {35,-40,41,-31};
Plane Surface(48) = {47};
Line Loop(49) = {32,42,39,41};
Plane Surface(50) = {49};
Line Loop(51) = {42,-38,37,-33};
Plane Surface(52) = {51};
Line Loop(53) = {33,34,31,32};
Plane Surface(54) = {53};
Surface Loop(55) = {44,52,-50,54,-46,48};
Volume(56) = {55,19};










