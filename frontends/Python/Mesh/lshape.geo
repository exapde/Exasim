element_size = 0.2;
Point(1) = {-1.0, -1.0, 0.0, element_size };
Point(2) = {0.0, -1.0, 0.0, element_size };
Point(3) = {0.0, 0.0, 0.0, element_size/10};
Point(4) = {1.0, 0.0, 0.0, element_size };
Point(5) = {1.0, 1.0, 0.0, element_size };
Point(6) = {-1.0, 1.0, 0.0, element_size };

Line(7) = {1,2};
Line(8) = {2,3};
Line(9) = {3,4};
Line(10) = {4,5};
Line(11) = {5,6};
Line(12) = {6,1};

Line Loop(13) = {7,8,9,10,11,12};

Plane Surface(14) = {13};

