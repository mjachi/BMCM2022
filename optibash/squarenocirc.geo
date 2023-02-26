// Gmsh project created on Sun Nov 13 06:18:42 2022
//+
Point(1) = {1, 0, 0, 1.0};
//+
Point(2) = {1, 1, 0, 1.0};
//+
Point(3) = {0, 1, 0, 1.0};
//+
Point(4) = {0, 0, 0, 1.0};
//+
Line(1) = {3, 2};
//+
Line(2) = {2, 1};
//+
Line(3) = {4, 1};
//+
Line(4) = {4, 3};
//+
Curve Loop(1) = {1, 2, -3, 4};
//+
Surface(1) = {1};
//+
Point(5) = {0.5, 1, 0, 1.0};
//+
Circle(5) = {2, 5, 3};
//+
Curve Loop(2) = {5, 1};
//+
Plane Surface(2) = {2};
//+
Point(6) = {0.5, 0, 0, 1.0};
//+
Circle(6) = {4, 6, 1};
//+
Curve Loop(3) = {6, -3};
//+
Plane Surface(3) = {3};
