// Gmsh project created on Sun Nov 13 12:33:51 2022
//+
Point(1) = {1, 0, -0, 1.0};
//+
Point(2) = {0, 0, -0, 1.0};
//+
Point(3) = {0.5, 0.8660254037844386, -0, 1.0};
//+
Point(4) = {0.5, 0, -0, 1.0};
//+
Line(1) = {3, 2};
//+
Line(2) = {3, 1};
//+
Circle(3) = {2, 4, 1};
//+
Curve Loop(1) = {2, -3, -1};
//+
Plane Surface(1) = {1};
