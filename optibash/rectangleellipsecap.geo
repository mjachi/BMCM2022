// Gmsh project created on Sun Nov 13 15:24:21 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 2, 0, 1.0};
//+
Point(4) = {0, 2, 0, 1.0};
//+
Point(5) = {0.5, 2, 0, 1.0};
//+
Point(6) = {0.5, 2.25, 0, 1.0};
//+
Ellipse(1) = {4, 5, 6, 3};
//+
Ellipse(1) = {6, 5, 4, 3};
//+
Ellipse(2) = {6, 5, 3, 4};
//+
Line(3) = {4, 1};
//+
Line(4) = {1, 2};
//+
Line(5) = {2, 3};
//+
Curve Loop(1) = {5, -1, 2, 3, 4};
//+
Plane Surface(1) = {1};
