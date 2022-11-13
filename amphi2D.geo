// Gmsh project created on Sat Nov 12 16:15:29 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0.25, 0, 0, 1.0};
//+
Point(3) = {-0.25, 0, 0, 1.0};
//+
Circle(1) = {3, 1, 2};
//+
Point(4) = {0, 0.25, 0, 1.0};
//+
Circle(1) = {3, 1, 4};
//+
Circle(2) = {2, 1, 4};
//+
Point(5) = {-0.25, -0.75, 0, 1.0};
//+
Point(6) = {0.25, -0.75, 0, 1.0};
//+
Line(3) = {3, 5};
//+
Line(4) = {5, 6};
//+
Line(5) = {6, 2};
//+
Curve Loop(1) = {5, 2, -1, 3, 4};
//+
Surface(1) = {1};
