// Gmsh project created on Sun Nov 13 12:19:20 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 0, 0, 1.0};
//+
Point(4) = {0.5, 0.8660254037844386, 0, 1.0};
//+
Line(1) = {4, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 4};
//+
Curve Loop(1) = {3, 1, 2};
//+
Plane Surface(1) = {1};
