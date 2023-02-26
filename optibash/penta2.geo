// Gmsh project created on Sun Nov 13 14:14:14 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {1, 0, 0, 1.0};
//+
Point(2) = {0.3090169943749475, 0.9510565162951536, 0, 1.0};
//+
Point(3) = {-0.8090169943749475,  0.5877852522924731, 0, 1.0};
//+
Point(4) = {-0.8090169943749477, -0.5877852522924729, 0, 1.0};
//+
Point(5) = {0.3090169943749473, -0.9510565162951538, 0, 1.0};//+
Line(1) = {2, 1};
//+
Line(2) = {1, 5};
//+
Line(3) = {5, 4};
//+
Line(4) = {2, 3};
//+
Line(5) = {3, 4};
//+
Curve Loop(1) = {2, 3, -5, -4, 1};
//+
Plane Surface(1) = {1};