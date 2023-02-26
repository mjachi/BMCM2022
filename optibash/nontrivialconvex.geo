// Gmsh project created on Sat Nov 12 23:04:26 2022
SetFactory("OpenCASCADE");
//+
Ellipse(1) = {0.0, 0.0, 0, 1, 0.5, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Plane Surface(1) = {1};
