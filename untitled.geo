//+
SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, 1, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Surface(1) = {1};
//+
Curve Loop(3) = {1};
//+
Surface(2) = {3};
//+
Curve Loop(5) = {1};
//+
Surface(3) = {5};
