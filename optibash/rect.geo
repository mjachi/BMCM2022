// Gmsh project created on Sun Nov 13 13:49:02 2022
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 2, 0, 1.0};
//+
Point(4) = {0, 2, 0, 1.0};
//+
Line(1) = {3, 2};
//+
Line(2) = {2, 1};
//+
Line(3) = {1, 4};
//+
Line(4) = {4, 3};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Point(5) = {1, 2, 0.2, 1.0};
//+
Point(6) = {0.5, 2, 0, 1.0};
//+
Circle(5) = {3, 6, 4};
//+
Recursive Delete {
  Curve{4}; 
}
//+
Recursive Delete {
  Curve{4}; 
}
//+
Curve Loop(2) = {1, 2, 3, -5};
//+
Plane Surface(2) = {2};
//+
Point(7) = {0.5, 0, 0, 1.0};
//+
Circle(6) = {1, 7, 2};
//+
Curve Loop(3) = {6, -1, 5, -3};
//+
Plane Surface(3) = {3};
