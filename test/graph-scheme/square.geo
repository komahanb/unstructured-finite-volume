// The unit square of n x n cells: quadrangles, or triangles where
// quadrangles = 0. gmsh -2 square.geo -setnumber n 8 -o square_8.msh
If (!Exists(n))
  n = 8;
EndIf
If (!Exists(quadrangles))
  quadrangles = 1;
EndIf
Point(1) = {0, 0, 0};
Point(2) = {1, 0, 0};
Point(3) = {1, 1, 0};
Point(4) = {0, 1, 0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Curve {1, 2, 3, 4} = n + 1;
Transfinite Surface {1};
If (quadrangles)
  Recombine Surface {1};
EndIf
Physical Curve("boundary") = {1, 2, 3, 4};
Physical Surface("square") = {1};
Mesh.MshFileVersion = 4.1;
