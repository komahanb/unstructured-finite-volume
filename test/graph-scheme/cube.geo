// The unit cube of n x n x n cells: hexahedra, or tetrahedra where
// hexahedra = 0. gmsh -3 cube.geo -setnumber n 4 -o cube_4.msh
If (!Exists(n))
  n = 4;
EndIf
If (!Exists(hexahedra))
  hexahedra = 1;
EndIf
Point(1) = {0, 0, 0};
Point(2) = {1, 0, 0};
Line(1) = {1, 2};
Transfinite Curve {1} = n + 1;
If (hexahedra)
  face[] = Extrude {0, 1, 0} { Curve{1}; Layers{n}; Recombine; };
  body[] = Extrude {0, 0, 1} { Surface{face[1]}; Layers{n}; Recombine; };
Else
  face[] = Extrude {0, 1, 0} { Curve{1}; Layers{n}; };
  body[] = Extrude {0, 0, 1} { Surface{face[1]}; Layers{n}; };
EndIf
Physical Surface("boundary") = {face[1], body[0], body[2], body[3], body[4], body[5]};
Physical Volume("cube") = {body[1]};
Mesh.MshFileVersion = 4.1;
