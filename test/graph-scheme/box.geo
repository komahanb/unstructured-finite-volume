// The periodic box of side 2 pi: 16 x 16 quadrilateral cells; the top
// edge is the translate of the bottom and the right of the left.
L = 2*Pi;
n = 16;
Point(1) = {0, 0, 0};
Point(2) = {L, 0, 0};
Point(3) = {L, L, 0};
Point(4) = {0, L, 0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {4, 3};
Line(4) = {1, 4};
Curve Loop(1) = {1, 2, -3, -4};
Plane Surface(1) = {1};
Transfinite Curve {1, 2, 3, 4} = n + 1;
Transfinite Surface {1};
Recombine Surface {1};
Periodic Curve {3} = {1} Translate {0, L, 0};
Periodic Curve {2} = {4} Translate {L, 0, 0};
Physical Curve("bottom") = {1};
Physical Curve("right")  = {2};
Physical Curve("top")    = {3};
Physical Curve("left")   = {4};
Physical Surface("box")  = {1};
Mesh.MshFileVersion = 4.1;
