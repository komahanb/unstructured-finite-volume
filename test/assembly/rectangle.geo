/*********************************************************************
 *
 *  Gmsh tutorial 1
 *
 *  Variables, elementary entities (points, curves, surfaces), physical
 *  entities (points, curves, surfaces)
 *
 *********************************************************************/

// The simplest construction in Gmsh's scripting language is the
// `affectation'. The following command defines a new variable `lc':

lc = 0.1;

// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is defined by a list of four numbers:
// three coordinates (X, Y and Z), and a characteristic length (lc) that sets
// the target element size at the point:

Point(1) = {0, 0, 0, lc};

// The distribution of the mesh element sizes is then obtained by interpolation
// of these characteristic lengths throughout the geometry. Another method to
// specify characteristic lengths is to use general mesh size Fields (see
// `t10.geo'). A particular case is the use of a background mesh (see `t7.geo').

// We can then define some additional points as well as our first curve.  Curves
// are Gmsh's second type of elementery entities, and, amongst curves, straight
// lines are the simplest. A straight line is defined by a list of point
// numbers. In the commands below, for example, the line 1 starts at point 1 and
// ends at point 2:

Point(2) = {1, 0,  0, lc} ;
Point(3) = {1, 1, 0, lc} ;
Point(4) = {0, 1, 0, lc} ;

Line(1) = {1,2} ;
Line(2) = {2,3} ;
Line(3) = {3,4} ;
Line(4) = {4,1} ;

// The third elementary entity is the surface. In order to define a simple
// rectangular surface from the four curves defined above, a curve loop has first
// to be defined. A curve loop is a list of connected curves, a sign being
// associated with each curve (depending on the orientation of the curve):

Curve Loop(1) = {1,2,3,4} ;

// We can then define the surface as a list of curve loops (only one here, since
// there are no holes--see `t4.geo'):

Plane Surface(1) = {1} ;

// At this level, Gmsh knows everything to display the rectangular surface 6 and
// to mesh it. An optional step is needed if we want to group elementary
// geometrical entities into more meaningful groups, e.g. to define some
// mathematical ("domain", "boundary"), functional ("left wing", "fuselage") or
// material ("steel", "carbon") properties.
//
// Such groups are called "Physical Groups" in Gmsh. By default, if physical
// groups are defined, Gmsh will export in output files only those elements that
// belong to at least one physical group. (To force Gmsh to save all elements,
// whether they belong to physical groups or not, set "Mesh.SaveAll=1;", or
// specify "-save_all" on the command line.)
//
// Here we define a physical curve that groups the left, bottom and right lines
// in a single group (with prescribed tag 5); and a physical surface with name
// "My surface" (with an automatic tag) containg the geometrical surface 1:

Physical Curve("BoundaryBottom") = {1} ;
Physical Curve("BoundaryRight") = {2} ;
Physical Curve("BoundaryTop") = {3} ;
Physical Curve("BoundaryLeft") = {4} ;
Physical Surface("Domain") = {1} ;