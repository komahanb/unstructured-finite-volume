#!/usr/bin/env python3
"""
generate the project meshes with the gmsh python api and write msh 4.1.

    python meshgen/generate.py <name> <outpath>

<name> is one of the meshes the configs and tests reference. the physical
group names match exactly what the solver looks up by name:

  box      front/back/top/bottom/left/right + fluid    (structured hex)
  square   BoundaryLeft/Right/Top/Bottom    + surface  (structured quad, 2d)
  sphere   boundary                         + fluid    (tets)
  cylinder "cylindrical wall"/inlet/outlet  + fluid    (tets)

geometry placement is honored so the analytic oracles hold: front/back are the
y=0/y=1 faces (1d conduction -> mean 2.5), BoundaryLeft/Right are the x=0/x=1
edges (-> u = x). box and square are structured (orthogonal), which is what
makes those linear fields exact.
"""

import sys
import gmsh

EPS = 1.0e-6


def name_entities(dim, classifier):
    """group dimension-`dim` entities into physical groups by the name the
    classifier returns for each entity's center of mass."""
    groups = {}
    for (d, tag) in gmsh.model.getEntities(dim):
        label = classifier(gmsh.model.occ.getCenterOfMass(d, tag))
        if label is not None:
            groups.setdefault(label, []).append(tag)

    for label, tags in groups.items():
        g = gmsh.model.addPhysicalGroup(dim, tags)
        gmsh.model.setPhysicalName(dim, g, label)


def whole_group(dim, label):
    """one physical group covering every dimension-`dim` entity."""
    tags = [t for (_, t) in gmsh.model.getEntities(dim)]
    g = gmsh.model.addPhysicalGroup(dim, tags)
    gmsh.model.setPhysicalName(dim, g, label)


def structured(dim, n):
    """turn the current model into a structured (recombined) grid of n cells
    per edge - hexes in 3d, quads in 2d."""
    for (_, t) in gmsh.model.getEntities(1):
        gmsh.model.mesh.setTransfiniteCurve(t, n + 1)

    for (_, t) in gmsh.model.getEntities(2):
        gmsh.model.mesh.setTransfiniteSurface(t)
        gmsh.model.mesh.setRecombine(2, t)

    if dim == 3:
        for (_, t) in gmsh.model.getEntities(3):
            gmsh.model.mesh.setTransfiniteVolume(t)


def box(n):
    """unit cube, structured hex, n cells per edge."""
    gmsh.model.add("box")
    gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
    gmsh.model.occ.synchronize()

    structured(3, n)

    def face(c):
        x, y, z = c
        if abs(x)     < EPS: return "left"
        if abs(x - 1) < EPS: return "right"
        if abs(y)     < EPS: return "front"
        if abs(y - 1) < EPS: return "back"
        if abs(z)     < EPS: return "bottom"
        if abs(z - 1) < EPS: return "top"
        return None

    name_entities(2, face)
    whole_group(3, "fluid")
    gmsh.model.mesh.generate(3)


def square(n):
    """unit square, structured quad, n cells per edge (2d)."""
    gmsh.model.add("square")
    gmsh.model.occ.addRectangle(0, 0, 0, 1, 1)
    gmsh.model.occ.synchronize()

    structured(2, n)

    def edge(c):
        x, y, _ = c
        if abs(x)     < EPS: return "BoundaryLeft"
        if abs(x - 1) < EPS: return "BoundaryRight"
        if abs(y)     < EPS: return "BoundaryBottom"
        if abs(y - 1) < EPS: return "BoundaryTop"
        return None

    name_entities(1, edge)
    whole_group(2, "surface")
    gmsh.model.mesh.generate(2)


def sphere(h):
    """unit ball, unstructured tets; whole surface named 'boundary'."""
    gmsh.model.add("sphere")
    gmsh.model.occ.addSphere(0, 0, 0, 1)
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)

    whole_group(2, "boundary")
    whole_group(3, "fluid")
    gmsh.model.mesh.generate(3)


def cylinder(h):
    """cylinder along z, unstructured tets; caps inlet/outlet, lateral wall."""
    height = 2.0
    gmsh.model.add("cylinder")
    gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, height, 0.5)
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)

    def surf(c):
        _, _, z = c
        if abs(z)          < EPS: return "inlet"
        if abs(z - height) < EPS: return "outlet"
        return "cylindrical wall"

    name_entities(2, surf)
    whole_group(3, "fluid")
    gmsh.model.mesh.generate(3)


RECIPES = {
    "box-3":           lambda: box(12),
    "box-36":          lambda: box(3),
    "box-fine":        lambda: box(24),
    "sphere":          lambda: sphere(0.25),
    "cylinder-coarse": lambda: cylinder(0.30),
    "square-10":       lambda: square(10),
    "square-20":       lambda: square(20),
    "square-40":       lambda: square(40),
    "square-80":       lambda: square(80),
    "square-160":      lambda: square(160),
}


def main():
    if len(sys.argv) != 3 or sys.argv[1] not in RECIPES:
        sys.exit("usage: generate.py <%s> <outpath>" % " | ".join(RECIPES))

    name, out = sys.argv[1], sys.argv[2]

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    RECIPES[name]()
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.write(out)
    gmsh.finalize()
    print("wrote %s (%s)" % (out, name))


if __name__ == "__main__":
    main()
