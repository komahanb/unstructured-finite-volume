#!/usr/bin/env python3
"""The periodic box [0, L]^3 of n x n x n hexahedral cells in the MSH 4.1
format: the nodes, the boundary quadrilaterals on the six faces of the
box with their physical names, the hexahedra, and the three periodic
links - the top the translate of the bottom, the back of the front,
the right of the left - each listing every node of its face. gmsh's
transfinite volume does not accept a face copied from its periodic
master, so the box is written here instead of meshed by gmsh.

    python3 box_3d.py [n] [L] [file]     defaults 8, 2 pi, box_3d.msh
"""
import math, sys

n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
L = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0 * math.pi
name = sys.argv[3] if len(sys.argv) > 3 else 'box_3d.msh'
m = n + 1

def node(i, j, k):
    return 1 + i + m * (j + m * k)

out = []
out.append('$MeshFormat\n4.1 0 8\n$EndMeshFormat')
faces = ['bottom', 'top', 'front', 'back', 'left', 'right']
out.append('$PhysicalNames\n7')
for s, f in enumerate(faces, 1):
    out.append('2 %d "%s"' % (s, f))
out.append('3 7 "box"\n$EndPhysicalNames')
# entities: no points or curves, six surfaces, one volume
out.append('$Entities\n0 0 6 1')
for s in range(1, 7):
    out.append('%d 0 0 0 %g %g %g 1 %d 0' % (s, L, L, L, s))
out.append('1 0 0 0 %g %g %g 1 7 6 1 2 3 4 5 6\n$EndEntities' % (L, L, L))
# nodes, one block on the volume
out.append('$Nodes\n1 %d 1 %d\n3 1 0 %d' % (m ** 3, m ** 3, m ** 3))
for k in range(m):
    for j in range(m):
        for i in range(m):
            out.append(str(node(i, j, k)))
for k in range(m):
    for j in range(m):
        for i in range(m):
            out.append('%.16g %.16g %.16g' % (L * i / n, L * j / n, L * k / n))
out.append('$EndNodes')
# elements: quads per face, hexahedra in the volume
quads = {
    'bottom': [(node(i, j, 0), node(i + 1, j, 0), node(i + 1, j + 1, 0), node(i, j + 1, 0)) for j in range(n) for i in range(n)],
    'top':    [(node(i, j, n), node(i + 1, j, n), node(i + 1, j + 1, n), node(i, j + 1, n)) for j in range(n) for i in range(n)],
    'front':  [(node(i, 0, k), node(i + 1, 0, k), node(i + 1, 0, k + 1), node(i, 0, k + 1)) for k in range(n) for i in range(n)],
    'back':   [(node(i, n, k), node(i + 1, n, k), node(i + 1, n, k + 1), node(i, n, k + 1)) for k in range(n) for i in range(n)],
    'left':   [(node(0, j, k), node(0, j + 1, k), node(0, j + 1, k + 1), node(0, j, k + 1)) for k in range(n) for j in range(n)],
    'right':  [(node(n, j, k), node(n, j + 1, k), node(n, j + 1, k + 1), node(n, j, k + 1)) for k in range(n) for j in range(n)],
}
num_elements = 6 * n * n + n ** 3
out.append('$Elements\n7 %d 1 %d' % (num_elements, num_elements))
tag = 0
for s, f in enumerate(faces, 1):
    out.append('2 %d 3 %d' % (s, n * n))
    for q in quads[f]:
        tag += 1
        out.append('%d %d %d %d %d' % ((tag,) + q))
out.append('3 1 5 %d' % n ** 3)
for k in range(n):
    for j in range(n):
        for i in range(n):
            tag += 1
            h = (node(i, j, k), node(i + 1, j, k), node(i + 1, j + 1, k), node(i, j + 1, k),
                 node(i, j, k + 1), node(i + 1, j, k + 1), node(i + 1, j + 1, k + 1), node(i, j + 1, k + 1))
            out.append('%d %d %d %d %d %d %d %d %d' % ((tag,) + h))
out.append('$EndElements')
# periodic links: every node of the image face paired with its master's
def affine(tx, ty, tz):
    return '16 1 0 0 %.16g 0 1 0 %.16g 0 0 1 %.16g 0 0 0 1' % (tx, ty, tz)
links = [
    ('2 2 1', affine(0, 0, L), [(node(i, j, n), node(i, j, 0)) for j in range(m) for i in range(m)]),
    ('2 4 3', affine(0, L, 0), [(node(i, n, k), node(i, 0, k)) for k in range(m) for i in range(m)]),
    ('2 6 5', affine(L, 0, 0), [(node(n, j, k), node(0, j, k)) for k in range(m) for j in range(m)]),
]
out.append('$Periodic\n3')
for head, aff, pairs in links:
    out.append(head)
    out.append(aff)
    out.append(str(len(pairs)))
    for a, b in pairs:
        out.append('%d %d' % (a, b))
out.append('$EndPeriodic')
open(name, 'w').write('\n'.join(out) + '\n')
print('%s: %d nodes, %d hexahedra, %d boundary quadrilaterals, 3 periodic links' % (name, m ** 3, n ** 3, 6 * n * n))
