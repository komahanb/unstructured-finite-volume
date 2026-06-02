#=====================================================================#
# Unstructured mesh handler.
#
# Author: Komahan Boopathy (komahan@gatech.edu)
#
# Index convention: connectivity arrays store 1-based IDs as values
# (cell/face/vertex numbers, as the Fortran does); array positions are
# 0-based. An ID indexes another array with (id - 1).
#=====================================================================#

import sys

import numpy as np

from .module_mesh_utils import (
    reverse_map, get_cell_faces, cross_product, find, distance,
)


class Mesh:
    """Mesh datatype. A collection of vertices, cells and faces."""

    #================================================================!
    # Constructor for mesh object using mesh loader
    #================================================================!
    def __init__(self, loader):
        self.max_print = 20
        self.initialized = False

        # Get the fundamental information needed
        (self.num_vertices, self.vertex_numbers, self.vertex_tags, self.vertices,
         self.num_edges, self.edge_numbers, self.edge_tags, self.edge_vertices, self.num_edge_vertices,
         self.num_faces, self.face_numbers, self.face_tags, self.face_vertices, self.num_face_vertices,
         self.num_cells, self.cell_numbers, self.cell_tags, self.cell_vertices, self.num_cell_vertices,
         self.num_tags, self.tag_numbers, self.tag_info) = loader.get_mesh_data()

        # Derived connectivity (set during initialize)
        self.vertex_cells = None
        self.num_vertex_cells = None
        self.vertex_faces = None
        self.num_vertex_faces = None
        self.vertex_edges = None
        self.num_vertex_edges = None
        self.num_cell_faces = None
        self.cell_faces = None
        self.num_face_cells = None
        self.face_cells = None

        # Sanity check (numbering continuous, may not start from one)
        if (self.num_vertices > 0 and
                self.vertex_numbers.max() - self.vertex_numbers.min() + 1 != self.num_vertices):
            raise SystemExit
        if (self.num_edges > 0 and
                self.edge_numbers.max() - self.edge_numbers.min() + 1 != self.num_edges):
            raise SystemExit
        if (self.num_faces > 0 and
                self.face_numbers.max() - self.face_numbers.min() + 1 != self.num_faces):
            raise SystemExit
        if (self.num_cells > 0 and
                self.cell_numbers.max() - self.cell_numbers.min() + 1 != self.num_cells):
            raise SystemExit

        # Perform initialization tasks and store the resulting flag
        self.initialized = self.initialize()
        if self.initialized is False:
            sys.stderr.write("Mesh.Construct: failed\n")
            raise SystemExit

    #===================================================================!
    # Build derived connectivity and geometry information.
    #===================================================================!
    def initialize(self):
        #---------------------------------------------------------------!
        # Find VertexCell conn. by inverting CellVertex conn.
        #---------------------------------------------------------------!
        print(" Inverting CellVertex Map...")
        self.vertex_cells, self.num_vertex_cells = reverse_map(
            self.cell_vertices, self.num_cell_vertices)

        if self.vertex_cells is not None:
            print("Vertex to cell info for", min(self.max_print, self.num_vertices),
                  " vertices out of ", self.num_vertices)
            for ivertex in range(1, min(self.max_print, self.num_vertices) + 1):
                print('vertex [', self.vertex_numbers[ivertex - 1], ']',
                      'num cells [', self.num_vertex_cells[ivertex - 1], ']',
                      'cells [', self.vertex_cells[0:self.num_vertex_cells[ivertex - 1], ivertex - 1], ']')
            # Sanity check
            if self.num_vertex_cells.min() < 1:
                sys.stderr.write('Error: There are vertices not mapped to a cell\n')
                raise SystemExit
        else:
            print("Vertex to cell info not computed")

        #---------------------------------------------------------------!
        # Find VertexFace conn. by inverting FaceVertex conn.
        #---------------------------------------------------------------!
        print(" Inverting FaceVertex Map...")
        self.vertex_faces, self.num_vertex_faces = reverse_map(
            self.face_vertices, self.num_face_vertices)

        if self.vertex_faces is not None:
            print("Vertex to face info for", min(self.max_print, self.num_vertices),
                  " vertices out of ", self.num_vertices)
            for ivertex in range(1, min(self.max_print, self.num_vertices) + 1):
                print('vertex [', self.vertex_numbers[ivertex - 1], ']',
                      'num_vertex_faces [', self.num_vertex_faces[ivertex - 1], ']',
                      'faces [', self.vertex_faces[0:self.num_vertex_faces[ivertex - 1], ivertex - 1], ']')
            # Sanity check
            if self.num_vertex_faces.min() < 1:
                sys.stderr.write('Error: There are vertices not mapped to a face\n')
                raise SystemExit
        else:
            print("Vertex to face info not computed")

        #---------------------------------------------------------------!
        # Find VertexEdge conn. by inverting EdgeVertex conn.
        #---------------------------------------------------------------!
        print(" Inverting EdgeVertex Map...")
        self.vertex_edges, self.num_vertex_edges = reverse_map(
            self.edge_vertices, self.num_edge_vertices)

        if self.vertex_edges is not None:
            print("Vertex to edge info for", min(self.max_print, self.num_vertices),
                  " vertices out of ", self.num_vertices)
            for ivertex in range(1, min(self.max_print, self.num_vertices) + 1):
                print('vertex [', self.vertex_numbers[ivertex - 1], ']',
                      'num_vertex_edges [', self.num_vertex_edges[ivertex - 1], ']',
                      'edges [', self.vertex_edges[0:self.num_vertex_edges[ivertex - 1], ivertex - 1], ']')
            if self.num_vertex_edges.min() < 1:
                sys.stderr.write('Error: There are vertices not mapped to a edge\n')
                raise SystemExit
        else:
            print("Vertex to edge info not computed")

        #---------------------------------------------------------------!
        # Find Cell Face conn. by combining two maps
        #---------------------------------------------------------------!
        print(" Combining CellVertex with VertexFace to get CellFace Map...")
        self.cell_faces, self.num_cell_faces = get_cell_faces(
            self.cell_vertices, self.vertex_faces, self.num_vertex_faces)

        if self.cell_faces is not None:
            print("Cell to face info for", min(self.max_print, self.num_cells),
                  " cells out of ", self.num_cells)
            for icell in range(1, min(self.max_print, self.num_cells) + 1):
                print('cell [', icell, ']',
                      'nfaces [', self.num_cell_faces[icell - 1], ']',
                      'faces [', self.cell_faces[0:self.num_cell_faces[icell - 1], icell - 1], ']')
        else:
            print("Vertex to edge info not computed")

        #---------------------------------------------------------------!
        # Find Face Cell conn. by inverting Cell Face conn.
        #---------------------------------------------------------------!
        self.face_cells, self.num_face_cells = reverse_map(
            self.cell_faces, self.num_cell_faces)

        for iface in range(1, min(self.max_print, self.num_faces) + 1):
            print('face [', self.face_numbers[iface - 1], ']',
                  'cells [', self.face_cells[0:self.num_face_cells[iface - 1], iface - 1], ']')

        if self.num_face_cells.min() < 1:
            sys.stderr.write('Error: There are faces not mapped to a cell\n')

        # tag_cells_vertices
        # Tag everything as domain
        self.cell_tags[:] = 0
        self.vertex_tags[:] = 0

        # Loop through boundary faces and tag cells and vertices
        for iface in range(1, self.num_faces + 1):
            # will work only at this point the domain is marked 0
            if self.face_tags[iface - 1] > 0:
                nfc = self.num_face_cells[iface - 1]
                # Copy face tags into corresponding cells
                cell_ids = self.face_cells[0:nfc, iface - 1]
                self.cell_tags[cell_ids - 1] = self.face_tags[iface - 1]
                # Copy face tags into corresponding vertices
                vert_ids = self.face_vertices[0:nfc, iface - 1]
                self.vertex_tags[vert_ids - 1] = self.face_tags[iface - 1]

        # Set the domain tags to maxval
        maxtagnum = self.tag_numbers.max()
        self.face_tags[self.face_tags == 0] = maxtagnum
        self.edge_tags[self.edge_tags == 0] = maxtagnum
        self.cell_tags[self.cell_tags == 0] = maxtagnum
        self.vertex_tags[self.vertex_tags == 0] = maxtagnum

        # number_tagged_faces
        self.num_tagged_faces = np.zeros(self.num_tags, dtype=int)

        mintag = self.tag_numbers.min()
        maxtag = self.tag_numbers.max()

        # Count the number of tagged faces of each kind
        for iface in range(1, self.num_faces + 1):
            for itag in range(mintag, maxtag + 1):
                ftag = self.face_tags[iface - 1]
                if itag == ftag:
                    self.num_tagged_faces[itag - 1] = 1 + self.num_tagged_faces[itag - 1]
                    break

        # Use the counted number to allocate and recount
        self.tagged_face_face = np.zeros((int(self.num_tagged_faces.max()), self.num_tags), dtype=int)
        self.num_tagged_faces[:] = 0
        for iface in range(1, self.num_faces + 1):
            for itag in range(mintag, maxtag + 1):
                ftag = self.face_tags[iface - 1]
                if itag == ftag:
                    self.num_tagged_faces[itag - 1] = 1 + self.num_tagged_faces[itag - 1]
                    self.tagged_face_face[self.num_tagged_faces[itag - 1] - 1, itag - 1] = iface
                    break

        #---------------------------------------------------------------!
        # Evaluate all geometric quantities needed for FVM assembly
        #---------------------------------------------------------------!
        print(' Calculating mesh geometry information')

        self.evaluate_cell_centers()
        self.evaluate_face_centers_areas()
        self.evaluate_face_tangents_normals()
        self.evaluate_cell_volumes()

        self.evaluate_centroidal_vector()
        self.evaluate_face_deltas()
        self.evaluate_face_weight()
        self.evaluate_vertex_weight()

        # Signal that all tasks are complete
        return True

    def evaluate_vertex_weight(self):
        print(' Evaluating face weights for interpolation from cells to vertex')

        cells = np.zeros(int(self.num_vertex_cells.max()), dtype=int)

        self.vertex_cell_weights = np.zeros(
            (int(self.num_vertex_cells.max()), self.num_vertices), dtype=np.float64)

        for ivertex in range(1, self.num_vertices + 1):
            nvc = self.num_vertex_cells[ivertex - 1]

            # actual cells numbers
            cells[0:nvc] = self.vertex_cells[0:nvc, ivertex - 1]

            total = 0.0

            for icell in range(1, nvc + 1):
                dcell = distance(self.cell_centers[:, cells[icell - 1] - 1],
                                 self.vertices[:, ivertex - 1])
                self.vertex_cell_weights[icell - 1, ivertex - 1] = 1.0 / dcell
                total = total + self.vertex_cell_weights[icell - 1, ivertex - 1]

            self.vertex_cell_weights[0:nvc, ivertex - 1] = \
                self.vertex_cell_weights[0:nvc, ivertex - 1] / total

        for ivertex in range(1, min(self.max_print, self.num_vertices) + 1):
            nvc = self.num_vertex_cells[ivertex - 1]
            print("vertex [", self.vertex_numbers[ivertex - 1], ']',
                  "weights [", self.vertex_cell_weights[0:nvc, ivertex - 1], ']')

    def evaluate_face_weight(self):
        print(' Evaluating face weights for interpolation from cells to face')
        self.face_cell_weights = np.zeros((2, self.num_faces), dtype=np.float64)
        for iface in range(1, self.num_faces + 1):
            # first cell is found for all faces
            cellindex1 = self.face_cells[0, iface - 1]
            xcellcenter1 = self.cell_centers[:, cellindex1 - 1]
            xfacecenter = self.face_centers[:, iface - 1]
            d1 = distance(xcellcenter1, xfacecenter)
            dinv1 = 1.0 / d1

            # Extract the second cell if this is not a boundary face/hole
            if self.num_face_cells[iface - 1] != 1:
                cellindex2 = self.face_cells[1, iface - 1]
                xcellcenter2 = self.cell_centers[:, cellindex2 - 1]
                d2 = distance(xcellcenter2, xfacecenter)
                dinv2 = 1.0 / d2
            else:
                dinv2 = 0.0

            weight = dinv1 / (dinv1 + dinv2)

            self.face_cell_weights[0:2, iface - 1] = [weight, 1.0 - weight]

        for iface in range(1, min(self.max_print, self.num_faces) + 1):
            print("face [", iface, "] ",
                  "weight [", self.face_cell_weights[0:2, iface - 1], "] ")

    def evaluate_face_deltas(self):
        print(" Evaluating face deltas")
        self.face_deltas = np.zeros(self.num_faces, dtype=np.float64)

        for gface in range(1, self.num_faces + 1):
            # First cell belonging to the face
            gcell = self.face_cells[0, gface - 1]

            # Face number in local numbering (find returns 1-based index)
            lface = find(self.cell_faces[:, gcell - 1], gface)

            # Index into normal array
            fn = self.cell_face_normals[:, lface - 1, gcell - 1]

            # Take absolute value of dot product
            self.face_deltas[gface - 1] = abs(np.dot(self.lvec[0:3, gface - 1], fn))

            print("face [", gface, "] ",
                  "delta [", self.face_deltas[gface - 1], "] ",
                  "skewness t.l [",
                  np.dot(self.lvec[0:3, gface - 1], self.cell_face_tangents[:, lface - 1, gcell - 1]), "] ",
                  "orthogonality t.n [",
                  np.dot(self.cell_face_tangents[:, lface - 1, gcell - 1],
                         self.cell_face_normals[:, lface - 1, gcell - 1]), "] ",
                  self.cell_face_normals[:, lface - 1, gcell - 1])

        # Check for negative volumes
        if abs(self.face_deltas.min()) < 1.0e-10:
            print('collinear faces/bad cell?')
            raise SystemExit

    def evaluate_centroidal_vector(self):
        print(" Evaluating centroidal vector...")

        self.lvec = np.zeros((3, self.num_faces), dtype=np.float64)

        for iface in range(1, self.num_faces + 1):
            cells = np.zeros(2, dtype=int)
            nfc = self.num_face_cells[iface - 1]
            cells[0:nfc] = self.face_cells[0:nfc, iface - 1]

            # boundary face if not the highest tag
            if self.face_tags[iface - 1] < self.tag_numbers.max():
                # Boundary faces .or. iface is in bfaces
                self.lvec[:, iface - 1] = \
                    self.face_centers[:, iface - 1] - self.cell_centers[:, cells[0] - 1]
            else:
                # Interior face; subtract neighbouring cell centers
                if self.num_face_cells[iface - 1] == 1:
                    # Accounts for internal holes
                    self.lvec[:, iface - 1] = \
                        self.face_centers[:, iface - 1] - self.cell_centers[:, cells[0] - 1]
                else:
                    self.lvec[:, iface - 1] = \
                        self.cell_centers[:, cells[1] - 1] - self.cell_centers[:, cells[0] - 1]

    def evaluate_cell_volumes(self):
        print(" Evaluating cell volumes...")

        self.cell_volumes = np.zeros(self.num_cells, dtype=np.float64)

        # V = \sum_f nx_f * xmid_f * A_f
        for lcell in range(1, self.num_cells + 1):
            self.cell_volumes[lcell - 1] = 0.0
            for lface in range(1, self.num_cell_faces[lcell - 1] + 1):
                # Global face index
                gface = self.cell_faces[lface - 1, lcell - 1]
                xmid = self.face_centers[0, gface - 1]
                nx = self.cell_face_normals[0, lface - 1, lcell - 1]
                area = self.face_areas[gface - 1]
                self.cell_volumes[lcell - 1] = self.cell_volumes[lcell - 1] + nx * xmid * area

        # Check for negative volumes
        if self.cell_volumes.min() < 0:
            print('negative volume encountered')
            raise SystemExit

    def evaluate_face_centers_areas(self):
        print(' Evaluating face centers and areas')

        self.face_areas = np.zeros(self.num_faces, dtype=np.float64)
        self.face_centers = np.zeros((3, self.num_faces), dtype=np.float64)

        for iface in range(1, self.num_faces + 1):
            facenodes = self.face_vertices[:, iface - 1]   # 1-based vertex ids

            # Compute the coordinates of face centers (this face has 2 edges)
            self.face_centers[0:3, iface - 1] = \
                self.vertices[0:3, facenodes - 1].sum(axis=1) / 2.0

            # Compute face areas
            v1 = self.vertices[:, facenodes[0] - 1]
            v2 = self.vertices[:, facenodes[1] - 1]
            self.face_areas[iface - 1] = distance(v1, v2)

        # Check for zero areas
        if abs(self.face_areas.min()) < 10.0 * np.finfo(np.float64).tiny:
            print('same points/bad face?')
            raise SystemExit

    def evaluate_cell_centers(self):
        print(' Evaluating cell centers')

        self.cell_centers = np.zeros((3, self.num_cells), dtype=np.float64)

        for icell in range(1, self.num_cells + 1):
            ncv = self.num_cell_vertices[icell - 1]
            vids = self.cell_vertices[0:ncv, icell - 1]   # 1-based vertex ids
            self.cell_centers[:, icell - 1] = \
                self.vertices[:, vids - 1].sum(axis=1) / float(ncv)

    def evaluate_face_tangents_normals(self):
        print(' Evaluating face tangents normals')

        maxcf = int(self.num_cell_faces.max())
        self.cell_face_normals = np.zeros((3, maxcf, self.num_cells), dtype=np.float64)
        self.cell_face_tangents = np.zeros((3, maxcf, self.num_cells), dtype=np.float64)

        # loop cells
        for icell in range(1, self.num_cells + 1):
            icv = self.cell_vertices[:, icell - 1]   # 1-based vertex ids (may include 0)

            # loop faces of each cell
            for iface in range(1, self.num_cell_faces[icell - 1] + 1):
                if iface == self.num_cell_faces[icell - 1]:
                    ifv1 = icv[iface - 1]
                    ifv2 = icv[0]
                else:
                    ifv1 = icv[iface - 1]
                    ifv2 = icv[iface + 1 - 1]

                # find the face vertex in cell order
                gface = self.cell_faces[iface - 1, icell - 1]

                t = self.vertices[:, ifv2 - 1] - self.vertices[:, ifv1 - 1]
                t = t / np.linalg.norm(t)

                # By anticlockwise convention
                n = np.zeros(3, dtype=np.float64)
                n[0] = t[1]
                n[1] = -t[0]
                n[2] = 0

                # Sanity check if the normal is facing out of the face
                tcn = np.zeros(3, dtype=np.float64)
                cross_product(n, t, tcn)
                if abs(tcn[2] - 1.0) > 1.0e-10:   # tangent x normal should give +k
                    print('face', gface, 'of cell', icell, 'has inward/non-unit normal', tcn[2])
                    raise SystemExit

                self.cell_face_normals[:, iface - 1, icell - 1] = n
                self.cell_face_tangents[:, iface - 1, icell - 1] = t

    #===================================================================!
    # Counts and returns the number of elements with supplied number of
    # nodes
    #===================================================================!
    def get_num_elems(self, num_elem_nodes):
        return int(np.count_nonzero(self.num_cell_vertices == num_elem_nodes))

    #===================================================================!
    # Prints the mesh topology and geometry information
    #===================================================================!
    def to_string(self):
        print(' Number of vertices :', self.num_vertices)
        print(' Number of cells    :', self.num_cells)
        print(' Number of faces    :', self.num_faces)

        for itag in range(1, self.num_tags + 1):
            print("tag number [", self.tag_numbers[itag - 1], "] ",
                  "info [", self.tag_info[itag - 1].str, "] ")

        if self.num_vertices > 0:
            print("Vertex info for ", min(self.max_print, self.num_vertices),
                  ' vertices out of ', self.num_vertices)
            print("number tag x y z")
            for ivertex in range(1, min(self.max_print, self.num_vertices) + 1):
                print(self.vertex_numbers[ivertex - 1],
                      self.vertex_tags[ivertex - 1],
                      self.vertices[:, ivertex - 1])

        if self.num_cells > 0:
            print("Cell info for ", min(self.max_print, self.num_cells),
                  ' cells out of ', self.num_cells)
            print("cno ctag ncv iverts")
            for icell in range(1, min(self.max_print, self.num_cells) + 1):
                print(self.cell_numbers[icell - 1],
                      self.cell_tags[icell - 1],
                      self.num_cell_vertices[icell - 1],
                      self.cell_vertices[0:self.num_cell_vertices[icell - 1], icell - 1])

        if self.num_faces > 0:
            print("Face info for ", min(self.max_print, self.num_faces),
                  ' faces out of ', self.num_faces)
            print("fno ftag nfv iverts")
            for iface in range(1, min(self.max_print, self.num_faces) + 1):
                print(self.face_numbers[iface - 1],
                      self.face_tags[iface - 1],
                      self.num_face_vertices[iface - 1],
                      self.face_vertices[0:self.num_face_vertices[iface - 1], iface - 1])

        if self.num_edges > 0:
            print("Edge info for ", min(self.max_print, self.num_edges),
                  ' edges out of ', self.num_edges)
            print("eno etag nev iverts")
            for iedge in range(1, min(self.max_print, self.num_edges) + 1):
                print(self.edge_numbers[iedge - 1],
                      self.edge_tags[iedge - 1],
                      self.num_edge_vertices[iedge - 1],
                      self.edge_vertices[0:self.num_edge_vertices[iedge - 1], iedge - 1])

        if self.initialized is True:
            print("Cell Geo. Data [index] [center] [volume] ")
            for icell in range(1, min(self.max_print, self.num_cells) + 1):
                print("local number [", self.cell_numbers[icell - 1], "] ",
                      "center [", self.cell_centers[:, icell - 1], "] ",
                      "volume [", self.cell_volumes[icell - 1], "] ")

            print("Face Data [index] [center] [area] ")
            for iface in range(1, min(self.max_print, self.num_faces) + 1):
                print("num [", iface, "] ",
                      "center [", self.face_centers[:, iface - 1], "] ",
                      "area [", self.face_areas[iface - 1], "] ",
                      "lvec [", self.lvec[:, iface - 1], "] ",
                      "delta [", self.face_deltas[iface - 1], "] ")
