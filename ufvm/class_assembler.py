#=====================================================================#
# Class responsible for matrix, right hand side assembly and boundary
# conditions.
#
# Index convention: cell/face state vectors (x, Ax, b, ss, phi*) are
# indexed by the 1-based cell/face/vertex id minus one.
#=====================================================================#

import numpy as np


class Assembler:

    #===================================================================!
    # Constructor for physics
    #
    # Fortran: interface assembler => construct(grid)
    #===================================================================!
    def __init__(self, grid):
        print("constructing assembler")

        # set symmetry to .true. for structured grid
        self.symmetry = True

        # Set mesh
        self.grid = grid
        self.grid.to_string()

        # Non symmetric jacobian
        self.symmetry = True

        # Number of state variables (in FVM, the number of cells)
        self.num_state_vars = self.grid.num_cells

        # Allocate the flux vector
        self.phi = np.zeros(self.num_state_vars, dtype=np.float64)

        # Matrix filters
        self.DIAGONAL = 0
        self.LOWER_TRIANGLE = -1
        self.UPPER_TRIANGLE = 1

    #===================================================================!
    # Assemble and return the full jacobian matrix
    #===================================================================!
    def get_jacobian(self, filter=None):
        A = np.zeros((self.num_state_vars, self.num_state_vars), dtype=np.float64)
        ex = np.zeros(self.num_state_vars, dtype=np.float64)

        if filter is not None:
            # Assemble only a part of the matrix
            for icol in range(1, self.num_state_vars + 1):
                ex[icol - 1] = 1.0
                self.get_jacobian_vector_product(A[:, icol - 1], ex, filter)
                ex[icol - 1] = 0.0
        else:
            # Assemble Full Matrix A = L + D + U
            for icol in range(1, self.num_state_vars + 1):
                ex[icol - 1] = 1.0
                self.get_jacobian_vector_product(A[:, icol - 1], ex)
                ex[icol - 1] = 0.0

        return A

    #===================================================================!
    # Assemble and return the full transpose jacobian matrix
    #===================================================================!
    def get_transpose_jacobian(self, filter=None):
        A = np.zeros((self.num_state_vars, self.num_state_vars), dtype=np.float64)
        ex = np.zeros(self.num_state_vars, dtype=np.float64)

        if filter is not None:
            for irow in range(1, self.num_state_vars + 1):
                ex[irow - 1] = 1.0
                self.get_jacobian_vector_product(A[irow - 1, :], ex, filter)
                ex[irow - 1] = 0.0
        else:
            for irow in range(1, self.num_state_vars + 1):
                ex[irow - 1] = 1.0
                self.get_jacobian_vector_product(A[irow - 1, :], ex)
                ex[irow - 1] = 0.0

        return A

    def get_jacobian_vector_product(self, Ax, x, filter=None):
        # laplace_normal
        grid = self.grid
        highest_tag = grid.tag_numbers.max()

        # Loop cells
        for icell in range(1, grid.num_cells + 1):
            ncf = grid.num_cell_faces[icell - 1]
            faces = grid.cell_faces[0:ncf, icell - 1]   # 1-based face ids

            # Loop faces
            Ax[icell - 1] = 0.0

            for iface in range(1, ncf + 1):
                gf = faces[iface - 1]
                ftag = grid.face_tags[gf - 1]
                fdelta = grid.face_deltas[gf - 1]
                farea = grid.face_areas[gf - 1]
                nfcells = grid.num_face_cells[gf - 1]

                # Add contribution from internal faces
                if ftag == highest_tag:

                    # Neighbour cell index
                    fcells = grid.face_cells[0:nfcells, gf - 1]

                    # Neighbour is the one with a different cell index
                    if fcells[0] == icell:
                        ncell = fcells[1]
                    else:
                        ncell = fcells[0]

                    if filter is not None:
                        # Assemble part of matrix
                        if filter == self.UPPER_TRIANGLE:
                            if ncell > icell:
                                Ax[icell - 1] = Ax[icell - 1] + farea * (x[ncell - 1]) / fdelta
                        elif filter == self.LOWER_TRIANGLE:
                            if ncell < icell:
                                Ax[icell - 1] = Ax[icell - 1] + farea * (x[ncell - 1]) / fdelta
                        elif filter == self.DIAGONAL:
                            Ax[icell - 1] = Ax[icell - 1] + farea * (-x[icell - 1]) / fdelta
                    else:
                        # Assemble Full of matrix
                        Ax[icell - 1] = Ax[icell - 1] + farea * (x[ncell - 1] - x[icell - 1]) / fdelta

                else:

                    # Diagonal self contributions (boundary faces)
                    if filter is not None:
                        if filter == self.DIAGONAL:
                            Ax[icell - 1] = Ax[icell - 1] + farea * (0.0 - x[icell - 1]) / fdelta
                    else:
                        Ax[icell - 1] = Ax[icell - 1] + farea * (0.0 - x[icell - 1]) / fdelta

    #===================================================================!
    # Compute vertex values by interpolating cell center values
    #===================================================================!
    def evaluate_vertex_flux(self, phic):
        grid = self.grid
        phiv = np.zeros(grid.num_vertices, dtype=np.float64)
        for ivertex in range(1, grid.num_vertices + 1):
            nvc = grid.num_vertex_cells[ivertex - 1]
            w = grid.vertex_cell_weights[0:nvc, ivertex - 1]
            icells = grid.vertex_cells[0:nvc, ivertex - 1]   # 1-based cell ids
            phiv[ivertex - 1] = np.dot(phic[icells - 1], w)
        return phiv

    #===================================================================!
    # Compute face center values by interpolating cell center values
    #===================================================================!
    def evaluate_face_flux(self, phic):
        grid = self.grid
        phif = np.zeros(grid.num_faces, dtype=np.float64)
        for iface in range(1, grid.num_faces + 1):
            nfc = grid.num_face_cells[iface - 1]
            w = grid.face_cell_weights[0:nfc, iface - 1]
            icells = grid.face_cells[0:nfc, iface - 1]   # 1-based cell ids
            phif[iface - 1] = np.dot(phic[icells - 1], w)
        return phif

    def get_transpose_jacobian_vector_product(self, Ax, x):
        Ax[:] = 0.0
        raise SystemExit('not implemented')

    #===================================================================!
    # Evaluate internal skew source based on the current cell states.
    # ss is written in place.
    #===================================================================!
    def get_skew_source(self, ss, phic):
        grid = self.grid

        # Evaluate nodal values of phi
        phiv = self.evaluate_vertex_flux(phic)

        # Make space for skew source terms
        ss[:] = 0

        highest_tag = grid.tag_numbers.max()

        for icell in range(1, grid.num_cells + 1):
            ncf = grid.num_cell_faces[icell - 1]
            faces = grid.cell_faces[0:ncf, icell - 1]

            for iface in range(1, ncf + 1):
                # Global face number
                gface = faces[iface - 1]

                # Ignore boundary faces from skew source evaluation
                if grid.face_tags[gface - 1] == highest_tag:

                    # Compute tangent.dot.lvector/delta
                    scale = np.dot(
                        grid.lvec[0:3, gface - 1],
                        grid.cell_face_tangents[0:3, iface - 1, icell - 1]
                    ) / grid.face_deltas[gface - 1]

                    # Get the vertices associated with this face
                    nfv = grid.num_face_vertices[gface - 1]
                    fvertices = grid.face_vertices[0:nfv, gface - 1]   # 1-based vertex ids

                    ss[icell - 1] = ss[icell - 1] + scale * (
                        phiv[fvertices[1] - 1] - phiv[fvertices[0] - 1])

    def get_source(self, b):
        grid = self.grid

        # Homogenous dirichlet boundary conditions
        phib = 0.0

        # add_boundary_terms
        highest_tag = grid.tag_numbers.max()
        for icell in range(1, grid.num_cells + 1):
            ncf = grid.num_cell_faces[icell - 1]
            faces = grid.cell_faces[0:ncf, icell - 1]

            b[icell - 1] = 0.0

            for iface in range(1, ncf + 1):
                gf = faces[iface - 1]
                ftag = grid.face_tags[gf - 1]
                fdelta = grid.face_deltas[gf - 1]
                farea = grid.face_areas[gf - 1]

                # Add contribution from internal faces
                if ftag != highest_tag:   # homogenous dirichlet
                    # Boundary faces (minus as we moved it to rhs)
                    b[icell - 1] = b[icell - 1] + farea * (-phib) / fdelta

        # cell_source
        for icell in range(1, grid.num_cells + 1):
            x = grid.cell_centers[:, icell - 1]
            cell_volume = grid.cell_volumes[icell - 1]
            b[icell - 1] = b[icell - 1] + evaluate_source(x) * cell_volume

    #===================================================================!
    # Write solution to file
    #===================================================================!
    def write_solution(self, filename, phic):
        grid = self.grid

        # Open resource
        path = filename.strip()
        try:
            f = open(path, 'w')
        except OSError:
            print("  >> Opening file ", path, " failed")
            return

        # Compute vertex values by interpolating cell center values
        phiv = self.evaluate_vertex_flux(phic)

        # Write header
        f.write(' TITLE = "FVM-Laplace"\n')
        f.write(' VARIABLES = "x" "y"  "T"\n')

        # Write Triangles/Quads (works only for homogeneous elements)
        if grid.num_cell_vertices.max() == 4:
            f.write(' ZONE T="Temperature", N= %d , E= %d , DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n'
                    % (grid.num_vertices, grid.num_cells))
        else:
            f.write(' ZONE T="Temperature", N= %d , E= %d , DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n'
                    % (grid.num_vertices, grid.num_cells))

        # Write vertices
        for i in range(1, grid.num_vertices + 1):
            f.write("   %s %s\n" % (
                " ".join(repr(v) for v in grid.vertices[0:2, i - 1]),
                repr(float(phiv[i - 1]))))

        # Write cell connectivities
        for i in range(1, grid.num_cells + 1):
            ncv = grid.num_cell_vertices[i - 1]
            f.write("   %s\n" % " ".join(str(int(v)) for v in grid.cell_vertices[0:ncv, i - 1]))

        # Close resource
        f.close()

    #===================================================================!
    # Create a state vector and sets values if a scalar is supplied
    #===================================================================!
    def create_vector(self, x=None, scalar=None):
        if x is not None:
            raise SystemExit("vector already allocated")
        x = np.zeros(self.num_state_vars, dtype=np.float64)
        if scalar is not None:
            x[:] = scalar
        return x


def evaluate_source(x):
    return 1.0 + 0.0 * x[0]   # -1.0 ! sin(x(1)) + cos(x(2))
