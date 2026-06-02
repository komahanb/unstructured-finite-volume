#=====================================================================#
# module_mesh_utils - free functions for mesh connectivity / geometry
#
# Index convention for the port: connectivity arrays store 1-based IDs
# as *values* (as the Fortran does), while array *positions* are 0-based.
# An ID is converted to a position with (id - 1); a 0-based loop position
# is converted back to a stored ID with (position + 1).
#=====================================================================#

import numpy as np


#===================================================================#
# Cross product for area computations
#===================================================================#
def cross_product(a, b, pdt):
    # use skew form and generalize to n-dimensions?
    pdt[0] = a[1] * b[2] - a[2] * b[1]
    pdt[1] = a[2] * b[0] - a[0] * b[2]
    pdt[2] = a[0] * b[1] - a[1] * b[0]


#===================================================================#
# Geometric distance between two points
#===================================================================#
def distanceX(X):
    # X(:,:) -> [[x,y,z], [1:2]]
    return np.sqrt(np.sum((X[:, 0] - X[:, 1]) ** 2))


def distanceAB(x, y):
    return np.sqrt(np.sum((np.asarray(x) - np.asarray(y)) ** 2))


def distance(*args):
    """Generic interface `distance` dispatching on argument count."""
    if len(args) == 1:
        return distanceX(args[0])
    else:
        return distanceAB(args[0], args[1])


#===================================================================#
# Invert a map. The keys of 'map' are values of 'inverse'. The
# values of 'map' are the keys of 'inverse'
#===================================================================#
def reverse_map(map, num_map_vals):
    map = np.asarray(map)
    num_map_vals = np.asarray(num_map_vals)

    if num_map_vals.size == 0:
        return None, None

    # Forward mapping size
    nkeysin = map.shape[1]
    nvalsin = map.shape[0]

    # Nothing to do (probably empty map)!
    if nkeysin == 0:
        return None, None

    # Find the maximum size of the inverse map values
    nkeysout = int(map.max())   # dangerous
    num_inverse_vals = np.zeros(nkeysout, dtype=int)
    for i in range(nkeysin):
        for j in range(num_map_vals[i]):
            value = map[j, i]
            num_inverse_vals[value - 1] = num_inverse_vals[value - 1] + 1
    nvalsout = int(num_inverse_vals.max())

    # Allocate the inverse map based on determined sizes
    inverse = np.zeros((nvalsout, nkeysout), dtype=int)

    # Point into the next available slot for each inverse_key
    ptr = np.zeros(nkeysout, dtype=int)
    for i in range(nkeysin):
        for j in range(num_map_vals[i]):
            value = map[j, i]
            ptr[value - 1] = ptr[value - 1] + 1
            inverse[ptr[value - 1] - 1, value - 1] = i + 1   # store key as 1-based id

    return inverse, num_inverse_vals


#===================================================================#
# Find the intersection of two arrays (move elsewhere)?
# Writes the first common element into c (in place), as in Fortran.
#===================================================================#
def intersection(a, b, c):
    sizea = np.asarray(a).size
    sizeb = np.asarray(b).size

    ctr = 0
    for i in range(sizea):
        for j in range(sizeb):
            # Copy entry to new list if equal
            if a[i] == b[j]:
                ctr = ctr + 1
                c[ctr - 1] = a[i]
                return


#===================================================================#
# Find if a target value is present in the array (move else where)?
# Returns the 1-based index of the value (or -1), matching Fortran.
#===================================================================#
def find(array, target_value):
    nentries = np.asarray(array).shape[0]

    for i in range(nentries):
        if array[i] == target_value:
            return i + 1

    return -1


#===================================================================#
# Checks if the first argument is a subset of the second argument
#===================================================================#
def is_subset(small, big):
    small = np.asarray(small)
    big = np.asarray(big)

    if small.size > big.size:
        return False

    # Create local copy of arrays
    sub = np.array(small)
    set_ = np.array(big)

    # Sort two arrays
    isort(sub)
    isort(set_)
    lensub = sub.size

    # Check if all entries are equal upto the length of the smallest array
    for i in range(lensub):
        if not (set_ == sub[i]).any():
            return False
    return True


#===================================================================#
# Sort an integer array ! move elsewhere?
#===================================================================#
def isort(array):
    n = array.shape[0]
    for j in range(n):
        for k in range(j + 1, n):
            if array[j] > array[k]:
                temp = array[k]
                array[k] = array[j]
                array[j] = temp


#===================================================================#
# Forms the cell faces from a pair of vertices belonging to cell.
#===================================================================#
def get_cell_faces(cell_vertices, vertex_faces, num_vertex_faces):
    cell_vertices = np.asarray(cell_vertices)
    vertex_faces = np.asarray(vertex_faces)
    num_vertex_faces = np.asarray(num_vertex_faces)

    nvertices = cell_vertices.shape[0]
    nfaces = nvertices
    ncells = cell_vertices.shape[1]

    # find how many faces are there based on nodes
    num_cell_faces = np.zeros(ncells, dtype=int)
    for icell in range(ncells):
        ctr = 0
        for iface in range(nvertices):
            if cell_vertices[iface, icell] != 0:
                ctr = ctr + 1
        num_cell_faces[icell] = ctr

    # Cell to face cell_vertices
    cell_faces = np.zeros((int(num_cell_faces.max()), ncells), dtype=int)
    for icell in range(ncells):
        face_ptr = 0
        for iface in range(num_cell_faces[icell]):
            # Get the first two vertices
            if iface == num_cell_faces[icell] - 1:
                v1 = cell_vertices[iface, icell]
                v2 = cell_vertices[0, icell]
            else:
                v1 = cell_vertices[iface, icell]
                v2 = cell_vertices[iface + 1, icell]

            face_ptr = face_ptr + 1
            # v1, v2 are 1-based vertex ids -> index with (v - 1)
            intersection(
                vertex_faces[0:num_vertex_faces[v2 - 1], v2 - 1],
                vertex_faces[0:num_vertex_faces[v1 - 1], v1 - 1],
                cell_faces[face_ptr - 1:face_ptr, icell])

    return cell_faces, num_cell_faces


#===================================================================#
# Determine if the face is a boundary face based on how many
# neighbouring cells it has.
#===================================================================#
def get_boundary_faces(num_face_cells):
    num_face_cells = np.asarray(num_face_cells)
    nfaces = num_face_cells.shape[0]

    # Boundary faces are the faces corresponding to just one cell
    nbfaces = 0
    for iface in range(nfaces):
        if num_face_cells[iface] == 1:
            nbfaces = nbfaces + 1

    boundary_faces = np.zeros(nbfaces, dtype=int)
    ctr = 0
    for iface in range(nfaces):
        if num_face_cells[iface] == 1:
            ctr = ctr + 1
            boundary_faces[ctr - 1] = iface + 1   # store 1-based face id

    return boundary_faces
