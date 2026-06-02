#=====================================================================#
# Concrete mesh_loader for GMSH format meshes.
#
# Parsing-heavy translation: Fortran 1-based token/line indices are kept
# as 1-based loop counters and converted to 0-based at every subscript
# (token(k) -> tokens[k-1], lines(a:b) -> lines[a-1:b]).
#=====================================================================#

import numpy as np

from .interface_mesh_loader import MeshLoader
from .class_file import File
from .class_string import String
from .class_set import Set
from .class_list import List


class GmshLoader(MeshLoader):

    #-------------------------------------------------------------------#
    # Interface to construct a mesh_loader for GMSH
    #-------------------------------------------------------------------#
    def __init__(self, filename):
        self.file = File(filename)   # mesh file

    #====================================================================#
    # Supply all information needed to create a mesh object
    #====================================================================#
    def get_mesh_data(self):
        # Load the mesh into memory
        print("Loading mesh file :" + self.file.filename)
        lines = self.file.read_lines()

        print("Identifying tags...")
        (idx_start_mesh, idx_end_mesh,
         idx_start_physical_names, idx_end_physical_names,
         idx_start_nodes, idx_end_nodes,
         idx_start_elements, idx_end_elements) = find_tags(lines)

        print("mesh           : ", idx_start_mesh, idx_end_mesh)
        print("physical names : ", idx_start_physical_names, idx_end_physical_names)
        print("nodes          : ", idx_start_nodes, idx_end_nodes)
        print("elements       : ", idx_start_elements, idx_end_elements)

        print("Reading mesh information...")

        # process_mesh_version
        mlines = lines[idx_start_mesh:idx_start_mesh + 1]   # lines(idx_start_mesh+1 : idx_start_mesh+1)
        num_tokens, tokens = mlines[0].tokenize(" ")
        if int(np.floor(tokens[0].asreal())) > 2:
            print("unsupported version of gmsh version ", tokens[0].str)
            raise SystemExit

        print("Reading physical tags...")

        # process_tags
        tag_lines = lines[idx_start_physical_names:idx_end_physical_names - 1]  # (+1 : -1)

        # Set the intent(out) variable for number of tags present
        num_tags = tag_lines[0].asinteger()

        # Allocate space for other two return variables
        tag_info = [None] * num_tags
        tag_numbers = np.zeros(num_tags, dtype=int)

        for iline in range(1, num_tags + 1):
            # Tokenize based on delimited space
            num_tokens, tokens = tag_lines[(iline + 1) - 1].tokenize(" ")

            # Second tag is the tag number
            tag_numbers[iline - 1] = tokens[2 - 1].asinteger()

            # Remove quotes on third tag
            length = len(tokens[3 - 1].str)
            tag_info[iline - 1] = String(tokens[3 - 1].str[2 - 1:length - 1])

        print("Reading vertices...")

        # process_nodes
        vlines = lines[idx_start_nodes + 1:idx_end_nodes - 1]   # (+2 : -1)

        # Set the number of vertices
        num_vertices = len(vlines)

        vertices = np.zeros((3, num_vertices), dtype=np.float64)
        vertex_numbers = np.zeros(num_vertices, dtype=int)
        vertex_tags = np.zeros(num_vertices, dtype=int)

        # Parse lines and store vertices
        for ivertex in range(1, num_vertices + 1):
            # Get the numbers of tokens and tokens
            num_tokens, tokens = vlines[ivertex - 1].tokenize(" ")

            # First token is the vertex number
            vertex_numbers[ivertex - 1] = ivertex

            # Second, third and fourth tokens are the coordinates
            vertices[:, ivertex - 1] = [t.asreal() for t in tokens[2 - 1:4]]

        print("Reading elements...")

        # elements
        elines = lines[idx_start_elements + 1:idx_end_elements - 1]   # (+2 : -1)

        num_lines = len(elines)

        # Zero counters for elements
        num_edges = 0
        num_faces = 0
        num_cells = 0

        # Count the number of cells present in elements
        for iline in range(1, num_lines + 1):
            num_tokens, tokens = elines[iline - 1].tokenize(" ")
            if tokens[2 - 1].asinteger() == 2:
                # Triangular element
                num_cells = num_cells + 1
            elif tokens[2 - 1].asinteger() == 3:
                # Quadrilateral element
                num_cells = num_cells + 1

        # Allocate space based on counted num of cells
        cell_numbers = np.zeros(num_cells, dtype=int)
        num_cell_vertices = np.zeros(num_cells, dtype=int)
        cell_vertices = np.zeros((4, num_cells), dtype=int)
        cell_tags = np.zeros(num_cells, dtype=int)

        # Create space for processing face information
        set_face_vertices = Set(2, 4 * num_cells)
        set_face_numbers = Set(1, 4 * num_cells)
        list_face_tags = List(1, 4 * num_cells)

        face_idx = 0
        cell_idx = 0

        for iline in range(1, num_lines + 1):
            num_tokens, tokens = elines[iline - 1].tokenize(" ")

            etype = tokens[2 - 1].asinteger()

            # Line element
            if etype == 1:
                # Carry out processing of physically tagged faces
                added = set_face_vertices.insert(
                    [t.asinteger() for t in tokens[6 - 1:7]])
                if added is True:
                    face_idx = face_idx + 1
                    set_face_numbers.add_entry([face_idx])
                    list_face_tags.add_entry([tokens[4 - 1].asinteger()])

            # Triangular element
            elif etype == 2:
                cell_idx = cell_idx + 1
                cell_numbers[cell_idx - 1] = cell_idx
                cell_tags[cell_idx - 1] = tokens[4 - 1].asinteger()
                cell_vertices[1 - 1:3, cell_idx - 1] = [t.asinteger() for t in tokens[6 - 1:8]]
                num_cell_vertices[cell_idx - 1] = 3

            # Quadrilateral element
            elif etype == 3:
                cell_idx = cell_idx + 1
                cell_numbers[cell_idx - 1] = cell_idx
                cell_tags[cell_idx - 1] = tokens[4 - 1].asinteger()
                cell_vertices[1 - 1:4, cell_idx - 1] = [t.asinteger() for t in tokens[6 - 1:9]]
                num_cell_vertices[cell_idx - 1] = 4

            # Node
            elif etype == 15:
                pass

            else:
                print('unknown element number', tokens[2 - 1].asinteger())
                raise SystemExit

        print("Finding faces...")

        # process_faces: Make ordered pair of vertices as faces in 2D
        for icell in range(1, num_cells + 1):
            for iverpair in range(1, num_cell_vertices[icell - 1] + 1):
                # Figure out a vertex pair that makes a face
                if iverpair == num_cell_vertices[icell - 1]:
                    idx = [iverpair, 1]
                else:
                    idx = [iverpair, iverpair + 1]

                # Add ordered pair of integers into set
                pair = cell_vertices[[idx[0] - 1, idx[1] - 1], icell - 1]
                added = set_face_vertices.insert(pair)
                if added is True:
                    face_idx = face_idx + 1
                    set_face_numbers.add_entry([face_idx])
                    list_face_tags.add_entry([0])

        num_faces = face_idx

        num_face_vertices = np.full(num_faces, 2, dtype=int)

        face_vertices = set_face_vertices.get_entries()

        # Set face numbers
        lface_numbers = set_face_numbers.get_entries()
        face_numbers = np.array(lface_numbers[0, :])

        # Set face tags
        lface_tags = list_face_tags.get_entries()
        face_tags = np.array(lface_tags[0, :])

        print("Processing edges...")

        # Edges are not in 2D
        edge_numbers = np.zeros(num_edges, dtype=int)
        edge_vertices = np.zeros((4, num_edges), dtype=int)
        num_edge_vertices = np.zeros(num_edges, dtype=int)
        edge_tags = np.zeros(num_edges, dtype=int)

        if num_faces > 4 * num_cells:
            raise SystemExit

        return (num_vertices, vertex_numbers, vertex_tags, vertices,
                num_edges, edge_numbers, edge_tags, edge_vertices, num_edge_vertices,
                num_faces, face_numbers, face_tags, face_vertices, num_face_vertices,
                num_cells, cell_numbers, cell_tags, cell_vertices, num_cell_vertices,
                num_tags, tag_numbers, tag_info)


#====================================================================#
# Locate the $MeshFormat / $PhysicalNames / $Nodes / $Elements section
# markers. Returns 1-based line indices, matching the Fortran.
#====================================================================#
def find_tags(lines):
    BEGIN_MESH = "$MeshFormat"
    END_MESH = "$EndMeshFormat"
    BEGIN_PHYSICAL_NAMES = "$PhysicalNames"
    END_PHYSICAL_NAMES = "$EndPhysicalNames"
    BEGIN_NODES = "$Nodes"
    END_NODES = "$EndNodes"
    BEGIN_ELEMENTS = "$Elements"
    END_ELEMENTS = "$EndElements"

    idx_start_mesh = idx_end_mesh = 0
    idx_start_physical_names = idx_end_physical_names = 0
    idx_start_nodes = idx_end_nodes = 0
    idx_start_elements = idx_end_elements = 0

    num_lines = len(lines)
    for iline in range(1, num_lines + 1):
        s = lines[iline - 1].str

        # Find mesh start and end
        if s.find(BEGIN_MESH) == 0:
            idx_start_mesh = iline
        if s.find(END_MESH) == 0:
            idx_end_mesh = iline

        # Find physical_names start and end
        if s.find(BEGIN_PHYSICAL_NAMES) == 0:
            idx_start_physical_names = iline
        if s.find(END_PHYSICAL_NAMES) == 0:
            idx_end_physical_names = iline

        # Find nodes start and end
        if s.find(BEGIN_NODES) == 0:
            idx_start_nodes = iline
        if s.find(END_NODES) == 0:
            idx_end_nodes = iline

        # Find elements start and end
        if s.find(BEGIN_ELEMENTS) == 0:
            idx_start_elements = iline
        if s.find(END_ELEMENTS) == 0:
            idx_end_elements = iline

    return (idx_start_mesh, idx_end_mesh,
            idx_start_physical_names, idx_end_physical_names,
            idx_start_nodes, idx_end_nodes,
            idx_start_elements, idx_end_elements)
