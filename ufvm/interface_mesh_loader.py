#=====================================================================#
# Abstract type for mesh loading
#=====================================================================#

from abc import ABC, abstractmethod


class MeshLoader(ABC):
    """Abstract type for mesh loading."""

    #-------------------------------------------------------------------#
    # Type bound procedure that returns all information needed for
    # mesh creation.
    #
    # Returns a tuple of all the out-arguments of the Fortran
    # get_mesh_data_interface, in declaration order:
    #
    #   (num_vertices, vertex_numbers, vertex_tags, vertices,
    #    num_edges, edge_numbers, edge_tags, edge_vertices, num_edge_vertices,
    #    num_faces, face_numbers, face_tags, face_vertices, num_face_vertices,
    #    num_cells, cell_numbers, cell_tags, cell_vertices, num_cell_vertices,
    #    num_tags, tag_numbers, tag_info)
    #-------------------------------------------------------------------#
    @abstractmethod
    def get_mesh_data(self):
        ...
