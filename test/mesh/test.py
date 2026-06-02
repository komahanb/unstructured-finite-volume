#=====================================================================#
# Test loading of mesh and mesh pre-processing
#
# Port of test/mesh/test.f90 (program test_mesh).
#=====================================================================#

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from ufvm import Mesh, GmshLoader, File, String


def test_gmsh_loader(filename):
    # Create a mesh loader for mesh file
    gmsh_loader_obj = GmshLoader(filename)
    mesh_obj = Mesh(gmsh_loader_obj)
    mesh_obj.to_string()


def main():
    #===================================================================#
    # Test the functionalities of Class GMSH_LOADER
    #===================================================================#
    files = [None] * 5
    files[1 - 1] = String('../rectangle.msh')
    files[2 - 1] = String('../square-10.msh')
    files[3 - 1] = String('../triangle.msh')
    files[4 - 1] = String('../frontal.msh')
    files[5 - 1] = String('../delaunay.msh')

    for ifile in range(1, len(files) + 1):
        print("Testing GMSH Loader with file ", files[ifile - 1].str)
        test_gmsh_loader(files[ifile - 1].str)

    #===================================================================#
    # Test the functionalities of Class String and Class File
    #===================================================================#
    filename = '../rectangle.msh'

    # Create a file object
    file_obj = File(filename)

    # Do a line by line read and print contents
    num_lines = file_obj.get_num_lines()
    lines = [None] * num_lines

    file_obj.open()
    for iline in range(1, num_lines + 1):
        lines[iline - 1] = file_obj.read_line()
    file_obj.close()
    for line in lines:
        line.print()

    # Read everything and print content
    lines = file_obj.read_lines()
    for line in lines:
        line.print()


if __name__ == '__main__':
    main()
