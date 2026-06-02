#=====================================================================#
# Contains a derived type 'file' and implemented procedures
#=====================================================================#

import os

from .class_string import String

# Module-level counter to emulate Fortran's free-unit discovery.
_next_unit = [100]


class File:
    """Derived type for file."""

    #===================================================================#
    # Construct a file object with filename
    #
    # Fortran: interface file => create(filename, line_width)
    #===================================================================#
    def __init__(self, filename, line_width=None):
        # Set the file name
        self.filename = filename

        # Line width
        if line_width is not None:
            self.buffer_size = line_width
        else:
            self.buffer_size = 100

        # Use an available handle for opening
        self.file_unit = _next_unit[0]
        _next_unit[0] += 1

        # Python file object backing the unit (None until opened)
        self._handle = None

    #===================================================================#
    # Open the file
    #===================================================================#
    def open(self):
        file_exists = os.path.exists(self.filename)
        if file_exists is False:
            print(' file does not exist ', self.filename)
            raise SystemExit
        self._handle = open(self.file_unit_to_path())

    def file_unit_to_path(self):
        # Helper: the unit maps to this file's path.
        return self.filename

    #===================================================================#
    # Close the file
    #===================================================================#
    def close(self):
        if self._handle is not None:
            self._handle.close()
            self._handle = None

    #===================================================================#
    # return the file unit number
    #===================================================================#
    def get_unit(self):
        return self.file_unit

    #===================================================================#
    # Utility function for get number of lines in mesh file
    #===================================================================#
    def get_num_lines(self):
        nlines = 0
        self.open()
        for _ in self._handle:
            nlines = nlines + 1
        self.close()
        return nlines

    #=================================================================#
    # Read one line and return a string object
    #=================================================================#
    def read_line(self):
        # Fortran reads a fixed-width buffer then trims trailing blanks.
        buffer = self._handle.readline()
        line = String(buffer.rstrip())
        return line

    #=================================================================#
    # Read all lines and return a string array
    #=================================================================#
    def read_lines(self):
        # Get number of lines in file to allocate space
        num_lines = self.get_num_lines()
        lines = [None] * num_lines

        # Loop through each line and read and store into lines
        self.open()
        for iline in range(num_lines):
            lines[iline] = self.read_line()
        self.close()
        return lines
