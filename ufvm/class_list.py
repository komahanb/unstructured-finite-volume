#=====================================================================#
# Basic implementation of an unordered list tuple
#=====================================================================#

import numpy as np


class List:
    """Datatype: unordered list of integer tuples.

    The Fortran `table(num_tuples, max_entries)` stores one tuple per
    column; this is preserved as a 2-D NumPy array of the same shape.
    """

    #===================================================================#
    # Constructor implementaion of list
    #
    # Fortran: interface list => create(num_tuples, max_entries)
    #===================================================================#
    def __init__(self, num_tuples, max_entries):
        self.num_tuples = num_tuples

        self.table = np.zeros((num_tuples, max_entries), dtype=int)
        self.num_entries = 0

    #===================================================================#
    # Add an entry into the list
    #===================================================================#
    def insert(self, tuple):
        # Check if tuple is in table
        self.num_entries = self.num_entries + 1
        self.table[:, self.num_entries - 1] = tuple
        return True

    #===================================================================#
    # Add an entry into the list
    #===================================================================#
    def add_entry(self, tuple):
        # Check if tuple is in table
        self.num_entries = self.num_entries + 1
        self.table[:, self.num_entries - 1] = tuple

    #===================================================================#
    # Get all the entries in the list as an array
    #===================================================================#
    def get_entries(self):
        entries = np.array(self.table[:, 0:self.num_entries])
        return entries
