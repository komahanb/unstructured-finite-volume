#=====================================================================#
# Basic implementation of an unordered set tuple
#=====================================================================#

import numpy as np


class Set:
    """Datatype: unordered set of integer tuples.

    The Fortran `table(num_tuples, max_entries)` stores one tuple per
    column; this is preserved as a 2-D NumPy array of the same shape.
    """

    #===================================================================#
    # Constructor implementaion of set
    #
    # Fortran: interface set => create(num_tuples, max_entries)
    #===================================================================#
    def __init__(self, num_tuples, max_entries):
        self.num_tuples = num_tuples

        self.table = np.zeros((num_tuples, max_entries), dtype=int)
        self.num_entries = 0

    #===================================================================#
    # Add an entry into the set
    #===================================================================#
    def insert(self, tuple):
        # Check if tuple is in table
        if self.contains(tuple) is False:
            self.num_entries = self.num_entries + 1
            self.table[:, self.num_entries - 1] = tuple
            return True
        else:
            return False

    #===================================================================#
    # Add an entry into the set
    #===================================================================#
    def add_entry(self, tuple):
        # Check if tuple is in table
        if self.contains(tuple) is False:
            self.num_entries = self.num_entries + 1
            self.table[:, self.num_entries - 1] = tuple

    #===================================================================#
    # Check if an entry is contains in the set
    #===================================================================#
    def contains(self, tuple):
        # loop through existing tuples and find if it exists
        for i in range(self.num_entries - 1, -1, -1):
            # Improve logic. This is expensive and unnnecessary
            if is_subset(tuple, self.table[:, i]) is True:
                return True
        return False

    #===================================================================#
    # Get all the entries in the set as an array
    #===================================================================#
    def get_entries(self):
        entries = np.array(self.table[:, 0:self.num_entries])
        return entries


#===================================================================#
# Checks if the first argument is a subset of the second argument
# (move elsewhere?)
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
        if not (set_ == sub[i]).any():    # any(set .eq. sub(i)) .eqv. .false.
            return False
    return True


#===================================================================#
# Sort an integer array ! move elsewhere?
#===================================================================#
def isort(array):
    # In-place bubble sort, mirroring the Fortran (array is intent(inout)).
    n = array.shape[0]
    for j in range(n):
        for k in range(j + 1, n):
            if array[j] > array[k]:
                temp = array[k]
                array[k] = array[j]
                array[j] = temp
