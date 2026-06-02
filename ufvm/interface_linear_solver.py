#=====================================================================#
# Interface for linear solvers to implement.
#=====================================================================#

from abc import ABC, abstractmethod


class LinearSolver(ABC):
    """Linear solver datatype (abstract)."""

    # real(dp) :: max_tol
    # integer  :: max_it

    #===================================================================!
    # type bound procedures
    #
    # Fortran: subroutine solve_interface(this, x) with x intent(out).
    # The solution vector is returned.
    #===================================================================!
    @abstractmethod
    def solve(self):
        ...
