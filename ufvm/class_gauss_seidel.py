#=====================================================================#
# Gauss-Seidel linear solver class that uses the functionalities of
# assembler class in the iterative solution process.
#=====================================================================#

import numpy as np

from .interface_linear_solver import LinearSolver


class GaussSeidel(LinearSolver):

    #===================================================================!
    # Constructor for linear solver
    #===================================================================!
    def __init__(self, FVAssembler, max_it, max_tol, print_level):
        self.FVAssembler = FVAssembler
        self.max_it = max_it
        self.max_tol = max_tol
        self.print_level = print_level

    #===================================================================!
    # Outer (skew-source) iterations
    #===================================================================!
    def solve(self):
        x = np.zeros(self.FVAssembler.num_state_vars, dtype=np.float64)
        self.FVAssembler.get_source(x)
        if np.linalg.norm(x) < np.finfo(np.float64).eps:
            print('zero rhs? stopping')
            raise SystemExit

        xold = np.zeros(self.FVAssembler.num_state_vars, dtype=np.float64)
        ss = np.zeros(self.FVAssembler.num_state_vars, dtype=np.float64)

        res_file = None
        if self.print_level == -1:
            res_file = open('gs.res', 'w')
            res_file.write(" iteration  residual\n")

        update_norm = np.finfo(np.float64).max
        iter = 1
        while update_norm > self.max_tol:

            xold[:] = x

            if iter == 1:
                ss[:] = 0.0
                num_inner_iters = self.iterate(x, ss)
            else:
                self.FVAssembler.get_skew_source(ss, x)
                num_inner_iters = self.iterate(x, ss)

            update_norm = np.linalg.norm(x - xold)
            if self.print_level == -1:
                res_file.write("%d %s\n" % (iter, repr(float(update_norm))))
            if self.print_level > 0:
                print("outer iter", iter, "num_inner_iters", num_inner_iters,
                      "update norm", update_norm)
            iter = iter + 1

        if self.print_level == -1:
            res_file.close()

        return x

    #===================================================================!
    # Gauss-Seidel inner iteration
    #===================================================================!
    def iterate(self, x, ss):
        b = np.zeros_like(x)
        Ux = np.zeros_like(x)
        D = np.zeros_like(x)
        identity = np.zeros_like(x)

        y = np.zeros_like(x)
        Ly = np.zeros_like(x)

        # Identity vector
        identity[:] = 1.0

        # Extract the diagonal entries
        self.FVAssembler.get_jacobian_vector_product(
            D, identity, filter=self.FVAssembler.DIAGONAL)

        # Assemble RHS of the linear system (source + boundary terms)
        self.FVAssembler.get_source(b)

        # Add the skew source terms if supplied
        b = b + ss

        bnorm = np.linalg.norm(b)

        # Homogeneous case (nothing to do)
        if bnorm <= self.max_tol:
            x[:] = 0.0
            return 1

        iter = 1
        tol = np.finfo(np.float64).max
        while (tol > self.max_tol) and (iter < self.max_it):

            # Form the residual (partial) after split
            self.FVAssembler.get_jacobian_vector_product(
                Ux, x, filter=self.FVAssembler.UPPER_TRIANGLE)

            R = b - Ux

            #--------------------------------------------------------------!
            # Solve the linear system: By=R ; (D+L)y=R ; Dy=R-Ly
            #--------------------------------------------------------------!
            # Initial guess is the current solution
            y[:] = x

            iter2 = 1
            tol2 = np.finfo(np.float64).max
            while (tol2 > self.max_tol) and (iter2 < self.max_it):

                self.FVAssembler.get_jacobian_vector_product(
                    Ly, y, filter=self.FVAssembler.LOWER_TRIANGLE)

                ynew = (R - Ly) / D
                tol2 = np.linalg.norm(y - ynew)

                if self.print_level > 2:
                    print("inner (2)", iter2, tol2)

                y[:] = ynew
                iter2 = iter2 + 1

            xnew = y.copy()
            tol = np.linalg.norm(x - xnew)

            if self.print_level > 1:
                print("inner (1)", iter, tol)

            x[:] = xnew
            iter = iter + 1

        return iter
