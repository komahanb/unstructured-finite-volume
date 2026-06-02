#=====================================================================#
# Gauss-Jacobi linear solver class that uses the functionalities of
# assembler class in the iterative solution process.
#=====================================================================#

import numpy as np

from .interface_linear_solver import LinearSolver


class GaussJacobi(LinearSolver):

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
        # Initial guess vector for the subspace is "b"
        x = np.zeros(self.FVAssembler.num_state_vars, dtype=np.float64)
        self.FVAssembler.get_source(x)
        if np.linalg.norm(x) < np.finfo(np.float64).eps:
            print('zero rhs? stopping')
            raise SystemExit

        xold = np.zeros(self.FVAssembler.num_state_vars, dtype=np.float64)
        ss = np.zeros(self.FVAssembler.num_state_vars, dtype=np.float64)

        res_file = None
        if self.print_level == -1:
            res_file = open('gj.res', 'w')
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
    # Gauss-Jacobi inner iteration
    #===================================================================!
    def iterate(self, x, ss):
        b = np.zeros_like(x)
        Ux = np.zeros_like(x)
        Lx = np.zeros_like(x)
        D = np.zeros_like(x)
        identity = np.zeros_like(x)

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

        #-----------------------------------------------------------------!
        # Apply Gauss Jacobi Iterative scheme until tolerance is achieved
        #-----------------------------------------------------------------!
        iter = 1
        tol = np.finfo(np.float64).max
        while (tol > self.max_tol) and (iter < self.max_it):

            # Form the residual (partial) after split
            self.FVAssembler.get_jacobian_vector_product(
                Ux, x, filter=self.FVAssembler.UPPER_TRIANGLE)

            self.FVAssembler.get_jacobian_vector_product(
                Lx, x, filter=self.FVAssembler.LOWER_TRIANGLE)

            R = b - Lx - Ux

            # Invert diagonal: D^{-1}(b - Lx - Ux)
            xnew = R / D
            tol = np.linalg.norm(x - xnew)

            if self.print_level > 1:
                print("inner (1)", iter, tol)

            x[:] = xnew
            iter = iter + 1

        return iter
