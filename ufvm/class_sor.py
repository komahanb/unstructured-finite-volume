#=====================================================================#
# Successive Over-Relaxation (SOR) linear solver class that uses the
# functionalities of assembler class in the iterative solution process.
#=====================================================================#

import numpy as np

from .interface_linear_solver import LinearSolver


class Sor(LinearSolver):

    #===================================================================!
    # Constructor for linear solver
    #
    # Fortran: interface sor => construct(FVAssembler, omega, max_it,
    #          max_tol, print_level)
    #===================================================================!
    def __init__(self, FVAssembler, omega, max_it, max_tol, print_level):
        self.FVAssembler = FVAssembler
        self.max_it = max_it
        self.max_tol = max_tol
        self.print_level = print_level
        self.omega = omega

    #===================================================================!
    # Estimate spectral radius using power iteration (maximum absolute
    # eigen value)
    #===================================================================!
    def estimate_spectral_radius(self, max_iter):
        # Create a random unit vector
        v = self.FVAssembler.create_vector()
        v[:] = np.random.random_sample(v.shape)
        v = v / np.linalg.norm(v)

        # A temp vector for processing
        w = self.FVAssembler.create_vector()

        mu = 0.0
        for iter in range(1, max_iter + 1):
            self.FVAssembler.get_jacobian_vector_product(w, v)
            wnorm = np.linalg.norm(w)
            v = w / wnorm
            mu = np.dot(v, w)
            print(iter, mu)

        return mu

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
            res_file = open('sor.res', 'w')
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
    # SOR inner iteration
    #===================================================================!
    def iterate(self, x, ss):
        b = np.zeros_like(x)
        Ux = np.zeros_like(x)
        D = np.zeros_like(x)
        Dx = np.zeros_like(x)
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

            # Form Ux
            self.FVAssembler.get_jacobian_vector_product(
                Ux, x, filter=self.FVAssembler.UPPER_TRIANGLE)

            # Form Dx
            self.FVAssembler.get_jacobian_vector_product(
                Dx, x, filter=self.FVAssembler.DIAGONAL)

            # R = w(b-Ux_k)+(1-w)Dx_k
            R = self.omega * (b - Ux) + (1.0 - self.omega) * Dx

            #--------------------------------------------------------------!
            # Solve the linear system: By=R ; (D+wL)y=R ; Dy=R-Ly
            #--------------------------------------------------------------!
            # Initial guess is the current solution
            y[:] = x

            iter2 = 1
            tol2 = np.finfo(np.float64).max
            while (tol2 > self.max_tol) and (iter2 < self.max_it):

                self.FVAssembler.get_jacobian_vector_product(
                    Ly, y, filter=self.FVAssembler.LOWER_TRIANGLE)

                ynew = (R - self.omega * Ly) / D
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
