#=====================================================================#
# Conjugate Gradient linear solver class that uses the functionalities
# of assembler class in the iterative solution process.
#=====================================================================#

import numpy as np

from .interface_linear_solver import LinearSolver


class ConjugateGradient(LinearSolver):

    #===================================================================!
    # Constructor for linear solver
    #
    # Fortran: interface conjugate_gradient => construct(FVAssembler,
    #          max_it, max_tol, print_level)
    #===================================================================!
    def __init__(self, FVAssembler, max_it, max_tol, print_level):
        self.FVAssembler = FVAssembler
        self.max_it = max_it
        self.max_tol = max_tol
        self.print_level = print_level

    #===================================================================!
    # Iterative linear solution using conjugate gradient method
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
            res_file = open('cg.res', 'w')
            res_file.write(" iteration  residual\n")

        update_norm = np.finfo(np.float64).max
        iter = 1
        while update_norm > self.max_tol:

            xold[:] = x

            # Inner iterations with CG
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
    # Iterative linear solution using conjugate gradient method
    #===================================================================!
    def iterate(self, x, ss):
        # Start the iteration counter
        iter = 1

        # Memory allocations
        b = np.zeros_like(x)
        p = np.zeros_like(x)
        r = np.zeros_like(x)
        w = np.zeros_like(x)
        Ax = np.zeros_like(x)
        tmp = np.zeros_like(x)

        rho = np.zeros(2, dtype=np.float64)

        # Norm of the right hand side
        self.FVAssembler.get_source(tmp)
        # Add the additional right hand side supplied
        tmp = tmp + ss
        if self.FVAssembler.symmetry is False:
            self.FVAssembler.get_transpose_jacobian_vector_product(b, tmp)
        else:
            b = tmp.copy()
        bnorm = np.linalg.norm(b)

        # Homogeneous case
        if bnorm <= self.max_tol:
            x[:] = 0.0
            return iter

        # Norm of the initial residual
        self.FVAssembler.get_jacobian_vector_product(tmp, x)
        if self.FVAssembler.symmetry is False:
            self.FVAssembler.get_transpose_jacobian_vector_product(Ax, tmp)
        else:
            Ax = tmp.copy()
        r = b - Ax
        rnorm = np.linalg.norm(r)
        tol = rnorm / bnorm
        rho[2 - 1] = rnorm * rnorm

        # Apply Iterative scheme until tolerance is achieved
        while (tol > self.max_tol) and (iter < self.max_it):

            # step (a) compute the descent direction
            if iter == 1:
                # steepest descent direction p
                p = r.copy()
            else:
                # take a conjugate direction
                beta = rho[2 - 1] / rho[1 - 1]
                p = r + beta * p

            # step (b) compute the solution update
            self.FVAssembler.get_jacobian_vector_product(tmp, p)
            if self.FVAssembler.symmetry is False:
                self.FVAssembler.get_transpose_jacobian_vector_product(w, tmp)
            else:
                w = tmp.copy()

            # step (c) compute the step size for update
            alpha = rho[2 - 1] / np.dot(p, w)

            # step (d) Add dx to the old solution
            x[:] = x + alpha * p

            # step (e) compute the new residual
            r = r - alpha * w

            # step(f) update values before next iteration
            rnorm = np.linalg.norm(r)
            tol = rnorm / bnorm

            if self.print_level > 1:
                print(iter, tol, rnorm, rho)

            iter = iter + 1

            rho[1 - 1] = rho[2 - 1]
            rho[2 - 1] = rnorm * rnorm

        return iter
