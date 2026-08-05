!=====================================================================!
! Restarted GMRES on the tower.
!
! For any operator, symmetric or not: build an orthonormal basis of
! the directions the operator itself suggests, and take the best
! answer that basis can express,
!
!      minimise || rhs - A x ||   over   x0 + span{r, Ar, A^2 r, ...}
!
! The basis grows one matvec at a time (Arnoldi, with the inner
! product keeping it orthonormal); a pair of rotations keeps the
! small problem triangular, so the residual norm is read off without
! solving anything until the end. After `restart` directions the
! basis is discarded and the walk begins again from the best point.
!
! Everything it needs is inherited: matvec, inner_product, norm.
! This file states the iteration and nothing else.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_gmres

  use iso_fortran_env  , only : dp => REAL64
  use graph_minimization, only : linear_solver

  implicit none

  private
  public :: gmres

  type, extends(linear_solver) :: gmres

     integer :: restart = 30

   contains

     procedure :: solve

  end type gmres

contains

  subroutine solve(this, rhs, x, achieved)

    class(gmres), intent(inout) :: this
    real(dp), intent(in)    :: rhs(:)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out)   :: achieved

    real(dp), allocatable :: basis(:,:), h(:,:), cs(:), sn(:), s(:)
    real(dp), allocatable :: r(:), w(:), y(:)
    real(dp) :: goal, beta, hik, radius
    integer :: n, m, outer, i, j, k

    n = size(x)
    m = max(this % restart, 1)

    allocate(basis(n, m + 1), h(m + 1, m), cs(m), sn(m), s(m + 1))

    goal = this % tolerance * (1.0_dp + this % norm(rhs))

    do outer = 1, this % max_iterations

       call this % matvec(x, r)
       r = rhs - r
       beta = this % norm(r)

       achieved = beta
       if (beta < goal) return

       basis = 0.0_dp
       h  = 0.0_dp
       s  = 0.0_dp
       s(1) = beta
       basis(:, 1) = r / beta

       do j = 1, m

          ! One more direction from the operator, kept orthonormal.
          call this % matvec(basis(:, j), w)
          do i = 1, j
             h(i, j) = this % inner_product(w, basis(:, i))
             w = w - h(i, j) * basis(:, i)
          end do
          h(j + 1, j) = this % norm(w)

          ! The rotations that came before, then the new one.
          do i = 1, j - 1
             hik        = cs(i) * h(i, j) + sn(i) * h(i + 1, j)
             h(i + 1, j) = -sn(i) * h(i, j) + cs(i) * h(i + 1, j)
             h(i, j)     = hik
          end do

          radius = sqrt(h(j, j)**2 + h(j + 1, j)**2)
          if (radius < tiny(1.0_dp)) then
             k = j
             exit
          end if
          cs(j) = h(j, j) / radius
          sn(j) = h(j + 1, j) / radius
          h(j, j)     = radius
          h(j + 1, j) = 0.0_dp

          s(j + 1) = -sn(j) * s(j)
          s(j)     =  cs(j) * s(j)

          k = j
          achieved = abs(s(j + 1))
          if (achieved < goal) exit
          if (h(j + 1, j) < tiny(1.0_dp)) exit
          if (j < m) basis(:, j + 1) = w / h(j + 1, j)

       end do

       ! The best answer this basis can express.
       allocate(y(k))
       do i = k, 1, -1
          y(i) = s(i)
          do j = i + 1, k
             y(i) = y(i) - h(i, j) * y(j)
          end do
          y(i) = y(i) / h(i, i)
       end do
       do i = 1, k
          x = x + y(i) * basis(:, i)
       end do
       deallocate(y)

       if (achieved < goal) then
          call this % matvec(x, r)
          achieved = this % norm(rhs - r)
          return
       end if

    end do

  end subroutine solve

end module class_graph_gmres
