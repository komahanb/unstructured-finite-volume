!=====================================================================!
! THE EIGENSYSTEM OF A SYMMETRIC MATRIX, A = V diag(lambda) V^T, by
! the cyclic Jacobi method.
!
! Each rotation in the plane (p, q) sets the entry a_pq to zero.
! A sweep rotates in every plane once.
! The sum of the squares of the entries off the diagonal falls at
! every rotation, and quadratically once the eigenvalues are apart.
! An entry off the diagonal is negligible where adding it to each of
! its two diagonal entries changes neither.
! The method stops after a sweep in which every entry was negligible.
! The eigenvalues are then accurate to rounding relative to the
! largest, and they are returned in increasing order with their
! vectors in the columns of V.
!
! Invalid input: a matrix that is not square. A matrix not reduced
! after as many sweeps as there are bits in a real stops the program.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module util_eigensystem

  use util_precision, only : dp

  implicit none

  private
  public :: symmetric_eigensystem

contains

  pure subroutine symmetric_eigensystem(a, values, vectors)

    real(dp)             , intent(in)  :: a(:,:)
    real(dp), allocatable, intent(out) :: values(:)
    real(dp), allocatable, intent(out) :: vectors(:,:)

    real(dp), allocatable :: s(:,:), column(:)
    real(dp) :: tau, t, c, sn, x, y, eigenvalue
    integer :: n, p, q, r, sweep, i, j
    logical :: rotated

    n = size(a, 1)
    if (size(a, 2) /= n) error stop 'util_eigensystem: the eigensystem is of a square matrix'
    s = (a + transpose(a)) / 2.0_dp
    allocate(vectors(n, n), source=0.0_dp)
    do p = 1, n
       vectors(p, p) = 1.0_dp
    end do

    rotated = n > 1
    sweep   = 0
    do while (rotated)
       sweep = sweep + 1
       if (sweep > digits(1.0_dp)) then
          error stop 'util_eigensystem: the matrix is not reduced after as many sweeps as there are bits in a real'
       end if
       rotated = .false.
       do p = 1, n - 1
          do q = p + 1, n
             if (s(p, q) == 0.0_dp) cycle
             if (abs(s(p, p)) + abs(s(p, q)) == abs(s(p, p)) .and. &
                  & abs(s(q, q)) + abs(s(p, q)) == abs(s(q, q))) then
                s(p, q) = 0.0_dp
                s(q, p) = 0.0_dp
                cycle
             end if
             rotated = .true.
             ! the rotation J with J_pp = J_qq = c, J_pq = -J_qp = sn,
             ! whose J^T S J has a zero at (p, q)
             tau = (s(q, q) - s(p, p)) / (2.0_dp * s(p, q))
             t   = sign(1.0_dp, tau) / (abs(tau) + hypot(1.0_dp, tau))
             c   = 1.0_dp / hypot(1.0_dp, t)
             sn  = t * c
             do r = 1, n
                x = s(r, p)
                y = s(r, q)
                s(r, p) = c * x - sn * y
                s(r, q) = sn * x + c * y
             end do
             do r = 1, n
                x = s(p, r)
                y = s(q, r)
                s(p, r) = c * x - sn * y
                s(q, r) = sn * x + c * y
             end do
             s(p, q) = 0.0_dp
             s(q, p) = 0.0_dp
             do r = 1, n
                x = vectors(r, p)
                y = vectors(r, q)
                vectors(r, p) = c * x - sn * y
                vectors(r, q) = sn * x + c * y
             end do
          end do
       end do
    end do

    values = [(s(i, i), i = 1, n)]
    do i = 2, n
       eigenvalue = values(i)
       column     = vectors(:, i)
       j = i - 1
       do while (j >= 1)
          if (values(j) <= eigenvalue) exit
          values(j + 1)     = values(j)
          vectors(:, j + 1) = vectors(:, j)
          j = j - 1
       end do
       values(j + 1)     = eigenvalue
       vectors(:, j + 1) = column
    end do

  end subroutine symmetric_eigensystem

end module util_eigensystem
