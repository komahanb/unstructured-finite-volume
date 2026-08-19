!=====================================================================!
! Dense direct solve: Gaussian elimination with partial pivoting
! as a concrete minimizer. The matrix is not passed in; it is
! assembled by applying the attached operation's matvec to each
! basis vector,
!
!      A(:,j) = matvec(e_j),
!
! and A x = rhs is then eliminated. matvec and the norm come from
! the attached operation, and the achieved residual is measured
! through them; no matrix representation is owned here, since
! class_graph_stencil already provides one.
!
! A direct solve is a single pass, so the tolerance and iteration
! budget inherited from the minimizer family are unused. The one
! numerical check is on the pivot: a pivot at or below
! singular_tolerance stops the program, because continuing would
! divide by a value indistinguishable from zero.
!
! solve_dense_matrix_with_dense_direct serves callers holding a
! plain dense array: the array is laid on a stencil_operator -
! each entry one weighted edge, column to row - a dense_direct is
! attached to that statement, and the solve runs through the same
! minimizer interface as every other solver. dense_matrix_of is
! the reverse direction: it assembles any operation's matrix
! column by column and returns it as a plain array.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_dense_direct

  use iso_fortran_env    , only : dp => REAL64
  use graph_operation_view, only : graph_operation
  use graph_directed_view , only : directed_graph
  use graph_field_calculus, only : graph_field
  use fractal_graph      , only : set_graph => graph
  use graph_minimization , only : minimizer
  use class_graph_stencil, only : stencil_operator
  use class_graph_field  , only : field

  implicit none

  private
  public :: dense_direct
  public :: solve_dense_matrix_with_dense_direct
  public :: dense_matrix_of

  type, extends(minimizer) :: dense_direct

     real(dp) :: singular_tolerance = 1.0e-14_dp

   contains

     procedure :: name  => dense_direct_name
     procedure :: solve => dense_direct_solve

  end type dense_direct

contains

  pure function dense_direct_name(this) result(name)

    class(dense_direct), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate

    name = 'dense direct'

  end function dense_direct_name

  !===================================================================!
  ! Solve A x = rhs where A(:, j) = matvec(e_j). Checks, each
  ! stopping the program: singular_tolerance must be positive,
  ! size(x) must equal size(rhs) because x is written in place,
  ! and every pivot must exceed singular_tolerance.
  !===================================================================!

  subroutine dense_direct_solve(this, rhs, x, achieved)

    class(dense_direct), intent(inout) :: this
    real(dp), intent(in)    :: rhs(:)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out)   :: achieved

    real(dp), allocatable :: a(:,:), basis(:), y(:), r(:), row(:)
    real(dp) :: swap_value, factor
    integer :: n, j, k, p, i

    if (this % singular_tolerance <= 0.0_dp) then
       error stop 'dense_direct: singular tolerance is positive'
    end if

    if (size(x) /= size(rhs)) then
       error stop 'dense_direct: solution size matches rhs'
    end if

    n = size(rhs)

    !----------------------------------------------------------------!
    ! Assemble the matrix: one matvec per basis vector, one column
    ! each.
    !----------------------------------------------------------------!

    allocate(a(n, n), basis(n))

    do j = 1, n
       basis    = 0.0_dp
       basis(j) = 1.0_dp
       call this % matvec(basis, y)
       a(:, j) = y
    end do

    !----------------------------------------------------------------!
    ! Gaussian elimination with partial pivoting, on the assembled
    ! matrix and a copy of the right-hand side.
    !----------------------------------------------------------------!

    r = rhs
    allocate(row(n))

    do k = 1, n

       p = k - 1 + maxloc(abs(a(k:n, k)), dim=1)

       if (abs(a(p, k)) <= this % singular_tolerance) then
          error stop 'dense_direct: pivot is nonsingular'
       end if

       if (p /= k) then
          row        = a(k, :)
          a(k, :)    = a(p, :)
          a(p, :)    = row
          swap_value = r(k)
          r(k)       = r(p)
          r(p)       = swap_value
       end if

       do i = k + 1, n
          factor    = a(i, k) / a(k, k)
          a(i, k:n) = a(i, k:n) - factor * a(k, k:n)
          r(i)      = r(i) - factor * r(k)
       end do

    end do

    do k = n, 1, -1
       x(k) = (r(k) - dot_product(a(k, k+1:n), x(k+1:n))) / a(k, k)
    end do

    !----------------------------------------------------------------!
    ! Measure the residual of the computed solution through the
    ! attached operation's matvec and norm.
    !----------------------------------------------------------------!

    call this % matvec(x, y)
    achieved = this % norm(rhs - y)

  end subroutine dense_direct_solve

  !===================================================================!
  ! Solve a plain dense array through the minimizer interface. The
  ! array must be square; a rectangular array stops the program.
  ! Each entry becomes one weighted edge (column to row, constants
  ! zero) of a stencil_operator, so assembling that statement's
  ! matrix reproduces the given entries exactly and the
  ! elimination sees precisely the caller's matrix.
  !===================================================================!

  subroutine solve_dense_matrix_with_dense_direct(a, b, singular_tolerance, &
       & x, achieved)

    real(dp)             , intent(in)  :: a(:,:)
    real(dp)             , intent(in)  :: b(:)
    real(dp)             , intent(in)  :: singular_tolerance
    real(dp), allocatable, intent(out) :: x(:)
    real(dp)             , intent(out) :: achieved

    type(stencil_operator) :: statement
    type(dense_direct)     :: solver

    integer , allocatable :: rows(:), columns(:)
    real(dp), allocatable :: weights(:), constant(:)
    integer :: n, i, j, e

    n = size(b)

    if (size(a, 1) /= n .or. size(a, 2) /= n) then
       error stop 'dense_direct: dense matrix is square'
    end if

    allocate(rows(n * n), columns(n * n), weights(n * n), constant(n))
    constant = 0.0_dp

    e = 0
    do j = 1, n
       do i = 1, n
          e          = e + 1
          rows(e)    = i
          columns(e) = j
          weights(e) = a(i, j)
       end do
    end do

    statement = stencil_operator(rows, columns, weights, constant, &
         & 'dense matrix')

    solver % singular_tolerance = singular_tolerance

    call solver % attach(statement, statement % pattern, &
         & statement % pattern % vertex_set(), &
         & statement % pattern % num_vertices())

    allocate(x(n))
    x = 0.0_dp

    call solver % solve(b, x, achieved)

  end subroutine solve_dense_matrix_with_dense_direct


  !===================================================================!
  ! Assemble the dense matrix of an operation's linear action, one
  ! apply per basis column, for callers that must hold the matrix
  ! whole (for example a transposed step Jacobian). width must be
  ! a positive whole multiple of the operation's domain size, and
  ! every apply must return exactly width values; both are checked
  ! and stop the program, because a mismatched column cannot be
  ! placed in the matrix.
  !===================================================================!

  subroutine dense_matrix_of(action, on, width, a)

    class(graph_operation), intent(in)  :: action
    class(directed_graph) , intent(in)  :: on
    integer               , intent(in)  :: width
    real(dp), allocatable , intent(out) :: a(:,:)

    type(field)     :: basis
    class(graph_field), allocatable :: output
    type(set_graph) :: dom
    real(dp), allocatable :: e(:), y(:)
    integer :: n_dom, ncomp, j

    call action % domain(on, dom, n_dom)

    if (n_dom <= 0) then
       error stop 'dense_direct: the operation''s domain is nonempty'
    end if
    if (width <= 0 .or. mod(width, n_dom) /= 0) then
       error stop 'dense_direct: the width carries a whole number per member'
    end if

    ncomp = width / n_dom

    allocate(a(width, width), e(width))

    do j = 1, width
       e    = 0.0_dp
       e(j) = 1.0_dp
       basis = field('basis', dom, n_dom, ncomp=ncomp)
       call basis % set_real_vector(e)
       call action % apply(on, [basis], output)
       call output % get_real_vector(y)
       if (size(y) /= width) then
          error stop 'dense_direct: the operation result matches the width'
       end if
       a(:, j) = y
    end do

  end subroutine dense_matrix_of

end module class_graph_dense_direct
