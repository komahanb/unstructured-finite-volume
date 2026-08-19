!=====================================================================!
! The dense direct solve on the tower.
!
! Gaussian elimination with partial pivoting as a concrete
! minimizer: the matrix is never handed in, it is PROBED from the
! inherited matvec, one basis direction per column - the same face
! jacobi probes by colour, asked densely -
!
!      A(:,j) = matvec(e_j)
!
! and the assembled system A x = rhs is eliminated exactly. This
! file states the elimination and nothing else: matvec and the
! norm are inherited from the attached operation, the achieved
! residual is measured through them, and no matrix representation
! is owned here - a matrix is a graph with weights on its edges,
! and class_graph_stencil already says so.
!
! A direct solve is one pass, not an iteration: the inherited
! tolerance and budget govern nothing here, and the one numerical
! law is the pivot's - a pivot at or below the singular tolerance
! is refused loudly, never limped past.
!
! The convenience adapter at the bottom is the door for callers
! holding a plain dense array - the GTI time drivers assemble
! small dense Jacobians column by column and need them solved. It
! is an ADAPTER INTO THE MINIMIZER TOWER, not a solver hierarchy:
! it lays the array on a stencil_operator - each entry one
! weighted edge, column to row - attaches a dense_direct to that
! statement, and solves through the same minimizer face every
! other solver answers. Sparse, iterative, and multigrid solves
! live in this same tower; nothing GTI-shaped lives here at all.
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
  public :: probed_dense_matrix

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
  ! Probe, eliminate, certify: the matrix from the inherited
  ! matvec, the solve by partial pivoting, and the achieved
  ! residual measured through the inherited norm.
  !===================================================================!

  subroutine dense_direct_solve(this, rhs, x, achieved)

    class(dense_direct), intent(inout) :: this
    real(dp), intent(in)    :: rhs(:)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out)   :: achieved

    real(dp), allocatable :: a(:,:), probe(:), y(:), r(:), row(:)
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
    ! Assemble by probing: one matvec per basis direction, each
    ! answer one column.
    !----------------------------------------------------------------!

    allocate(a(n, n), probe(n))

    do j = 1, n
       probe    = 0.0_dp
       probe(j) = 1.0_dp
       call this % matvec(probe, y)
       a(:, j) = y
    end do

    !----------------------------------------------------------------!
    ! Gaussian elimination with partial pivoting, on the probed
    ! copy and a copy of the right-hand side.
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
    ! The truth either way: what the attached operation says about
    ! the answer, in the tower's own norm.
    !----------------------------------------------------------------!

    call this % matvec(x, y)
    achieved = this % norm(rhs - y)

  end subroutine dense_direct_solve

  !===================================================================!
  ! The convenience adapter into the minimizer tower: a plain dense
  ! array becomes a stencil_operator - each entry one weighted
  ! edge, column to row, constants zero - a dense_direct attaches
  ! to that statement exactly as any minimizer attaches to any
  ! operation, and the solve runs through the same face. Probing
  ! the stencil reassembles the given entries exactly, so the
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
  ! The compiled-dense door: probe any operation into the dense
  ! matrix of its linear action, one apply per basis column - the
  ! same probing this seat's solve performs against an attached
  ! matvec, offered raw to callers that must hold the matrix whole:
  ! a transposed step Jacobian, a small compiled block. A door, not
  ! an algebra - the numbers leave as one plain array.
  !===================================================================!

  subroutine probed_dense_matrix(action, on, width, a)

    class(graph_operation), intent(in)  :: action
    class(directed_graph) , intent(in)  :: on
    integer               , intent(in)  :: width
    real(dp), allocatable , intent(out) :: a(:,:)

    type(field)     :: probe
    class(graph_field), allocatable :: answer
    type(set_graph) :: dom
    real(dp), allocatable :: e(:), y(:)
    integer :: n_dom, ncomp, j

    call action % domain(on, dom, n_dom)

    if (n_dom <= 0) then
       error stop 'dense_direct: the probed operation''s domain is empty'
    end if
    if (width <= 0 .or. mod(width, n_dom) /= 0) then
       error stop 'dense_direct: the probed width carries a whole number per member'
    end if

    ncomp = width / n_dom

    allocate(a(width, width), e(width))

    do j = 1, width
       e    = 0.0_dp
       e(j) = 1.0_dp
       probe = field('probe', dom, n_dom, ncomp=ncomp)
       call probe % set_real_vector(e)
       call action % apply(on, [probe], answer)
       call answer % get_real_vector(y)
       if (size(y) /= width) then
          error stop 'dense_direct: the probed operation answers its own width'
       end if
       a(:, j) = y
    end do

  end subroutine probed_dense_matrix

end module class_graph_dense_direct
