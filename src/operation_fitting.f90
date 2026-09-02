!=====================================================================!
! LEVEL 2 OF THE STRATIFICATION . THE FITTING FAMILY
!
! Two abstractions, from the form-and-coefficients principle: an
! expansion has a FORM - which functions are in the basis - and
! COEFFICIENTS - the weights computed within that form. The two
! change at different rates and belong to different types:
!
!      fit ·············· the coefficient sector, fast: given a form
!                         and a set of positions, find the
!                         weights that reproduce the target exactly
!                         on the form's span. A fit is an OPERATION
!                         on the point set, and its solve is a
!                         minimization solved by the level's own
!                         solver - never by an explicit formula.
!
!      form_optimizer ··· the form sector, slow: it CONTROLS fits,
!                         adjusting which basis members are active
!                         - removing members the points cannot
!                         resolve, adding members the residual
!                         requires. A form change is a re-typing
!                         event; the fit then operates within the
!                         new form.
!
! The form is defined one level DOWN, as its own type: a fit STORES a
! form as an operator stores coefficients. Polynomial or wave, the
! fit is independent of the basis - it evaluates the form it was
! given and computes the coefficients. One concrete fit; the
! variation is defined on the form.
!
!      B(m,j) = basis_m(x_j)        r(m) = scale * d(basis_m)/dn |at
!      (B B^T) lambda = r           w = B^T lambda
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_fitting

  use util_precision  , only : dp, spacing_at_one
  use operation_action, only : operation, contract
  use operation_action, only : binding, bound_real_vector
  use operation_action, only : emit
  use view_directed, only : directed_graph
  use field_calculus, only : field, FIELD_REAL
  use graph_fractal      , only : graph
  use field_forms        , only : form
  use field_stored  , only : stored_field
  use operation_stencil, only : stencil
  use operation_conjugate_gradient, only : conjugate_gradient

  implicit none

  private
  public :: fit
  public :: form_optimizer

  !===================================================================!
  ! The fit: one concrete. The target is stored as components, the shape
  ! is STORED - a level-1 form whose membership the form sector writes
  ! and this fit respects: the form IS a support, and its members
  ! specify which table entries are active.
  !===================================================================!

  type, extends(operation) :: fit

     class(form), allocatable :: shape

     real(dp), allocatable :: at(:)
     real(dp) :: direction(3) = [1.0_dp, 0.0_dp, 0.0_dp]
     real(dp) :: scale        = 1.0_dp

   contains

     procedure :: name   => fit_name
     procedure :: apply  => fit_apply

  end type fit

  interface fit
     module procedure create_fit
  end interface fit

  !===================================================================!
  ! The form optimizer: it stores no state of its own; it evaluates a
  ! form on a point set and adjusts the member set.
  !===================================================================!

  type, abstract :: form_optimizer

   contains

     procedure(adapt_interface), deferred :: adapt

  end type form_optimizer

  abstract interface

     subroutine adapt_interface(this, shape, positions)
       import :: form_optimizer, form, dp
       class(form_optimizer), intent(in) :: this
       class(form), intent(inout) :: shape
       real(dp), intent(in) :: positions(:)
     end subroutine adapt_interface

  end interface

!=====================================================================!
! The pruner: the form family's first concretion.
!
! The simplest form decision, made a type: a basis member the
! points cannot resolve - one whose column of values vanishes on the
! whole point set - is removed from the member set before any fit
! runs. What was previously a pivot test inside a solve is now a
! form change, owned by the sector that changes forms, at the lower
! rate of form changes.
!
!=====================================================================!

  public :: pruner

  type, extends(form_optimizer) :: pruner

     real(dp) :: threshold = spacing_at_one

   contains

     procedure :: adapt

  end type pruner

contains

  type(fit) function create_fit(shape, at, direction, scale) result(this)

    class(form), intent(in)        :: shape
    real(dp), intent(in)           :: at(:)
    real(dp), intent(in)           :: direction(3)
    real(dp), intent(in), optional :: scale

    allocate(this % shape, source=shape)
    this % at        = at
    this % direction = direction
    if (present(scale)) this % scale = scale

    ! ONE ARGUMENT: THE FIELD THE FORM IS FITTED TO. Its entries are
    ! the points of the point set and its components are their
    ! coordinates, so the count is the dimension of the space the
    ! form spans and not one. A form over a line takes one component
    ! and a form over a plane takes two, and stating one for both
    ! rejects every fit above a line.
    call this % declare_arguments(1, [contract(FIELD_REAL, shape % dimension())])

  end function create_fit

  pure function fit_name(this) result(name)

    class(fit), intent(in) :: this
    character(len=:), allocatable :: name

    name = 'fit over ' // merge('a stored form', 'no form yet  ', &
         & allocated(this % shape))

  end function fit_name
  !===================================================================!
  ! Positions in, weights out. The conditions respect the member set;
  ! the dual is passed to the level's own solver.
  !===================================================================!

  subroutine fit_apply(this, input_graph, inputs, output)

    class(fit), intent(in)                         :: this
    class(directed_graph), intent(in)                       :: input_graph
    type(binding), intent(in), optional       :: inputs(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field)   :: out
    type(stencil) :: dual
    type(conjugate_gradient) :: solver
    real(dp), allocatable :: positions(:), w(:), b(:,:), bw(:,:)
    real(dp), allocatable :: g(:,:), r(:), lam(:), point_weight(:)
    integer , allocatable :: member_list(:)
    logical , allocatable :: active(:)
    real(dp) :: achieved, d2, nearest
    integer :: npts, nc, i, j, d

    npts = input_graph % num_vertices()

    allocate(w(npts))
    w = 0.0_dp

    if (present(inputs)) then
       call bound_real_vector(inputs, this % argument(1), positions)

       d = this % shape % dimension()
       if (size(this % at) /= d .or. size(positions) /= d * npts) then
          error stop 'fitting: the form, the target and the positions read one dimension'
       end if

       nc = this % shape % num_members()
       allocate(b(nc, npts), g(nc, nc), r(nc), lam(nc))

       ! The distance weight: a point's weight is the inverse of its
       ! squared distance from the target; the target's own point, when
       ! it is a member, is weighted as the nearest neighbour.
       allocate(point_weight(npts))
       nearest = huge(1.0_dp)
       do j = 1, npts
          call this % shape % values(positions(d * (j - 1) + 1 : d * j), &
               & this % at, b(:, j))
          d2 = sum((positions(d * (j - 1) + 1 : d * j) - this % at)**2)
          point_weight(j) = d2
          if (d2 > 0.0_dp) nearest = min(nearest, d2)
       end do
       do j = 1, npts
          point_weight(j) = 1.0_dp / max(point_weight(j), nearest)
       end do

       call this % shape % slopes(this % at, this % at, &
            & this % direction, r)
       r = this % scale * r

       ! Membership is the member set: a table entry outside the form's
       ! member set contributes no condition and no right-hand side.
       call this % shape % members(member_list)
       allocate(active(nc))
       active = .false.
       do i = 1, size(member_list)
          if (member_list(i) >= 1 .and. member_list(i) <= nc) then
             active(member_list(i)) = .true.
          end if
       end do
       do i = 1, nc
          if (.not. active(i)) then
             b(i, :) = 0.0_dp
             r(i)    = 0.0_dp
          end if
       end do

       allocate(bw(nc, npts))
       do j = 1, npts
          bw(:, j) = b(:, j) * point_weight(j)
       end do
       g = matmul(bw, transpose(b))
       do i = 1, nc
          if (.not. active(i)) g(i, i) = 1.0_dp
       end do

       dual = stencil(g, label='fitting dual')

       call solver % attach(dual, dual % pattern, dual % pattern % vertex_set(), &
            & dual % pattern % num_vertices())
       solver % tolerance      = spacing_at_one
       solver % max_iterations = 50

       lam = 0.0_dp
       call solver % solve(r, lam, achieved)

       do j = 1, npts
          w(j) = 0.0_dp
          do i = 1, nc
             w(j) = w(j) + b(i, j) * lam(i)
          end do
          w(j) = w(j) * point_weight(j)
       end do

    end if

    out = stored_field('fit weights', input_graph % vertex_set(), input_graph % num_vertices())
    call out % set_real_vector(w)

    call emit(out, output)

  end subroutine fit_apply


  !===================================================================!
  ! Remove the members the points cannot resolve. The constant member
  ! is always retained: its value is nonzero at every point.
  !===================================================================!

  subroutine adapt(this, shape, positions)

    class(pruner), intent(in) :: this
    class(form), intent(inout) :: shape
    real(dp), intent(in) :: positions(:)

    real(dp), allocatable :: phi(:), column_norm(:), centroid(:)
    integer :: nc, npts, j, m, d

    d    = shape % dimension()
    nc   = shape % num_members()
    npts = size(positions) / d
    allocate(centroid(d))

    ! The members are evaluated about the point set's centroid, so the
    ! members the points cannot resolve do not depend on the point
    ! set's position.
    centroid = 0.0_dp
    do j = 1, npts
       centroid = centroid + positions(d * (j - 1) + 1 : d * j)
    end do
    centroid = centroid / real(max(npts, 1), dp)

    allocate(phi(nc), column_norm(nc))
    column_norm = 0.0_dp

    do j = 1, npts
       call shape % values(positions(d * (j - 1) + 1 : d * j), centroid, phi)
       do m = 1, nc
          column_norm(m) = column_norm(m) + phi(m) * phi(m)
       end do
    end do

    ! Membership is the member set: the pruned form is restricted to
    ! the retained table entries. A set needs no second list to specify
    ! its members, and the form sets its own - this procedure decides
    ! the members and does not write them directly.
    call shape % restrict(pack([(m, m = 1, nc)], column_norm > this % threshold * maxval(column_norm)))

  end subroutine adapt


end module operation_fitting
