!=====================================================================!
! LEVEL 2 OF THE STRATIFICATION . THE FITTING FAMILY
!
! Two abstractions, one sentence from the form-and-coefficients
! doctrine: an expansion has a FORM - which functions stand in the
! basis - and COEFFICIENTS - the weights found inside that form.
! The two move at different speeds and belong to different citizens:
!
!      fit ·············· the coefficient sector, fast: given a form
!                         and a constellation of positions, find the
!                         weights that reproduce the target exactly
!                         on the form's span. A fit is an OPERATION
!                         on the point set, and its solve is a
!                         minimization answered by the level's own
!                         solver - never by hand.
!
!      form_optimizer ··· the form sector, slow: it GOVERNS fits,
!                         adjusting which basis members stand active
!                         - pruning what the points cannot see,
!                         admitting what the residual demands. A
!                         form change is a re-typing event; the fit
!                         then works inside the new form.
!
! The form lives one level DOWN, as its own citizen: a fit HOLDS a
! form the way an operator holds coefficients. Polynomial or wave,
! the fit does not care - it evaluates the shape it was handed and
! finds the coefficients. One concrete fit; the variation lives on
! the form, where it belongs.
!
!      B(m,j) = basis_m(x_j)        r(m) = scale * d(basis_m)/dn |at
!      (B B^T) lambda = r           w = B^T lambda
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_fitting

  use iso_fortran_env    , only : dp => REAL64
  use graph_grammar      , only : graph, graph_operation
  use graph_field_calculus, only : graph_field
  use fractal_graph      , only : set_graph => graph
  use graph_forms        , only : form
  use class_graph_field  , only : field
  use class_graph_stencil, only : stencil_operator
  use class_graph_conjugate_gradient, only : conjugate_gradient

  implicit none

  private
  public :: fit
  public :: form_optimizer

  !===================================================================!
  ! The fit: one concrete. The target rides as components, the shape
  ! is HELD - a level-1 form whose membership the form sector writes
  ! and this fit honours: the form IS a support, and its members say
  ! which table entries stand.
  !===================================================================!

  type, extends(graph_operation) :: fit

     class(form), allocatable :: shape

     real(dp) :: at(3)        = 0.0_dp
     real(dp) :: direction(3) = [1.0_dp, 0.0_dp, 0.0_dp]
     real(dp) :: scale        = 1.0_dp

   contains

     procedure :: name   => fit_name
     procedure :: domain => fit_domain
     procedure :: apply  => fit_apply

  end type fit

  interface fit
     module procedure create_fit
  end interface fit

  !===================================================================!
  ! The form optimizer: it holds no machinery of its own; it reads a
  ! form against a constellation and adjusts the roster.
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

contains

  type(fit) function create_fit(shape, at, direction, scale) result(this)

    class(form), intent(in)        :: shape
    real(dp), intent(in)           :: at(3)
    real(dp), intent(in)           :: direction(3)
    real(dp), intent(in), optional :: scale

    allocate(this % shape, source=shape)
    this % at        = at
    this % direction = direction
    if (present(scale)) this % scale = scale

  end function create_fit

  pure function fit_name(this) result(name)

    class(fit), intent(in) :: this
    character(len=:), allocatable :: name

    name = 'fit over ' // merge('a held form', 'no form yet', &
         & allocated(this % shape))

  end function fit_name

  subroutine fit_domain(this, input_graph, domain, nentries)

    class(fit), intent(in)                 :: this
    class(graph), intent(in)               :: input_graph
    type(set_graph), intent(out) :: domain
    integer        , intent(out) :: nentries

    associate (u1 => this); end associate

    domain   = input_graph % all_vertices()
    nentries = input_graph % num_vertices()

  end subroutine fit_domain

  !===================================================================!
  ! Positions in, weights out. The conditions honour the roster; the
  ! dual is handed to the level's own solver.
  !===================================================================!

  subroutine fit_apply(this, input_graph, input_data, output)

    class(fit), intent(in)                         :: this
    class(graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(field)   :: out
    type(stencil_operator) :: dual
    type(conjugate_gradient) :: solver
    real(dp), allocatable :: positions(:), w(:), b(:,:), bw(:,:)
    real(dp), allocatable :: g(:,:), r(:), lam(:), entries(:), price(:)
    integer , allocatable :: rows(:), columns(:), standing(:)
    logical , allocatable :: stands(:)
    real(dp) :: achieved, d2, nearest
    integer :: npts, nc, i, j, v

    npts = input_graph % num_vertices()

    allocate(w(npts))
    w = 0.0_dp

    if (present(input_data)) then

       call input_data(1) % get_real_vector(positions)

       nc = this % shape % size_of()
       allocate(b(nc, npts), g(nc, nc), r(nc), lam(nc))

       ! The distance metric: a point's share is priced by how far
       ! it stands from the target; the target's own point, when it
       ! is a member, is priced as the nearest neighbour.
       allocate(price(npts))
       nearest = huge(1.0_dp)
       do j = 1, npts
          call this % shape % values(positions(3 * j - 2 : 3 * j), &
               & this % at, b(:, j))
          d2 = sum((positions(3 * j - 2 : 3 * j) - this % at)**2)
          price(j) = d2
          if (d2 > 0.0_dp) nearest = min(nearest, d2)
       end do
       do j = 1, npts
          price(j) = 1.0_dp / max(price(j), nearest)
       end do

       call this % shape % slopes(this % at, this % at, &
            & this % direction, r)
       r = this % scale * r

       ! Membership is the roster: a table entry outside the form's
       ! member set carries no condition and no demand.
       call this % shape % members(standing)
       allocate(stands(nc))
       stands = .false.
       do i = 1, size(standing)
          if (standing(i) >= 1 .and. standing(i) <= nc) then
             stands(standing(i)) = .true.
          end if
       end do
       do i = 1, nc
          if (.not. stands(i)) then
             b(i, :) = 0.0_dp
             r(i)    = 0.0_dp
          end if
       end do

       allocate(bw(nc, npts))
       do j = 1, npts
          bw(:, j) = b(:, j) * price(j)
       end do
       g = matmul(bw, transpose(b))
       do i = 1, nc
          if (.not. stands(i)) g(i, i) = 1.0_dp
       end do

       allocate(rows(nc * nc), columns(nc * nc), entries(nc * nc))
       do i = 1, nc
          do j = 1, nc
             rows((i - 1) * nc + j)    = i
             columns((i - 1) * nc + j) = j
             entries((i - 1) * nc + j) = g(i, j)
          end do
       end do

       dual = stencil_operator(rows, columns, entries, &
            & [(0.0_dp, i = 1, nc)], label='fitting dual')

       call solver % attach(dual, dual % pattern, dual % pattern % vertex_set(), &
            & dual % pattern % num_vertices())
       solver % tolerance      = 1.0d-14
       solver % max_iterations = 50

       lam = 0.0_dp
       call solver % solve(r, lam, achieved)

       do j = 1, npts
          w(j) = 0.0_dp
          do i = 1, nc
             w(j) = w(j) + b(i, j) * lam(i)
          end do
          w(j) = w(j) * price(j)
       end do

    end if

    out = field('fit weights', input_graph % vertex_set(), input_graph % num_vertices())
    call out % set_real_vector(w)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine fit_apply

end module graph_fitting
