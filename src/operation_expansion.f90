!=====================================================================!
! THE EXPANSION ALONG THE DESIGN.
!
! A jet is the derivatives of a quantity along a coordinate at one
! point, to an order: the coefficients of its Taylor expansion there,
! the expansion not summed. The design coordinate nu is a point of
! the design space, fixed at a value nu_0, and its discrete image is
! the expansion of an order at that value, so that the discrete
! manifold is the instants (and cells) times the orders of the
! expansion:
!
!             nu_0 ------>  nu  (the design, one point, expanded to order 3)
!
!   t_n  |  [u, u_t, u_tt]  [u, u_t, u_tt]'  [u, u_t, u_tt]''  [u, u_t, u_tt]'''
!    :   |        :
!   t_2  |  [u, u_t, u_tt]  [u, u_t, u_tt]'  [u, u_t, u_tt]''  [u, u_t, u_tt]'''
!   t_1  |  [u, u_t, u_tt]  [u, u_t, u_tt]'  [u, u_t, u_tt]''  [u, u_t, u_tt]'''
!        |
!   time v      order 0          order 1          order 2           order 3
!              (the state)     (d/dnu)          (d^2/dnu^2)       (d^3/dnu^3)
!
! Every point stores the jet of the unknowns along time, which the
! families connect; column 0 is the solution, column 1 its
! derivative along nu at every instant, column 2 the second
! derivative, and so on. The derivative of a discrete field along nu
! is the columns from the next on. The multipliers have the same
! columns: the equations' at every instant, the reactions', the
! gauges', and the design condition's, whose column m is minus the
! (m + 1)-th derivative of the objective.
!
! The columns follow from the equations holding at every nu: with
! F(u(nu), nu) = 0 and A the jacobian of F in u, A u' + dF/dnu = 0
! gives column 1 by one linear solve with A, A u'' + (terms in u',
! u) = 0 gives column 2 by one more solve with the same A, and so
! on; the jet arithmetic of util_derivative_terms evaluates F on the
! jets of u and nu and returns the jet of F, whose coefficient at
! order m is the right side of the m-th solve. Nothing is
! differentiated by hand and nothing is approximated by differences.
!
! SEVERAL DESIGN COORDINATES. An expansion is along a seed, a vector
! over the design coordinates: the derivative of each coordinate
! along the expansion. The expansion along a coordinate has that
! coordinate's unit vector as its seed; the expansion along a
! direction, the direction. Stated by its order alone, the discrete
! image is the expansion along each coordinate; stated with seeds,
! along each seed given.
!
! THE COLUMNS. With m the order and the expansions numbered from one,
!
!      column 0                 the value, of every expansion
!      column (e - 1) m + o     the derivative of order o, 1 to m,
!                               along the expansion e
!
! so that every quantity of the discrete solution - an unknown, a
! multiplier, a functional - is its columns, 1 + m times the
! expansions in number, and the problem without a design is the one
! column zero.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_expansion

  use util_precision, only : dp

  implicit none

  private
  public :: expansion

  !===================================================================!
  ! The expansion of an order at the design values, along its seeds:
  ! one column of `along` per expansion, the derivative of each design
  ! coordinate along it; `at` the value of each design coordinate and
  ! `name` its name. A manifold without a design has the expansion of
  ! order zero with no seed, no value and no name.
  !===================================================================!

  type :: expansion

     integer :: order = 0
     real(dp), allocatable :: along(:,:)
     real(dp), allocatable :: at(:)
     character(len=32), allocatable :: name(:)

   contains

     procedure :: num_designs
     procedure :: num_expansions
     procedure :: seed
     procedure :: values
     procedure :: num_columns
     procedure :: column

  end type expansion

  interface expansion
     module procedure create_expansion
     module procedure create_expansion_along
     module procedure create_expansion_along_several
  end interface expansion

contains

  !===================================================================!
  ! THE CONSTRUCTORS: of an order, along each design coordinate when
  ! no seed is given, along the seeds given otherwise, one per column,
  ! each a vector over the design coordinates. Invalid input: a
  ! negative order, no seed, or a seed that is the zero vector.
  !===================================================================!

  function create_expansion(order) result(this)

    integer, intent(in) :: order
    type(expansion) :: this

    character(len=250) :: message

    if (order < 0) then
       write(message,'(a,i0)') 'operation_expansion: the expansion along the design is of order zero at &
            &least; order = ', order
       error stop trim(message)
    end if
    this % order = order

  end function create_expansion

  function create_expansion_along(order, along) result(this)

    integer , intent(in) :: order
    real(dp), intent(in) :: along(:)
    type(expansion) :: this

    this = create_expansion_along_several(order, reshape(along, [size(along), 1]))

  end function create_expansion_along

  function create_expansion_along_several(order, along) result(this)

    integer , intent(in) :: order
    real(dp), intent(in) :: along(:,:)
    type(expansion) :: this

    integer :: e

    this = create_expansion(order)
    if (size(along, 2) < 1) then
       error stop 'operation_expansion: an expansion along directions names one direction at least'
    end if
    do e = 1, size(along, 2)
       if (all(along(:, e) == 0.0_dp)) then
          error stop 'operation_expansion: a direction of the expansion is a nonzero vector over the design coordinates'
       end if
    end do
    this % along = along

  end function create_expansion_along_several

  !===================================================================!
  ! The design coordinates the expansion is at, and the expansions:
  ! none of either without a design.
  !===================================================================!

  pure integer function num_designs(this)
    class(expansion), intent(in) :: this
    num_designs = 0
    if (allocated(this % at)) num_designs = size(this % at)
  end function num_designs

  pure integer function num_expansions(this)
    class(expansion), intent(in) :: this
    num_expansions = 0
    if (allocated(this % along)) num_expansions = size(this % along, 2)
  end function num_expansions

  ! the seed of the e-th expansion: the derivative of every design
  ! coordinate along it
  pure function seed(this, e) result(along)
    class(expansion), intent(in) :: this
    integer         , intent(in) :: e
    real(dp), allocatable :: along(:)
    if (.not. allocated(this % along)) then
       allocate(along(1), source=0.0_dp)
       return
    end if
    if (e < 1 .or. e > size(this % along, 2)) then
       error stop 'operation_expansion: the expansion named does not exist'
    end if
    along = this % along(:, e)
  end function seed

  ! the design vector: the value of every design coordinate, or the
  ! one value zero without a design, which no expression reads
  pure function values(this) result(v)
    class(expansion), intent(in) :: this
    real(dp), allocatable :: v(:)
    if (.not. allocated(this % at)) then
       allocate(v(1), source=0.0_dp)
       return
    end if
    v = this % at
  end function values

  !===================================================================!
  ! THE COLUMNS: their number, and the column of the derivative of
  ! order o along the expansion e - column zero at order zero, the
  ! value every expansion begins at.
  !===================================================================!

  pure integer function num_columns(this)
    class(expansion), intent(in) :: this
    num_columns = 1 + this % num_expansions() * this % order
  end function num_columns

  pure integer function column(this, e, o)
    class(expansion), intent(in) :: this
    integer         , intent(in) :: e, o
    column = 0
    if (o > 0) column = (e - 1) * this % order + o
  end function column

end module operation_expansion
