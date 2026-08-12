!=====================================================================!
! LEVEL 2 OF THE STRATIFICATION . THE MINIMIZATION
!
! The first level with a goal: drive a residual toward zero. This
! module holds the minimizer base - ONE family for one story:
! attach a statement, drive its residual to zero. Linear solvers,
! newton, and whatever else minimizes a residual are its
! concretions; their differences are governance inside the family,
! never a second taxonomy. The solver vocabulary is defined here as
! thin delegations to engine seats, so the engine never learns a
! solver word and the solver never says apply or measure:
!
!      matvec ········· the operation applied, minus its constant
!      inner_product ·· a sum reduction with the second field as
!                       the measure
!      norm ··········· the norm reduction
!      sweep_order ···· the colouring walk
!      diagonal ······· matvec probed by colour: applied to one
!                       colour's indicator, the answer at a member
!                       IS its diagonal entry, because no two
!                       neighbours share a colour
!      constant ······· the affine part of the attached operation -
!                       what boundary values and sources contribute
!                       at zero state; the assembled right hand side
!                       is its negative
!
! A concrete solver works in plain arrays - fetched once, worked,
! written back once, as the field banner orders - and states only
! its iteration. Everything else is inherited from here.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_minimization

  use iso_fortran_env       , only : dp => REAL64
  use graph_grammar         , only : graph, graph_field, graph_operation
  use graph_carrier         , only : member_set, counted_set
  use graph_calculus        , only : graph_functional
  use class_graph_field     , only : field
  use class_graph_reduction , only : reduction, REDUCE_SUM, REDUCE_NORM
  use class_graph_walk      , only : walk, WALK_COLOURING

  implicit none

  private
  public :: minimizer

  !===================================================================!
  ! The base: an attached operation, the graph it reads, and the
  ! tolerances every iteration honours.
  !===================================================================!

  type, abstract, extends(graph_operation) :: minimizer

     class(graph_operation), allocatable :: action
     class(graph)          , allocatable :: on

     type(counted_set) :: cells

     ! A second seat, one member per NUMBER rather than per cell.
     ! The pairings live here: a measure carries one weight per
     ! entry, so a dot product over wide entries must be taken on
     ! the values themselves - the calculus says as much in its own
     ! banner. With one number per cell the two seats are the same
     ! set, and nothing changes.
     type(counted_set) :: numbers

     ! How wide an entry is. One number per cell is the common case
     ! and the default; a state with several numbers per cell - a
     ! complex plane point, a species vector, a whole spatial field
     ! standing at one instant - says so at attach, and every word
     ! below then measures the entire vector instead of its first
     ! stripe.
     integer :: ncomp = 1

     real(dp), allocatable :: affine(:)

     integer  :: max_iterations = 1000
     real(dp) :: tolerance      = 1.0d-10

   contains

     procedure :: attach
     procedure :: matvec
     procedure :: inner_product
     procedure :: norm
     procedure :: sweep_order
     procedure :: diagonal
     procedure :: constant

     ! The operation face: a solver IS an operation - the one that
     ! answers the attached statement. apply solves from zero, so a
     ! solver composes wherever operations go; a preconditioner is
     ! exactly this face of an inner solver.
     procedure :: domain => solver_domain
     procedure :: apply  => solver_apply

     procedure(solve_interface), deferred :: solve

  end type minimizer

  abstract interface

     !----------------------------------------------------------------!
     ! Drive || rhs - matvec(x) || under the tolerance, within the
     ! iteration budget. The achieved norm reports the truth either
     ! way.
     !----------------------------------------------------------------!

     subroutine solve_interface(this, rhs, x, achieved)
       import :: minimizer, dp
       class(minimizer), intent(inout) :: this
       real(dp), intent(in)    :: rhs(:)
       real(dp), intent(inout) :: x(:)
       real(dp), intent(out)   :: achieved
     end subroutine solve_interface

  end interface

contains

  !===================================================================!
  ! Take the operation and the graph it reads. The affine part is
  ! measured here, once: the operation applied to nothing is what
  ! its boundary values and sources say by themselves.
  !===================================================================!

  subroutine attach(this, action, on, ncomp)

    class(minimizer)  , intent(inout) :: this
    class(graph_operation), intent(in)    :: action
    class(graph)          , intent(in)    :: on
    integer, intent(in), optional         :: ncomp

    real(dp), allocatable :: zero(:)
    integer :: v, nv

    if (allocated(this % action)) deallocate(this % action)
    allocate(this % action, source=action)
    if (allocated(this % on)) deallocate(this % on)
    allocate(this % on, source=on)

    this % ncomp = 1
    if (present(ncomp)) this % ncomp = max(ncomp, 1)

    nv = on % num_vertices()
    ! The solver's fields ride the statement graph's own carrier, so
    ! every kernel downstream recognizes their domain by identity.
    this % cells   = on % vertex_set()
    this % numbers = counted_set('numbers', nv * this % ncomp)

    allocate(zero(nv * this % ncomp))
    zero = 0.0_dp
    call raw_apply(this, zero, this % affine)

  end subroutine attach

  !===================================================================!
  ! The operation applied as it stands, affine part and all.
  !===================================================================!

  subroutine raw_apply(this, x, y)

    class(minimizer), intent(in)   :: this
    real(dp), intent(in)               :: x(:)
    real(dp), allocatable, intent(out) :: y(:)

    type(field) :: state
    class(graph_field), allocatable :: answer

    state = field('state', this % cells, ncomp=this % ncomp)
    call state % set_real_vector(x)

    call this % action % apply(this % on, [state], answer)
    call answer % get_real_vector(y)

  end subroutine raw_apply

  !===================================================================!
  ! The solver's words, each a delegation.
  !===================================================================!

  subroutine matvec(this, x, y)

    class(minimizer), intent(in)   :: this
    real(dp), intent(in)               :: x(:)
    real(dp), allocatable, intent(out) :: y(:)

    call raw_apply(this, x, y)
    y = y - this % affine

  end subroutine matvec

  real(dp) function inner_product(this, u, v) result(prod)

    class(minimizer), intent(in) :: this
    real(dp), intent(in) :: u(:), v(:)

    type(reduction) :: total
    type(field) :: uf, vf
    class(graph_functional), allocatable :: answer
    real(dp), allocatable :: got(:)

    uf = field('u', this % numbers)
    call uf % set_real_vector(u)
    vf = field('v', this % numbers)
    call vf % set_real_vector(v)

    total = reduction(REDUCE_SUM)
    call total % reduce(uf, answer, measure=vf)

    call answer % get_real_vector(got)
    prod = got(1)

  end function inner_product

  real(dp) function norm(this, u) result(length)

    class(minimizer), intent(in) :: this
    real(dp), intent(in) :: u(:)

    type(reduction) :: measure_of
    type(field) :: uf
    class(graph_functional), allocatable :: answer
    real(dp), allocatable :: got(:)

    uf = field('u', this % numbers)
    call uf % set_real_vector(u)

    measure_of = reduction(REDUCE_NORM)
    call measure_of % reduce(uf, answer)

    call answer % get_real_vector(got)
    length = got(1)

  end function norm

  subroutine sweep_order(this, colours)

    class(minimizer), intent(in)  :: this
    integer, allocatable, intent(out) :: colours(:)

    type(walk) :: colouring
    class(graph_field), allocatable :: answer

    colouring = walk(WALK_COLOURING)
    call colouring % apply(this % on, output=answer)
    call answer % get_integer_vector(colours)

  end subroutine sweep_order

  !===================================================================!
  ! The diagonal, probed by colour. Applied to the indicator of one
  ! colour class, the answer at a member is that member's diagonal
  ! entry, because none of its neighbours is in the class. One
  ! matvec per colour, a handful in all.
  !===================================================================!

  subroutine diagonal(this, d)

    class(minimizer), intent(in)   :: this
    real(dp), allocatable, intent(out) :: d(:)

    integer , allocatable :: colours(:)
    real(dp), allocatable :: indicator(:), y(:)
    integer :: nv, col, v

    if (this % ncomp > 1) then
       ! The probe reads one answer per cell, and a wide entry has
       ! several. A block probe is the honest generalization and no
       ! citizen has asked for one yet.
       error stop 'diagonal: the coloured probe answers one number per cell'
    end if

    nv = size(this % affine)
    allocate(d(nv), indicator(nv))
    d = 0.0_dp

    call this % sweep_order(colours)

    do col = 1, maxval(colours)

       indicator = 0.0_dp
       do v = 1, nv
          if (colours(v) == col) indicator(v) = 1.0_dp
       end do

       call this % matvec(indicator, y)

       do v = 1, nv
          if (colours(v) == col) d(v) = y(v)
       end do

    end do

  end subroutine diagonal

  !===================================================================!
  ! The affine part, for the caller assembling an equation: the
  ! statement action(q) = 0 reads matvec(q) = -constant.
  !===================================================================!

  subroutine constant(this, g)

    class(minimizer), intent(in)   :: this
    real(dp), allocatable, intent(out) :: g(:)

    g = this % affine

  end subroutine constant

  !===================================================================!
  ! The operation face.
  !===================================================================!

  subroutine solver_domain(this, input_graph, domain)

    class(minimizer), intent(in)       :: this
    class(graph), intent(in)               :: input_graph
    class(member_set), allocatable, intent(out) :: domain

    associate (u1 => this); end associate

    call input_graph % all_vertices(domain)

  end subroutine solver_domain

  subroutine solver_apply(this, input_graph, input_data, output)

    class(minimizer), intent(in)                   :: this
    class(graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    class(minimizer), allocatable :: worker
    type(field) :: out
    real(dp), allocatable :: rhs(:), x(:)
    real(dp) :: achieved

    associate (u1 => input_graph); end associate

    allocate(x(size(this % affine)))
    x = 0.0_dp

    if (present(input_data)) then
       call input_data(1) % get_real_vector(rhs)
       allocate(worker, source=this)
       call worker % solve(rhs, x, achieved)
    end if

    out = field('solution', this % cells, ncomp=this % ncomp)
    call out % set_real_vector(x)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine solver_apply

end module graph_minimization
