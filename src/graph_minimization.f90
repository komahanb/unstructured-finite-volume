!=====================================================================!
! LEVEL 7 OF THE NEW TOWER . THE MINIMIZATION
!
! The first level with a goal, and its law in one line: given a
! residual map R : U -> Y, vary the values on the UNKNOWN domain U
! to drive the values on the RESIDUAL domain Y toward zero. U and Y
! are set graph identities - never assumed to be anyone's
! vertices; the graph argument survives only as the legacy
! operation host the compatibility apply() signature still wants. This
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
  use graph_operation_view  , only : graph_operation
  use graph_ordinary_view   , only : ordinary_graph
  use graph_field_calculus  , only : graph_field
  use fractal_graph      , only : set_graph => graph
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

     ! THE EXECUTION CONTEXT. The graph handed to the action when it
     ! is applied, and nothing more. It is whatever the action needs
     ! to compute with - a mesh, a compatibility host, a conduit -
     ! and it carries no authority over the solver's own structure.
     class(ordinary_graph)          , allocatable :: on

     ! THE DEPENDENT-VARIABLE COUPLING. Which unknowns feed which:
     ! the stencil the structural algorithms need, and the ONLY thing
     ! sweep_order is permitted to colour.
     !
     ! It is OPTIONAL AT ATTACH AND HAS NO FALLBACK. A structure-free
     ! minimizer - gmres, conjugate gradient, newton - never asks for
     ! it and need never supply one. A structured one - jacobi,
     ! gauss-seidel - is handed it by a caller that knows which object
     ! owns the dependent axis, and fails loudly if it was not.
     !
     ! `coupling := on` would be exactly the mistake this seat exists
     ! to prevent: the graph an action happens to execute over is not
     ! evidence about which unknowns are coupled. Where the two really
     ! are the same graph, the CALLER says so, at its own call site.
     class(ordinary_graph)          , allocatable :: coupling

     ! The unknown domain U: where the answer lives, explicit at
     ! attach, identity preserved - never inferred from the host.
     type(set_graph) :: unknown_domain
     integer         :: n_unknown_domain = 0

     ! The residual domain Y: what the action answers on, asked of
     ! the action itself at attach.
     type(set_graph) :: residual_domain
     integer         :: n_residual_domain = 0

     ! A second seat, one member per NUMBER rather than per cell.
     ! The pairings live here: a measure carries one weight per
     ! entry, so a dot product over wide entries must be taken on
     ! the values themselves - the calculus says as much in its own
     ! banner. With one number per cell the two seats are the same
     ! set, and nothing changes.
     type(set_graph) :: numbers
     integer         :: n_numbers = 0

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

  subroutine attach(this, action, on, unknown_domain, n_unknown_domain, &
       & ncomp, coupling)

    class(minimizer)  , intent(inout) :: this
    class(graph_operation), intent(in)    :: action
    class(ordinary_graph)          , intent(in)    :: on
    type(set_graph)       , intent(in)    :: unknown_domain
    integer               , intent(in)    :: n_unknown_domain
    integer, intent(in), optional         :: ncomp
    class(ordinary_graph)  , intent(in), optional  :: coupling

    real(dp), allocatable :: zero(:)
    integer :: n

    if (allocated(this % action)) deallocate(this % action)
    allocate(this % action, source=action)
    if (allocated(this % on)) deallocate(this % on)
    allocate(this % on, source=on)

    ! The dependent-variable coupling arrives EXPLICIT or not at all.
    ! No fallback to the execution context: a solver that needs
    ! structure and was given none says so when it reaches for it,
    ! rather than colouring whatever graph happened to be nearby.
    if (allocated(this % coupling)) deallocate(this % coupling)
    if (present(coupling)) allocate(this % coupling, source=coupling)

    this % ncomp = 1
    if (present(ncomp)) this % ncomp = max(ncomp, 1)

    ! The unknown domain arrives EXPLICIT and identity-preserving.
    ! No hidden fallback to the host's vertices: a caller that
    ! means vertices says so at its own call site.
    this % unknown_domain   = unknown_domain
    this % n_unknown_domain = n_unknown_domain

    ! The residual domain is the action's own answer.
    call action % domain(on, this % residual_domain, this % n_residual_domain)

    !----------------------------------------------------------------!
    ! attach is re-enterable - Newton calls it once per iteration - and
    ! a graph signs ONCE. The old counted_set constructor minted a
    ! fresh number domain on every attach, so a fresh one is minted
    ! here too, by resetting the component to an unsigned graph before
    ! declaring it. Signing the same variable twice is refused, and
    ! rightly; this says which of the two meanings was intended.
    !----------------------------------------------------------------!

    n = this % n_unknown_domain

    block
      type(set_graph) :: unsigned
      this % numbers = unsigned
    end block
    call this % numbers % declare()

    this % n_numbers = n * this % ncomp

    allocate(zero(n * this % ncomp))
    zero = 0.0_dp
    call raw_apply(this, zero, this % affine)

    ! Classification is not admissibility: U and Y stay distinct
    ! identities, but THIS solver family is square - the scalar
    ! dimensions must agree. A rectangular least-squares family may
    ! earn R^n -> R^m later; it has not yet.
    if (size(this % affine) /= n * this % ncomp) then
       error stop 'minimization: the current solver family requires equal &
            &unknown and residual value dimensions'
    end if

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

    state = field('state', this % unknown_domain, this % n_unknown_domain, ncomp=this % ncomp)
    call state % set_real_vector(x)

    call this % action % apply(this % on, [state], answer)

    block
      type(set_graph) :: got
      integer         :: n_got
      got   = answer % domain()
      n_got = answer % num_entries()
      if (.not. got % same_as(this % residual_domain)) then
         error stop 'minimization: the action must answer on its stated residual domain'
      end if
    end block

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

    uf = field('u', this % numbers, this % n_numbers)
    call uf % set_real_vector(u)
    vf = field('v', this % numbers, this % n_numbers)
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

    uf = field('u', this % numbers, this % n_numbers)
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

    ! THE COLOURING IS OF THE UNKNOWNS' COUPLING, never of the
    ! execution context. Two unknowns may share a colour only when
    ! nothing couples them, and the graph an action runs over knows
    ! nothing about that.
    if (.not. allocated(this % coupling)) then
       error stop 'minimization: a sweep needs the dependent-variable &
            &coupling - attach it with coupling='
    end if

    colouring = walk(WALK_COLOURING)
    call colouring % apply(this % coupling, output=answer)
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

  subroutine solver_domain(this, input_graph, domain, nentries)

    class(minimizer), intent(in)       :: this
    class(ordinary_graph), intent(in)               :: input_graph
    type(set_graph), intent(out) :: domain
    integer        , intent(out) :: nentries

    associate (u1 => input_graph); end associate

    ! The solver's answer is a solution on U.
    domain   = this % unknown_domain
    nentries = this % n_unknown_domain

  end subroutine solver_domain

  subroutine solver_apply(this, input_graph, input_data, output)

    class(minimizer), intent(in)                   :: this
    class(ordinary_graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    class(minimizer), allocatable :: worker
    type(field) :: out
    real(dp), allocatable :: rhs(:), x(:)
    real(dp) :: achieved

    associate (u1 => input_graph); end associate

    ! x IS a state on the unknown domain; say so.
    allocate(x(this % n_unknown_domain * this % ncomp))
    x = 0.0_dp

    if (present(input_data)) then
       block
         type(set_graph) :: got
         integer         :: n_got
         got   = input_data(1) % domain()
         n_got = input_data(1) % num_entries()
         if (.not. got % same_as(this % residual_domain)) then
            error stop 'minimization: a right-hand side lives on the residual domain'
         end if
       end block
       call input_data(1) % get_real_vector(rhs)
       allocate(worker, source=this)
       call worker % solve(rhs, x, achieved)
    end if

    out = field('solution', this % unknown_domain, this % n_unknown_domain, ncomp=this % ncomp)
    call out % set_real_vector(x)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine solver_apply

end module graph_minimization
