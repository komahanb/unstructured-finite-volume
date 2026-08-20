!=====================================================================!
! LEARNING TOWER . LEVEL 7 . MINIMIZATION
!
! The level answers one question: GIVEN a supplied residual map
! R : Theta -> Y, can the existing residual minimizer vary the
! trainable parameter until the residual vanishes? This is the
! FIRST rung where the parameter CHANGES:
!
!      Theta = { w } c--> V     the trainable domain, |Theta| = 1
!      Y     = { r }            the residual row,     |Y|     = 1
!      host                     a SEVEN-vertex legacy operation
!                               host, deliberately the wrong size
!                               for everything - operation-host
!                               compatibility debt, not a domain
!
! The residual formula r(w) = x_data*w - y_data is TEST-LOCAL
! ORACLE DATA supplied from above the frontier: Level 7 does not
! derive it from T_flow, and predict/error remain lawless members
! of O until Level 8. The solver is the ordinary GMRES citizen -
! attach, constant, solve - which discovers the linear action
! 2w from the opaque residual through the minimizer's own affine
! split R(w) - R(0). This is exact parameter fitting for an
! attainable square residual equation: one parameter, one residual.
! It is NOT gradient descent, NOT SGD, NOT backpropagation - no
! gradient of any loss exists here, and Level 6's structural
! J_Theta pattern is deliberately NOT consumed: structural sparsity
! and numerical operator action stay unconflated. The same fixture,
! handed (x,y) = (4,8) instead of (2,6), fits w = 2: the machinery
! solves x*w = y; it never recites 3.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module learning_residual_fixture

  ! The test-local oracle: R : field(Theta) -> field(Y), affine in
  ! the one trainable parameter. The data (x_data, y_data) live HERE,
  ! parametrically - never in production, never hard-coded in the
  ! map. This fixture never inspects T_flow, OP_PREDICT or OP_ERROR:
  ! it is opaque to the solver and to Level 7 alike.

  use iso_fortran_env  , only : dp => REAL64
  use learning_assert  , only : SLOT_W
  use graph_fractal        , only : graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map
  use map_set_representation, only : set_representation
  use map_inclusion  , only : inclusion_map, declared_subobject
  use operation_action, only : operation
  use view_directed, only : directed_graph
  use field_calculus, only : field
  use field_stored, only : stored_field

  implicit none

  private
  public :: affine_learning_residual, ROW_R

  integer, parameter :: ROW_R = 1

  type, extends(operation) :: affine_learning_residual
     ! Identity, count, and the action's OWN coordinates: a
     ! representation carries no identity, so holding one keeps this
     ! type free of any caller's map.
     type(graph) :: theta, y
     integer         :: n_theta = 0, n_y = 0
     class(set_representation), allocatable :: c_theta, c_y
     real(dp)                       :: x_data, y_data
   contains
     procedure :: name   => oracle_name
     procedure :: domain => oracle_domain
     procedure :: apply  => oracle_apply
  end type affine_learning_residual

  interface affine_learning_residual
     module procedure create_oracle
  end interface affine_learning_residual

contains

  type(affine_learning_residual) function create_oracle(theta, y, sets, &
       & x_data, y_data) result(this)
    type(graph), intent(in) :: theta, y
    type(set_map)  , intent(in) :: sets
    real(dp)       , intent(in) :: x_data, y_data
    this % theta = theta ; this % n_theta = sets % num_members_of(theta)
    this % y     = y     ; this % n_y     = sets % num_members_of(y)
    call sets % extent_of(theta, this % c_theta)
    call sets % extent_of(y    , this % c_y)
    this % x_data = x_data
    this % y_data = y_data
  end function create_oracle

  pure function oracle_name(this) result(name)
    class(affine_learning_residual), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'affine learning residual oracle'
  end function oracle_name

  subroutine oracle_domain(this, input_graph, domain, num_entries)
    class(affine_learning_residual), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    type(graph), intent(out) :: domain
    integer        , intent(out) :: num_entries
    integer         :: n_y
    associate (u1 => input_graph); end associate
    domain   = this % y
    num_entries = this % n_y
  end subroutine oracle_domain

  subroutine oracle_apply(this, input_graph, input_data, output)

    class(affine_learning_residual), intent(in)    :: this
    class(directed_graph), intent(in)                       :: input_graph
    class(field), intent(in), optional       :: input_data(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field)                    :: out
    type(graph) :: dom
    real(dp), allocatable          :: q(:), r(:)
    real(dp)                       :: w
    integer         :: n_y

    associate (u1 => input_graph); end associate

    dom = input_data(1) % domain()
    if (.not. dom % same_as(this % theta)) then
       error stop 'oracle: the state must live on the trainable domain'
    end if
    call input_data(1) % real_vector(q)

    ! The parameter, through Theta's own enumeration - never by
    ! assuming member integer = vector position.
    w = q(this % c_theta % local_index(SLOT_W))

    allocate(r(1))
    r(this % c_y % local_index(ROW_R)) = this % x_data * w - this % y_data

    out = stored_field('residual', this % y, this % n_y)
    call out % set_real_vector(r)
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine oracle_apply

end module learning_residual_fixture

program learning_level_7

  use iso_fortran_env  , only : dp => REAL64
  use learning_assert  , only : report, verdict, SLOT_W
  use field_calculus, only : field
  use field_stored, only : stored_field
  use view_directed_stored      , only : stored_directed_graph
  use operation_gmres, only : gmres
  use graph_fractal        , only : graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map
  use map_inclusion  , only : inclusion_map, declared_subobject
  use learning_residual_fixture, only : affine_learning_residual, ROW_R

  implicit none

  type(graph)     :: v, y, hv
  type(graph)      :: theta
  type(stored_directed_graph)    :: host
  type(affine_learning_residual) :: oracle, witness
  type(gmres)           :: solver, fitter
  type(stored_field)           :: state
  class(field), allocatable :: resid
  type(graph)  :: dom
  real(dp), allocatable :: g(:), rhs(:), w_state(:), rr(:)
  real(dp)              :: achieved
  integer               :: nfail
  type(inclusion_map)     :: inclusions
  integer         :: n_dom
  type(set_map)     :: sets

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "learning tower . level 7 . minimization"
  write(*,'(1x,a)') "============================================="

  call v % declare()
  call sets % bind(v, counted_set_representation(5))
  call theta % declare()
  call sets       % bind(theta, listed_set_representation([SLOT_W]))
  call inclusions % include_in(theta, v)
  call y % declare()
  call sets % bind(y, counted_set_representation(1))

  ! Seven vertices: the wrong size for everything, on purpose.
  host = stored_directed_graph(7, tails=[1,2,3,4,5,6], heads=[2,3,4,5,6,7])

  ! The primary training problem: (x, y) = (2, 6), supplied as data.
  oracle = affine_learning_residual(theta, y, sets, 2.0_dp, 6.0_dp)

  call report(.not. theta % same_as(y), &
       & "Theta and Y are distinct domains, cardinality " // &
       & "notwithstanding", nfail)
  hv = host % vertex_set()
  call report(.not. hv % same_as(theta) .and. &
       &      host % num_vertices() /= sets % num_members_of(theta), &
       & "the host's seven vertices are nobody's trainables", nfail)

  call solver % attach(oracle, host, theta, sets % num_members_of(theta))
  solver % tolerance      = 1.0d-12
  solver % max_iterations = 50

  call solver % domain(host, dom, n_dom)
  call report(dom % same_as(theta), &
       & "the solver's answer domain is Theta, by identity", nfail)
  call oracle % domain(host, dom, n_dom)
  call report(dom % same_as(y), &
       & "and the oracle answers on Y", nfail)

  ! The affine split, agreed on both sides: R(0) = -6, rhs = 6.
  call solver % constant(g)
  call report(size(g) .eq. 1 .and. abs(g(1) + 6.0_dp) < 1.0d-12, &
       & "the affine constant is -6: R at the untrained parameter", &
       & nfail)
  rhs = -g
  call report(abs(rhs(1) - 6.0_dp) < 1.0d-12, &
       & "so the equation to satisfy is matvec(w) = 6", nfail)

  ! w0 = 0: Level 5's initial parameter, carried forward in meaning.
  w_state = [0.0_dp]
  call solver % solve(rhs, w_state, achieved)

  call report(achieved < 1.0d-10, &
       & "the residual is driven below tolerance", nfail)
  call report(abs(w_state(sets % index_in(theta, SLOT_W)) - 3.0_dp) &
       &      < 1.0d-9, &
       & "learned w = 3, read through Theta's enumeration", nfail)
  call report(abs(w_state(sets % index_in(theta, SLOT_W))) > 1.0_dp, &
       & "and the parameter actually MOVED: w = 3 is not w0 = 0", &
       & nfail)

  ! Direct confirmation through the oracle itself: R(w_learned) = 0.
  state = stored_field('state', theta, sets % num_members_of(theta))
  call state % set_real_vector(w_state)
  call oracle % apply(host, [state], resid)
  call resid % real_vector(rr)
  call report(abs(rr(1)) < 1.0d-9, &
       & "the supplied residual vanishes at the learned parameter", &
       & nfail)

  ! The anti-hardcode witness: the SAME fixture, the SAME solver
  ! path, handed (x, y) = (4, 8) - and it fits w = 2, because the
  ! machinery solves x*w = y; it never recites 3.
  witness = affine_learning_residual(theta, y, sets, 4.0_dp, 8.0_dp)
  call fitter % attach(witness, host, theta, sets % num_members_of(theta))
  fitter % tolerance      = 1.0d-12
  fitter % max_iterations = 50

  call fitter % constant(g)
  rhs = -g
  w_state = [0.0_dp]
  call fitter % solve(rhs, w_state, achieved)

  call report(achieved < 1.0d-10 .and. &
       &      abs(w_state(sets % index_in(theta, SLOT_W)) - 2.0_dp) &
       &      < 1.0d-9, &
       & "(x,y) = (4,8) fits w = 2: fitted, never recited", nfail)

  call verdict(nfail, "level 7")

end program learning_level_7
