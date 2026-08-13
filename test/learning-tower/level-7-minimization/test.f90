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
! derive it from R_flow, and predict/error remain lawless members
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
  ! map. This fixture never inspects R_flow, OP_PREDICT or OP_ERROR:
  ! it is opaque to the solver and to Level 7 alike.

  use iso_fortran_env  , only : dp => REAL64
  use learning_assert  , only : SLOT_W
  use graph_carrier    , only : member_set
  use graph_grammar    , only : graph, graph_field, graph_operation
  use class_graph_field, only : field

  implicit none

  private
  public :: affine_learning_residual, ROW_R

  integer, parameter :: ROW_R = 1

  type, extends(graph_operation) :: affine_learning_residual
     class(member_set), allocatable :: theta, y
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

  type(affine_learning_residual) function create_oracle(theta, y, &
       & x_data, y_data) result(this)
    class(member_set), intent(in) :: theta, y
    real(dp)         , intent(in) :: x_data, y_data
    allocate(this % theta, source=theta)
    allocate(this % y    , source=y)
    this % x_data = x_data
    this % y_data = y_data
  end function create_oracle

  pure function oracle_name(this) result(name)
    class(affine_learning_residual), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'affine learning residual oracle'
  end function oracle_name

  subroutine oracle_domain(this, input_graph, domain)
    class(affine_learning_residual), intent(in) :: this
    class(graph), intent(in) :: input_graph
    class(member_set), allocatable, intent(out) :: domain
    associate (u1 => input_graph); end associate
    allocate(domain, source=this % y)
  end subroutine oracle_domain

  subroutine oracle_apply(this, input_graph, input_data, output)

    class(affine_learning_residual), intent(in)    :: this
    class(graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(field)                    :: out
    class(member_set), allocatable :: dom
    real(dp), allocatable          :: q(:), r(:)
    real(dp)                       :: w

    associate (u1 => input_graph); end associate

    call input_data(1) % domain(dom)
    if (.not. dom % same_as(this % theta)) then
       error stop 'oracle: the state must live on the trainable domain'
    end if
    call input_data(1) % get_real_vector(q)

    ! The parameter, through Theta's own enumeration - never by
    ! assuming member integer = vector position.
    w = q(this % theta % local_index(SLOT_W))

    allocate(r(1))
    r(this % y % local_index(ROW_R)) = this % x_data * w - this % y_data

    out = field('residual', this % y)
    call out % set_real_vector(r)
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine oracle_apply

end module learning_residual_fixture

program learning_level_7

  use iso_fortran_env  , only : dp => REAL64
  use learning_assert  , only : report, verdict, SLOT_W
  use graph_carrier    , only : counted_set, subset_set, member_set
  use graph_grammar    , only : graph_field
  use class_graph_field, only : field
  use class_graph      , only : stored_graph
  use class_graph_gmres, only : gmres
  use learning_residual_fixture, only : affine_learning_residual, ROW_R

  implicit none

  type(counted_set)     :: v, y, hv
  type(subset_set)      :: theta
  type(stored_graph)    :: host
  type(affine_learning_residual) :: oracle, witness
  type(gmres)           :: solver, fitter
  type(field)           :: state
  class(graph_field), allocatable :: resid
  class(member_set), allocatable  :: dom
  real(dp), allocatable :: g(:), rhs(:), w_state(:), rr(:)
  real(dp)              :: achieved
  integer               :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "learning tower . level 7 . minimization"
  write(*,'(1x,a)') "============================================="

  v     = counted_set('value-slots'  , 5)
  theta = subset_set('trainable', v, [SLOT_W])
  y     = counted_set('residual-rows', 1)

  ! Seven vertices: the wrong size for everything, on purpose.
  host = stored_graph(7, tails=[1,2,3,4,5,6], heads=[2,3,4,5,6,7])

  ! The primary training problem: (x, y) = (2, 6), supplied as data.
  oracle = affine_learning_residual(theta, y, 2.0_dp, 6.0_dp)

  call report(.not. theta % same_as(y), &
       & "Theta and Y are distinct domains, cardinality " // &
       & "notwithstanding", nfail)
  hv = host % vertex_set()
  call report(.not. hv % same_as(theta) .and. &
       &      host % num_vertices() /= theta % size(), &
       & "the host's seven vertices are nobody's trainables", nfail)

  call solver % attach(oracle, host, theta)
  solver % tolerance      = 1.0d-12
  solver % max_iterations = 50

  call solver % domain(host, dom)
  call report(dom % same_as(theta), &
       & "the solver's answer domain is Theta, by identity", nfail)
  call oracle % domain(host, dom)
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
  call report(abs(w_state(theta % local_index(SLOT_W)) - 3.0_dp) &
       &      < 1.0d-9, &
       & "learned w = 3, read through Theta's enumeration", nfail)
  call report(abs(w_state(theta % local_index(SLOT_W))) > 1.0_dp, &
       & "and the parameter actually MOVED: w = 3 is not w0 = 0", &
       & nfail)

  ! Direct confirmation through the oracle itself: R(w_learned) = 0.
  state = field('state', theta)
  call state % set_real_vector(w_state)
  call oracle % apply(host, [state], resid)
  call resid % get_real_vector(rr)
  call report(abs(rr(1)) < 1.0d-9, &
       & "the supplied residual vanishes at the learned parameter", &
       & nfail)

  ! The anti-hardcode witness: the SAME fixture, the SAME solver
  ! path, handed (x, y) = (4, 8) - and it fits w = 2, because the
  ! machinery solves x*w = y; it never recites 3.
  witness = affine_learning_residual(theta, y, 4.0_dp, 8.0_dp)
  call fitter % attach(witness, host, theta)
  fitter % tolerance      = 1.0d-12
  fitter % max_iterations = 50

  call fitter % constant(g)
  rhs = -g
  w_state = [0.0_dp]
  call fitter % solve(rhs, w_state, achieved)

  call report(achieved < 1.0d-10 .and. &
       &      abs(w_state(theta % local_index(SLOT_W)) - 2.0_dp) &
       &      < 1.0d-9, &
       & "(x,y) = (4,8) fits w = 2: fitted, never recited", nfail)

  call verdict(nfail, "level 7")

end program learning_level_7
