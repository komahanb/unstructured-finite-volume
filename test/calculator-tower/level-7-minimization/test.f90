!=====================================================================!
! CALCULATOR TOWER . LEVEL 7 . MINIMIZATION
!
! The level answers one question: GIVEN a residual map R : U -> Y,
! can ordinary minimization drive it to zero - without assuming the
! unknowns are anyone's vertices?
!
!      U = { e, c }  c--> X      the unknowns, DEPENDENCY order
!      Y = { r_c, r_e }          the residual rows of Level 6
!      host                     a SEVEN-vertex legacy operation
!                               host, deliberately the wrong size
!                               for everything - operation-host
!                               compatibility debt, not a domain
!
! The residual formulas r_c = q(c) - 5 and r_e = q(e) - 4 q(c) are
! TEST DATA supplied from above the frontier: Level 7 does not
! derive 5, 4, +, or multiplication, and their calculator meaning
! is not certified until Level 8. The solver is the ordinary GMRES
! citizen: attach, constant, solve - no calculator API, no Newton,
! no manual Jacobian, no arithmetic knowledge anywhere in
! production.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module affine_residual_fixture

  ! The test-local oracle: R : field(U) -> field(Y). These formulas
  ! are supplied data. This fixture never inspects T_flow, OP_PLUS
  ! or OP_TIMES - it is opaque to the solver and to Level 7 alike.

  use iso_fortran_env  , only : dp => REAL64
  use calculator_assert, only : SLOT_C, SLOT_E
  use fractal_graph        , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use graph_set_map        , only : set_map
  use graph_set_representation, only : set_representation
  use graph_inclusion_map  , only : inclusion_map, declared_subobject
  use graph_operation_view, only : graph_operation
  use graph_directed_view, only : directed_graph
  use graph_field_calculus, only : graph_field
  use class_graph_field, only : field

  implicit none

  private
  public :: affine_residual, ROW_C, ROW_E

  integer, parameter :: ROW_C = 1
  integer, parameter :: ROW_E = 2

  type, extends(graph_operation) :: affine_residual
     !----------------------------------------------------------------!
     ! Identity, count, and the action's OWN coordinates - copied in
     ! at construction, the way a CSR relation copies the numbering it
     ! will execute against. A representation is not a map: it carries
     ! no identity and answers only where a member stands, so holding
     ! one keeps this type free of any caller's interpretation.
     !----------------------------------------------------------------!
     type(set_graph) :: u, y
     integer         :: n_u = 0, n_y = 0
     class(set_representation), allocatable :: u_coords, y_coords
   contains
     procedure :: name   => oracle_name
     procedure :: domain => oracle_domain
     procedure :: apply  => oracle_apply
  end type affine_residual

  interface affine_residual
     module procedure create_oracle
  end interface affine_residual

contains

  type(affine_residual) function create_oracle(u, y, sets) result(this)
    type(set_graph), intent(in) :: u, y
    type(set_map)  , intent(in) :: sets
    this % u = u ; this % n_u = sets % size_of(u)
    this % y = y ; this % n_y = sets % size_of(y)
    call sets % extent_of(u, this % u_coords)
    call sets % extent_of(y, this % y_coords)
  end function create_oracle

  pure function oracle_name(this) result(name)
    class(affine_residual), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'affine residual oracle'
  end function oracle_name

  subroutine oracle_domain(this, input_graph, domain, nentries)
    class(affine_residual), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    type(set_graph), intent(out) :: domain
    integer        , intent(out) :: nentries
    integer         :: n_y
    associate (u1 => input_graph); end associate
    domain   = this % y
    nentries = this % n_y
  end subroutine oracle_domain

  subroutine oracle_apply(this, input_graph, input_data, output)

    class(affine_residual), intent(in)             :: this
    class(directed_graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(field)                    :: out
    type(set_graph) :: dom
    real(dp), allocatable          :: q(:), r(:)
    real(dp)                       :: qc, qe
    !----------------------------------------------------------------!
    ! The action's OWN coordinates, compiled here from the identities
    ! and counts it stored. U and Y both enumerate 1..n, so a counted
    ! representation reproduces them exactly and no caller's map has
    ! to be held in this type.
    !----------------------------------------------------------------!

    associate (u1 => input_graph); end associate

    dom = input_data(1) % domain()
    if (.not. dom % same_as(this % u)) then
       error stop 'oracle: the state must live on the unknown domain'
    end if
    call input_data(1) % get_real_vector(q)

    ! Every access through U's own enumeration: U is { e, c }.
    qc = q(this % u_coords % local_index(SLOT_C))
    qe = q(this % u_coords % local_index(SLOT_E))

    allocate(r(2))
    r(this % y_coords % local_index(ROW_C)) = qc - 5.0_dp
    r(this % y_coords % local_index(ROW_E)) = qe - 4.0_dp * qc

    out = field('residual', this % y, this % n_y)
    call out % set_real_vector(r)
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine oracle_apply

end module affine_residual_fixture

program calculator_level_7

  use iso_fortran_env  , only : dp => REAL64
  use calculator_assert, only : report, verdict, SLOT_C, SLOT_E
  use class_graph      , only : directed_stored_graph
  use class_graph_gmres, only : gmres
  use fractal_graph        , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use graph_set_map        , only : set_map
  use graph_set_representation, only : set_representation
  use graph_inclusion_map  , only : inclusion_map, declared_subobject
  use affine_residual_fixture, only : affine_residual, ROW_C, ROW_E

  implicit none

  type(set_graph)     :: x, y, hv
  type(set_graph)      :: u
  type(directed_stored_graph)    :: host
  type(affine_residual) :: oracle
  type(gmres)           :: solver
  type(set_graph) :: dom
  real(dp), allocatable :: g(:), rhs(:), q(:)
  real(dp)              :: achieved
  integer               :: nfail
  type(set_map)     :: sets
  type(inclusion_map)     :: inclusions
  integer         :: n_dom

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "calculator tower . level 7 . minimization"
  write(*,'(1x,a)') "============================================="

  call x % declare()
  call sets % bind(x, counted_set_representation(5))
  call u % declare()
  call sets       % bind(u, listed_set_representation([SLOT_E, SLOT_C]))
  call inclusions % include_in(u, x)
  call y % declare()
  call sets % bind(y, counted_set_representation(2))

  ! Seven vertices: the wrong size for everything, on purpose.
  host = directed_stored_graph(7, tails=[1,2,3,4,5,6], heads=[2,3,4,5,6,7])

  oracle = affine_residual(u, y, sets)

  call report(.not. u % same_as(y), &
       & "U and Y are distinct domains, cardinality notwithstanding", &
       & nfail)
  hv = host % vertex_set()
  call report(.not. hv % same_as(u) .and. &
       &      host % num_vertices() /= sets % size_of(u), &
       & "the host's seven vertices are nobody's unknowns", nfail)

  call solver % attach(oracle, host, u, sets % size_of(u))
  solver % tolerance      = 1.0d-12
  solver % max_iterations = 50

  call solver % domain(host, dom, n_dom)
  call report(dom % same_as(u), &
       & "the solver's answer domain is U, by identity", nfail)
  call oracle % domain(host, dom, n_dom)
  call report(dom % same_as(y), &
       & "and the oracle answers on Y", nfail)

  call solver % constant(g)
  rhs = -g

  ! A deliberately wrong start, stored in U enumeration order {e,c}.
  q = [-7.0_dp, 11.0_dp]
  call solver % solve(rhs, q, achieved)

  call report(achieved < 1.0d-10, &
       & "the residual is driven below tolerance", nfail)
  call report(abs(q(sets % index_in(u, SLOT_C)) - 5.0_dp) < 1.0d-9, &
       & "q(c) = 5, read through U's enumeration", nfail)
  call report(abs(q(sets % index_in(u, SLOT_E)) - 20.0_dp) < 1.0d-9, &
       & "q(e) = 20 - and the host never mattered", nfail)

  call verdict(nfail, "level 7")

end program calculator_level_7
