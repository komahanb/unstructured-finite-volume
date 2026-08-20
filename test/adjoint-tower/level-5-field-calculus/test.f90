!=====================================================================!
! ADJOINT TOWER . LEVEL 5 . THE FIRST NUMBERS
!
! The level answers one question: WHAT VALUES EXIST BEFORE ANYTHING
! IS SOLVED. Two fields, and no more:
!
!      p  = [2]      on P = { p }        the parameter, given
!      q0 = [0, 0]   on Q = { u, v }     the state, deliberately WRONG
!
! The initial state is not an answer and does not pretend to be one:
! the true state is [2, 4] and Level 7 must travel to it. What is
! NOT here matters as much as what is. There is no residual field,
! no response field, and above all
!
!      no lambda = 0
!
! merely because Level 7 will one day need an adjoint. A zero
! adjoint is not an empty adjoint - it is a claim about a solution
! that has not been computed, and fabricating one would let a
! later bug pass by never overwriting it. An initial guess belongs
! to the solve that makes it, not to the rung that declares
! domains.
!
! Roles remain domains: both fields are the ordinary production
! field, and only their domains say which is a parameter and which
! is a state. No parameter_field, no state_field, no adjoint_field.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program adjoint_level_5

  use iso_fortran_env  , only : dp => REAL64
  use adjoint_assert   , only : report, verdict
  use adjoint_assert   , only : VAR_P, VAR_U, VAR_V
  use graph_fractal        , only : set_graph => graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map
  use map_inclusion  , only : inclusion_map, declared_subobject
  use class_graph_field, only : field

  implicit none

  type(set_graph) :: v
  type(set_graph)  :: p_dom, q_dom
  type(field)       :: p_field, q0
  integer           :: nfail
  type(set_map)     :: sets
  type(inclusion_map)     :: inclusions

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "adjoint tower . level 5 . fields"
  write(*,'(1x,a)') "============================================="

  call v % declare()
  call sets % bind(v, counted_set_representation(3))
  call p_dom % declare()
  call sets       % bind(p_dom, listed_set_representation([VAR_P]))
  call inclusions % include_in(p_dom, v)
  call q_dom % declare()
  call sets       % bind(q_dom, listed_set_representation([VAR_U, VAR_V]))
  call inclusions % include_in(q_dom, v)

  call check_parameter(nfail)
  call check_initial_state(nfail)
  call check_roles_by_domain(nfail)

  ! No residual field, no response field, no adjoint field is built
  ! here - and no zero stands in for any of them.

  call verdict(nfail, "level 5")

contains

  !===================================================================!
  ! The parameter: one entry, read through P's own enumeration.
  !===================================================================!

  subroutine check_parameter(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: dom
    real(dp), allocatable          :: val(:)

    p_field = field('parameter', p_dom, sets % size_of(p_dom))
    call p_field % set_real_vector([2.0_dp])

    dom = p_field % domain()
    call report(dom % same_as(p_dom), &
         & "the parameter field's domain is P, by identity", nfail)
    call report(p_field % num_entries() .eq. 1, &
         & "one parameter entry", nfail)

    call p_field % get_real_vector(val)
    call report(abs(val(sets % index_in(p_dom, VAR_P)) - 2.0_dp) &
         &      < 1.0d-14, &
         & "p = 2, read through the domain map", nfail)

  end subroutine check_parameter

  !===================================================================!
  ! The initial state: [0, 0] on { u, v }, in declaration order -
  ! and wrong on purpose. Nothing here knows the true [2, 4].
  !===================================================================!

  subroutine check_initial_state(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: dom
    real(dp), allocatable          :: val(:)

    q0 = field('initial state', q_dom, sets % size_of(q_dom))
    call q0 % set_real_vector([0.0_dp, 0.0_dp])

    dom = q0 % domain()
    call report(dom % same_as(q_dom), &
         & "the state field's domain is Q, by identity", nfail)
    call report(q0 % num_entries() .eq. 2 .and. &
         &      q0 % num_components() .eq. 1, &
         & "two state entries, one component", nfail)

    call q0 % get_real_vector(val)
    call report(abs(val(sets % index_in(q_dom, VAR_U))) < 1.0d-14 .and. &
         &      abs(val(sets % index_in(q_dom, VAR_V))) < 1.0d-14, &
         & "q0 = [0, 0]: a starting point, not an answer", nfail)

    call report(sets % index_in(q_dom, VAR_U) .eq. 1 .and. &
         &      sets % index_in(q_dom, VAR_V) .eq. 2 .and. &
         &      sets % index_in(v, VAR_V) .eq. 3, &
         & "storage follows Q's enumeration, not the variable ids", &
         & nfail)

  end subroutine check_initial_state

  !===================================================================!
  ! Parameter and state are told apart by their DOMAINS - one field
  ! type, two roles, no subclass anywhere.
  !===================================================================!

  subroutine check_roles_by_domain(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: dp_, dq

    dp_ = p_field % domain()
    dq = q0 % domain()

    call report(.not. dp_ % same_as(q_dom) .and. &
         &      .not. dq % same_as(p_dom), &
         & "parameter and state are distinguished by domain alone", &
         & nfail)
    call report(sets % size_of(dp_) .ne. sets % size_of(dq), &
         & "and here even their sizes differ - though sizes are " // &
         & "never what settles it", nfail)

  end subroutine check_roles_by_domain

end program adjoint_level_5
