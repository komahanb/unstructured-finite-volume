!=====================================================================!
! The higher-order chain-rule assembly suite: one assembler, one
! point, per-argument derivative channels, and the total
! derivatives of the toy triple
!
!      R_i(q, xi) = q_i^2 + xi
!      F(q, xi)   = 1/2 sum q_i^2 + xi
!      M(q, xi)   = xi sum q_i
!
! assembled to degrees 0, 1, 2 through the public form contract
! alone:
!
!      degree 0 :  Phi
!      degree 1 :  D Phi [x^(1)]
!      degree 2 :  D Phi [x^(2)] + D^2 Phi [x^(1), x^(1)]
!
! along the path with
!
!      q^(1) = [1, 0, 2]      xi^(1) = 0.25
!      q^(2) = [.5, 1, 1.5]   xi^(2) = 0.75
!
! M is the witness that ordered pairs are assembled: its symmetric
! mixed partial enters twice, once as (q, xi) and once as (xi, q).
! An absent second derivative reads as zero and the curvature
! terms still apply.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_gti_chain_rule_assembly

  use iso_fortran_env        , only : dp => REAL64
  use fractal_graph          , only : graph
  use class_graph_field      , only : field
  use gti_value_buffers      , only : gti_value_buffer
  use gti_evaluation_points  , only : gti_evaluation_point
  use gti_form_interface     , only : GTI_ARG_STATE, GTI_ARG_DESIGN
  use gti_chain_rule_assemblies, only : gti_chain_rule_assembler, &
       & gti_chain_channel, gti_chain_input
  use gti_toy_forms          , only : toy_residual_form, toy_functional_form, &
       & toy_mixed_form

  implicit none

  type(graph) :: states, designs
  type(field) :: q_field, xi_field

  type(gti_evaluation_point)     :: point
  type(gti_chain_rule_assembler) :: assembler
  type(gti_chain_input)          :: input, sparse

  type(toy_residual_form)   :: r_form
  type(toy_functional_form) :: f_form
  type(toy_mixed_form)      :: m_form

  type(gti_value_buffer) :: out

  real(dp), allocatable :: rv(:)
  integer :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "gti higher-order chain-rule assembly suite"
  write(*,'(1x,a)') "============================================="

  !-------------------------------------------------------------------!
  ! The point: q = [1,2,3] on the state set, xi = [0.5] on the
  ! design set.
  !-------------------------------------------------------------------!

  call states  % declare()
  call designs % declare()

  q_field = field('q', states, 3)
  call q_field % set_real_vector([1.0_dp, 2.0_dp, 3.0_dp])

  xi_field = field('xi', designs, 1)
  call xi_field % set_real_vector([0.5_dp])

  allocate(point % state % component(1))
  point % state % component(1) % value = q_field

  allocate(point % design % component(1))
  point % design % component(1) % value = xi_field

  r_form % nequations = 3

  !-------------------------------------------------------------------!
  ! The full input: a q channel and an xi channel, each carrying
  ! x^(1) and x^(2).
  !-------------------------------------------------------------------!

  allocate(input % channel(2))

  allocate(input % channel(1) % derivative(2))
  call input % channel(1) % derivative(1) % values % set_real([1.0_dp, 0.0_dp, 2.0_dp])
  call input % channel(1) % derivative(2) % values % set_real([0.5_dp, 1.0_dp, 1.5_dp])

  allocate(input % channel(2) % derivative(2))
  input % channel(2) % derivative(1) % argument_kind = GTI_ARG_DESIGN
  call input % channel(2) % derivative(1) % values % set_real([0.25_dp])
  input % channel(2) % derivative(2) % argument_kind = GTI_ARG_DESIGN
  call input % channel(2) % derivative(2) % values % set_real([0.75_dp])

  !-------------------------------------------------------------------!
  ! The sparse input: same first derivatives, no second anywhere -
  ! the q channel by an empty seat, the xi channel by no seat.
  !-------------------------------------------------------------------!

  allocate(sparse % channel(2))

  allocate(sparse % channel(1) % derivative(2))
  call sparse % channel(1) % derivative(1) % values % set_real([1.0_dp, 0.0_dp, 2.0_dp])

  allocate(sparse % channel(2) % derivative(1))
  sparse % channel(2) % derivative(1) % argument_kind = GTI_ARG_DESIGN
  call sparse % channel(2) % derivative(1) % values % set_real([0.25_dp])

  !-------------------------------------------------------------------!
  ! The input structure answers shape and occupancy separately.
  !-------------------------------------------------------------------!

  call report(input % size() == 2 .and. sparse % size() == 2, &
       & "the input counts its argument channels", nfail)

  call report(input % channel(1) % has_degree(1) .and. &
       &      input % channel(1) % has_degree(2) .and. &
       &      input % channel(2) % has_degree(2), &
       & "the full channels carry both path derivatives", nfail)

  call report(sparse % channel(1) % has_degree(1) .and. &
       & .not. sparse % channel(1) % has_degree(2) .and. &
       & .not. sparse % channel(2) % has_degree(2), &
       & "an empty seat and a missing seat both read as absent", nfail)

  !-------------------------------------------------------------------!
  ! Degree 0: the value itself.
  !-------------------------------------------------------------------!

  call assembler % assemble(r_form, point, 0, input, out)
  call out % get_real(rv)
  call report(matches(rv, [1.5_dp, 4.5_dp, 9.5_dp]), &
       & "R degree 0: R_i = q_i^2 + xi", nfail)

  call assembler % assemble(f_form, point, 0, input, out)
  call out % get_real(rv)
  call report(matches(rv, [7.5_dp]), &
       & "F degree 0: F = 1/2 sum q^2 + xi", nfail)

  call assembler % assemble(m_form, point, 0, input, out)
  call out % get_real(rv)
  call report(matches(rv, [3.0_dp]), &
       & "M degree 0: M = xi sum q", nfail)

  !-------------------------------------------------------------------!
  ! Degree 1: D Phi [x^(1)] summed over the channels.
  !-------------------------------------------------------------------!

  call assembler % assemble(r_form, point, 1, input, out)
  call out % get_real(rv)
  call report(matches(rv, [2.25_dp, 0.25_dp, 12.25_dp]), &
       & "R degree 1: Dq R [q1] + Dxi R [xi1]", nfail)

  call assembler % assemble(f_form, point, 1, input, out)
  call out % get_real(rv)
  call report(matches(rv, [7.25_dp]), &
       & "F degree 1: Dq F [q1] + Dxi F [xi1]", nfail)

  call assembler % assemble(m_form, point, 1, input, out)
  call out % get_real(rv)
  call report(matches(rv, [3.0_dp]), &
       & "M degree 1: xi sum q1 + xi1 sum q", nfail)

  !-------------------------------------------------------------------!
  ! Degree 2: transport of x^(2) plus curvature over ordered
  ! pairs of x^(1).
  !-------------------------------------------------------------------!

  call assembler % assemble(r_form, point, 2, input, out)
  call out % get_real(rv)
  call report(matches(rv, [3.75_dp, 4.75_dp, 17.75_dp]), &
       & "R degree 2: Dq R [q2] + Dxi R [xi2] + D2 R [q1,q1]", nfail)

  call assembler % assemble(f_form, point, 2, input, out)
  call out % get_real(rv)
  call report(matches(rv, [12.75_dp]), &
       & "F degree 2: Dq F [q2] + Dxi F [xi2] + D2 F [q1,q1]", nfail)

  call assembler % assemble(m_form, point, 2, input, out)
  call out % get_real(rv)
  call report(matches(rv, [7.5_dp]), &
       & "M degree 2 counts both ordered cross terms (q,xi) and (xi,q)", nfail)

  !-------------------------------------------------------------------!
  ! Missing second derivatives read as zero; the curvature terms
  ! still apply.
  !-------------------------------------------------------------------!

  call assembler % assemble(r_form, point, 2, sparse, out)
  call out % get_real(rv)
  call report(matches(rv, [2.0_dp, 0.0_dp, 8.0_dp]), &
       & "absent x2 is zero: R degree 2 keeps only D2 R [q1,q1]", nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all chain-rule assembly checks passed"
  else
     error stop
  end if

contains

  pure function matches(values, expected) result(ok)

    real(dp), intent(in) :: values(:), expected(:)
    logical :: ok

    ok = size(values) == size(expected)
    if (ok) ok = all(abs(values - expected) < 1.0e-14_dp)

  end function matches

  subroutine report(ok, label, nfail)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    integer, intent(inout) :: nfail
    if (ok) then
       write(*,'(1x,a,a)') "PASS : ", label
    else
       write(*,'(1x,a,a)') "FAIL : ", label
       nfail = nfail + 1
    end if
  end subroutine report

end program test_gti_chain_rule_assembly
