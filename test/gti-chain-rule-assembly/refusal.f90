!=====================================================================!
! The refusals that must die at the assembly seat, one per
! invocation:
!
!      negdeg    a derivative degree of -1
!      highdeg   a derivative degree of 3, past phase-2 support
!      dupq      two channels both naming the state component q
!      dupxi     two channels both naming the design xi
!      shifty    a form whose output_signature answer shifts
!                between calls, so a term that passed the
!                evaluator's shape law no longer matches the
!                running accumulation
!
! The shifting form is the witness for the accumulation guard:
! every term is held to the form's declaration by the evaluator,
! but the declaration is re-read per call - only the accumulation
! check pins the shape across calls.
!
! Every case must error stop before any wrong number is produced;
! a case that returns is a failure of the suite.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_shifting_forms

  use iso_fortran_env      , only : dp => REAL64
  use gti_value_buffers    , only : gti_value_buffer
  use gti_evaluation_points, only : gti_evaluation_point
  use gti_form_interface   , only : gti_differentiable_form, &
       & gti_partial_request, gti_direction_bundle, GTI_ARG_STATE

  implicit none

  private
  public :: shifting_form

  !===================================================================!
  ! The module counter the pure signature reads and the impure
  ! action advances: the declaration is [2,1] until the first
  ! partial action runs, then [1,1].
  !===================================================================!

  integer :: shift_count = 0

  type, extends(gti_differentiable_form) :: shifting_form

   contains

     procedure :: name             => shifting_name
     procedure :: input_signature  => shifting_input_signature
     procedure :: output_signature => shifting_output_signature
     procedure :: max_degree       => shifting_max_degree
     procedure :: value            => shifting_value
     procedure :: partial_action   => shifting_partial_action

  end type shifting_form

contains

  pure function shifting_name(this) result(name)

    class(shifting_form), intent(in) :: this
    character(len=:), allocatable :: name

    name = 'shifting form'

  end function shifting_name

  pure function shifting_input_signature(this) result(signature)

    class(shifting_form), intent(in) :: this
    integer, allocatable :: signature(:)

    signature = [GTI_ARG_STATE]

  end function shifting_input_signature

  pure function shifting_output_signature(this) result(signature)

    class(shifting_form), intent(in) :: this
    integer, allocatable :: signature(:)

    if (shift_count == 0) then
       signature = [2, 1]
    else
       signature = [1, 1]
    end if

  end function shifting_output_signature

  pure function shifting_max_degree(this) result(degree)

    class(shifting_form), intent(in) :: this
    integer :: degree

    degree = 2

  end function shifting_max_degree

  subroutine shifting_value(this, point, output)

    class(shifting_form)      , intent(in)    :: this
    type(gti_evaluation_point), intent(in)    :: point
    type(gti_value_buffer)    , intent(inout) :: output

    ! the shifter reads nothing of its point
    associate(unread => point)
    end associate

    call output % set_real([0.0_dp, 0.0_dp])

  end subroutine shifting_value

  subroutine shifting_partial_action(this, point, request, directions, output)

    class(shifting_form)      , intent(in)    :: this
    type(gti_evaluation_point), intent(in)    :: point
    type(gti_partial_request) , intent(in)    :: request
    type(gti_direction_bundle), intent(in)    :: directions(:)
    type(gti_value_buffer)    , intent(inout) :: output

    ! the shifter reads neither point, request, nor directions
    associate(p => point, r => request, d => directions)
    end associate

    shift_count = shift_count + 1
    call output % set_real([0.0_dp])

  end subroutine shifting_partial_action

end module gti_shifting_forms

program refusal

  use iso_fortran_env        , only : dp => REAL64
  use gti_value_buffers      , only : gti_value_buffer
  use gti_evaluation_points  , only : gti_evaluation_point
  use gti_form_interface     , only : GTI_ARG_DESIGN
  use gti_chain_rule_assemblies, only : gti_chain_rule_assembler, gti_chain_input
  use gti_toy_forms          , only : toy_residual_form
  use gti_shifting_forms     , only : shifting_form

  implicit none

  type(gti_chain_rule_assembler) :: assembler
  type(gti_evaluation_point)     :: point
  type(gti_chain_input)          :: input
  type(toy_residual_form)        :: r_form
  type(shifting_form)            :: shifter
  type(gti_value_buffer)         :: out
  character(len=32) :: which

  call get_command_argument(1, which)

  select case (trim(which))

  case ('negdeg')

     call assembler % assemble(r_form, point, -1, input, out)

  case ('highdeg')

     call assembler % assemble(r_form, point, 3, input, out)

  case ('dupq')

     allocate(input % channel(2))
     allocate(input % channel(1) % derivative(1))
     call input % channel(1) % derivative(1) % values % set_real([1.0_dp, 0.0_dp, 2.0_dp])
     allocate(input % channel(2) % derivative(1))
     call input % channel(2) % derivative(1) % values % set_real([1.0_dp, 0.0_dp, 2.0_dp])
     call assembler % assemble(r_form, point, 1, input, out)

  case ('dupxi')

     allocate(input % channel(2))
     allocate(input % channel(1) % derivative(1))
     input % channel(1) % derivative(1) % argument_kind = GTI_ARG_DESIGN
     call input % channel(1) % derivative(1) % values % set_real([0.25_dp])
     allocate(input % channel(2) % derivative(1))
     input % channel(2) % derivative(1) % argument_kind = GTI_ARG_DESIGN
     call input % channel(2) % derivative(1) % values % set_real([0.75_dp])
     call assembler % assemble(r_form, point, 1, input, out)

  case ('shifty')

     allocate(input % channel(1))
     allocate(input % channel(1) % derivative(1))
     call input % channel(1) % derivative(1) % values % set_real([1.0_dp, 0.0_dp, 2.0_dp])
     call assembler % assemble(shifter, point, 1, input, out)

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(which)

end program refusal
