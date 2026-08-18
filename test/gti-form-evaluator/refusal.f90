!=====================================================================!
! The refusals that must die at the driver seat, one per
! invocation. The specimen is a lying form - it answers the
! abstract contract, and lies:
!
!      siglen        output_signature is [1], not [nentries,ncomp]
!      negent        output_signature claims nentries = -1
!      zerocomp      output_signature claims ncomp = 0
!      valueshape    the declaration says [2,1], value fills one
!      partialshape  the declaration says [2,1], the partial
!                    action fills one
!      storage       the filled shape matches the declaration, but
!                    the storage holds too few scalars
!
! Every case must error stop inside the evaluator's shape law; a
! case that returns is a failure of the suite.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_lying_forms

  use iso_fortran_env      , only : dp => REAL64
  use gti_value_buffers    , only : gti_value_buffer
  use gti_evaluation_points, only : gti_evaluation_point
  use gti_form_interface   , only : gti_differentiable_form, &
       & gti_partial_request, gti_direction_bundle

  implicit none

  private
  public :: lying_form

  !===================================================================!
  ! One configurable liar: it declares `claimed` as its output
  ! signature, fills `value_fills` and `partial_fills` scalars
  ! through the lawful setter, and - when starving - writes the
  ! buffer's fields directly so shape and storage disagree.
  !===================================================================!

  type, extends(gti_differentiable_form) :: lying_form

     integer, allocatable :: claimed(:)
     integer :: value_fills   = 1
     integer :: partial_fills = 1
     logical :: starves_storage = .false.

   contains

     procedure :: name             => lying_name
     procedure :: input_signature  => lying_input_signature
     procedure :: output_signature => lying_output_signature
     procedure :: max_degree       => lying_max_degree
     procedure :: value            => lying_value
     procedure :: partial_action   => lying_partial_action

  end type lying_form

contains

  pure function lying_name(this) result(name)

    class(lying_form), intent(in) :: this
    character(len=:), allocatable :: name

    name = 'lying form'

  end function lying_name

  pure function lying_input_signature(this) result(signature)

    class(lying_form), intent(in) :: this
    integer, allocatable :: signature(:)

    signature = [integer ::]

  end function lying_input_signature

  pure function lying_output_signature(this) result(signature)

    class(lying_form), intent(in) :: this
    integer, allocatable :: signature(:)

    signature = this % claimed

  end function lying_output_signature

  pure function lying_max_degree(this) result(degree)

    class(lying_form), intent(in) :: this
    integer :: degree

    degree = 2

  end function lying_max_degree

  subroutine lying_value(this, point, output)

    class(lying_form)         , intent(in)    :: this
    type(gti_evaluation_point), intent(in)    :: point
    type(gti_value_buffer)    , intent(inout) :: output

    ! the liar reads nothing of its point: its answer is its lie
    associate(unread => point)
    end associate

    if (this % starves_storage) then
       output % nentries = this % claimed(1)
       output % ncomp    = this % claimed(2)
       output % rvals    = [0.0_dp]
    else
       call output % set_real(spread(0.0_dp, dim=1, ncopies=this % value_fills))
    end if

  end subroutine lying_value

  subroutine lying_partial_action(this, point, request, directions, output)

    class(lying_form)         , intent(in)    :: this
    type(gti_evaluation_point), intent(in)    :: point
    type(gti_partial_request) , intent(in)    :: request
    type(gti_direction_bundle), intent(in)    :: directions(:)
    type(gti_value_buffer)    , intent(inout) :: output

    ! the liar reads neither point, request, nor directions
    associate(p => point, r => request, d => directions)
    end associate

    call output % set_real(spread(0.0_dp, dim=1, ncopies=this % partial_fills))

  end subroutine lying_partial_action

end module gti_lying_forms

program refusal

  use gti_value_buffers    , only : gti_value_buffer
  use gti_evaluation_points, only : gti_evaluation_point
  use gti_form_interface   , only : gti_partial_request, gti_direction_bundle
  use gti_form_evaluators  , only : gti_form_evaluator
  use gti_lying_forms      , only : lying_form

  implicit none

  type(gti_form_evaluator)   :: evaluator
  type(lying_form)           :: liar
  type(gti_evaluation_point) :: point
  type(gti_partial_request)  :: request
  type(gti_direction_bundle) :: no_directions(0)
  type(gti_value_buffer)     :: out
  character(len=32) :: which

  call get_command_argument(1, which)

  select case (trim(which))

  case ('siglen')

     liar % claimed = [1]
     call evaluator % value(liar, point, out)

  case ('negent')

     liar % claimed = [-1, 1]
     call evaluator % value(liar, point, out)

  case ('zerocomp')

     liar % claimed = [1, 0]
     call evaluator % value(liar, point, out)

  case ('valueshape')

     liar % claimed = [2, 1]
     liar % value_fills = 1
     call evaluator % value(liar, point, out)

  case ('partialshape')

     liar % claimed = [2, 1]
     liar % partial_fills = 1
     request = gti_partial_request()
     call evaluator % partial_action(liar, point, request, no_directions, out)

  case ('storage')

     liar % claimed = [2, 1]
     liar % starves_storage = .true.
     call evaluator % value(liar, point, out)

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(which)

end program refusal
