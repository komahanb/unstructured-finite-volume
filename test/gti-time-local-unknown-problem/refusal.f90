!=====================================================================!
! The refusals that must die at the unknown seat, one per
! invocation:
!
!      nosamples   an empty sample set
!      idx0        unknown index 0, below the samples
!      idxhigh     unknown index past the samples
!      missingq    an unknown sample whose bundle holds no q
!      novalues    a trial buffer holding no real values
!      entries     a trial with the wrong entry count
!      comp        a trial with the wrong component count
!      storage     a trial whose storage contradicts its own shape
!      badfield    an unknown q of a foreign dynamic type
!
! The foreign field is a test-local mock: it answers the abstract
! field contract with the right shape, so validation admits it,
! and the mint refuses it - the one place a dynamic type is
! trusted is the one place it is checked.
!
! Every case must error stop before any residual is evaluated; a
! case that returns is a failure of the suite.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_mock_fields

  use iso_fortran_env     , only : dp => REAL64
  use graph_field_calculus, only : graph_field, set_graph, GRAPH_FIELD_REAL

  implicit none

  private
  public :: mock_field

  !===================================================================!
  ! A field of a dynamic type the unknown seat does not support:
  ! it claims shape (3, 1) and values [1, 2, 4], enough to pass
  ! every shape law, and is not class_graph_field::field.
  !===================================================================!

  type, extends(graph_field) :: mock_field

   contains

     procedure :: name   => mock_name
     procedure :: units  => mock_units
     procedure :: domain => mock_domain

     procedure :: num_components => mock_num_components
     procedure :: num_entries    => mock_num_entries
     procedure :: value_kind     => mock_value_kind

     procedure :: get_integer_vector   => mock_get_integer
     procedure :: set_integer_vector   => mock_set_integer
     procedure :: get_real_vector      => mock_get_real
     procedure :: set_real_vector      => mock_set_real
     procedure :: get_complex_vector   => mock_get_complex
     procedure :: set_complex_vector   => mock_set_complex
     procedure :: get_logical_vector   => mock_get_logical
     procedure :: set_logical_vector   => mock_set_logical
     procedure :: get_character_vector => mock_get_character
     procedure :: set_character_vector => mock_set_character

  end type mock_field

contains

  pure function mock_name(this) result(name)
    class(mock_field), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'mock'
  end function mock_name

  pure function mock_units(this) result(units)
    class(mock_field), intent(in) :: this
    character(len=:), allocatable :: units
    units = '-'
  end function mock_units

  type(set_graph) function mock_domain(this) result(domain)
    class(mock_field), intent(in) :: this
    ! an undeclared domain: never asked in this suite
  end function mock_domain

  pure integer function mock_num_components(this)
    class(mock_field), intent(in) :: this
    mock_num_components = 1
  end function mock_num_components

  pure integer function mock_num_entries(this)
    class(mock_field), intent(in) :: this
    mock_num_entries = 3
  end function mock_num_entries

  pure integer function mock_value_kind(this)
    class(mock_field), intent(in) :: this
    mock_value_kind = GRAPH_FIELD_REAL
  end function mock_value_kind

  pure subroutine mock_get_integer(this, values)
    class(mock_field), intent(in) :: this
    integer, allocatable, intent(out) :: values(:)
    allocate(values(0))
  end subroutine mock_get_integer

  pure subroutine mock_set_integer(this, values)
    class(mock_field), intent(inout) :: this
    integer, intent(in) :: values(:)
    ! the mock keeps nothing: the values are ignored by design
    associate(unkept => values)
    end associate
  end subroutine mock_set_integer

  pure subroutine mock_get_real(this, values)
    class(mock_field), intent(in) :: this
    real(dp), allocatable, intent(out) :: values(:)
    values = [1.0_dp, 2.0_dp, 4.0_dp]
  end subroutine mock_get_real

  pure subroutine mock_set_real(this, values)
    class(mock_field), intent(inout) :: this
    real(dp), intent(in) :: values(:)
    ! the mock keeps nothing: the values are ignored by design
    associate(unkept => values)
    end associate
  end subroutine mock_set_real

  pure subroutine mock_get_complex(this, values)
    class(mock_field), intent(in) :: this
    complex(dp), allocatable, intent(out) :: values(:)
    allocate(values(0))
  end subroutine mock_get_complex

  pure subroutine mock_set_complex(this, values)
    class(mock_field), intent(inout) :: this
    complex(dp), intent(in) :: values(:)
    ! the mock keeps nothing: the values are ignored by design
    associate(unkept => values)
    end associate
  end subroutine mock_set_complex

  pure subroutine mock_get_logical(this, values)
    class(mock_field), intent(in) :: this
    logical, allocatable, intent(out) :: values(:)
    allocate(values(0))
  end subroutine mock_get_logical

  pure subroutine mock_set_logical(this, values)
    class(mock_field), intent(inout) :: this
    logical, intent(in) :: values(:)
    ! the mock keeps nothing: the values are ignored by design
    associate(unkept => values)
    end associate
  end subroutine mock_set_logical

  pure subroutine mock_get_character(this, values)
    class(mock_field), intent(in) :: this
    character(len=:), allocatable, intent(out) :: values(:)
    allocate(character(len=1) :: values(0))
  end subroutine mock_get_character

  pure subroutine mock_set_character(this, values)
    class(mock_field), intent(inout) :: this
    character(len=*), intent(in) :: values(:)
    ! the mock keeps nothing: the values are ignored by design
    associate(unkept => values)
    end associate
  end subroutine mock_set_character

end module gti_mock_fields

program refusal

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use class_graph_field    , only : field
  use gti_value_buffers    , only : gti_value_buffer
  use gti_time_local_schemes, only : gti_time_sample
  use gti_time_local_unknown_problems, only : gti_time_local_unknown_problem
  use gti_mock_fields      , only : mock_field

  implicit none

  type(graph) :: states
  type(field) :: q0_field, q1_field
  type(mock_field) :: foreign

  type(gti_time_sample) :: samples(2)
  type(gti_time_sample) :: none(0)
  type(gti_time_sample), allocatable :: work(:)

  type(gti_time_local_unknown_problem) :: problem
  type(gti_value_buffer) :: trial
  character(len=32) :: which

  call get_command_argument(1, which)

  call states % declare()

  q0_field = field('q0', states, 3)
  call q0_field % set_real_vector([1.0_dp, 2.0_dp, 4.0_dp])
  q1_field = field('q1', states, 3)
  call q1_field % set_real_vector([3.0_dp, 5.0_dp, 8.0_dp])

  allocate(samples(1) % state % component(1))
  samples(1) % state % component(1) % value = q0_field
  allocate(samples(2) % state % component(1))
  samples(2) % state % component(1) % value = q1_field

  call trial % set_real([3.0_dp, 5.0_dp, 8.0_dp])

  select case (trim(which))

  case ('nosamples')

     call problem % inject_trial_q(none, 1, trial, work)

  case ('idx0')

     call problem % inject_trial_q(samples, 0, trial, work)

  case ('idxhigh')

     call problem % inject_trial_q(samples, 3, trial, work)

  case ('missingq')

     deallocate(samples(2) % state % component)
     call problem % inject_trial_q(samples, 2, trial, work)

  case ('novalues')

     call trial % clear()
     call problem % inject_trial_q(samples, 2, trial, work)

  case ('entries')

     call trial % set_real([3.0_dp, 5.0_dp])
     call problem % inject_trial_q(samples, 2, trial, work)

  case ('comp')

     call trial % set_real([3.0_dp, 5.0_dp, 8.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], ncomp=2)
     call problem % inject_trial_q(samples, 2, trial, work)

  case ('storage')

     trial % nentries = 3
     trial % ncomp    = 1
     trial % rvals    = [3.0_dp, 5.0_dp]
     call problem % inject_trial_q(samples, 2, trial, work)

  case ('badfield')

     samples(2) % state % component(1) % value = foreign
     call problem % inject_trial_q(samples, 2, trial, work)

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(which)

end program refusal
