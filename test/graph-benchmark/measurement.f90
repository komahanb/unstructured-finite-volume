!=====================================================================!
! One timer, one allocation counter and one record line for every
! scaling benchmark. A record is one line of key=value tokens, read by
! scaling.py; a phase measures wall seconds, allocation calls and
! requested bytes between begin and end.
!=====================================================================!

module benchmark_measurement

  use iso_fortran_env, only : int64
  use iso_c_binding  , only : c_int64_t
  use util_precision , only : dp

  implicit none

  private
  public :: phase, begin_phase, end_phase, record, peak_rss_kilobytes
  public :: argument_integer, argument_string, argument_real, pseudo_random_permutation

  interface
     subroutine allocation_begin() bind(c, name='allocation_begin')
     end subroutine allocation_begin
     subroutine allocation_end(calls, bytes) bind(c, name='allocation_end')
       import c_int64_t
       integer(c_int64_t), intent(out) :: calls, bytes
     end subroutine allocation_end
     function peak_rss_kilobytes() bind(c, name='peak_rss_kilobytes') result(kilobytes)
       import c_int64_t
       integer(c_int64_t) :: kilobytes
     end function peak_rss_kilobytes
  end interface

  type :: phase
     integer(int64) :: started = 0
     real(dp) :: seconds = 0.0_dp
     integer(c_int64_t) :: calls = 0, bytes = 0
  end type phase

contains

  subroutine begin_phase(this)
    type(phase), intent(out) :: this
    integer(int64) :: rate
    call allocation_begin()
    call system_clock(this % started, rate)
  end subroutine begin_phase

  subroutine end_phase(this)
    type(phase), intent(inout) :: this
    integer(int64) :: finished, rate
    call system_clock(finished, rate)
    call allocation_end(this % calls, this % bytes)
    this % seconds = real(finished - this % started, dp) / real(rate, dp)
  end subroutine end_phase

  ! One record line: the case tokens, then the phase name and its measures.
  subroutine record(case_tokens, name, this)
    character(len=*), intent(in) :: case_tokens, name
    type(phase), intent(in) :: this
    write(*, '(a,1x,a,1x,a,a,1x,a,es16.8e3,1x,a,i0,1x,a,i0)') 'record', trim(case_tokens), &
         & 'phase=', trim(name), 'seconds=', this % seconds, 'calls=', this % calls, 'bytes=', this % bytes
  end subroutine record

  integer function argument_integer(position, default) result(value)
    integer, intent(in) :: position, default
    character(len=64) :: argument
    value = default
    if (command_argument_count() < position) return
    call get_command_argument(position, argument)
    if (len_trim(argument) > 0) read(argument, *) value
  end function argument_integer

  real(dp) function argument_real(position, default) result(value)
    integer, intent(in) :: position
    real(dp), intent(in) :: default
    character(len=64) :: argument
    value = default
    if (command_argument_count() < position) return
    call get_command_argument(position, argument)
    if (len_trim(argument) > 0) read(argument, *) value
  end function argument_real

  function argument_string(position, default) result(value)
    integer, intent(in) :: position
    character(len=*), intent(in) :: default
    character(len=:), allocatable :: value
    character(len=64) :: argument
    value = default
    if (command_argument_count() < position) return
    call get_command_argument(position, argument)
    if (len_trim(argument) > 0) value = trim(argument)
  end function argument_string

  ! A reproducible permutation of 1..n from a linear congruential
  ! sequence with a fixed seed: the same input on every run.
  function pseudo_random_permutation(n, seed) result(values)
    integer, intent(in) :: n, seed
    integer, allocatable :: values(:)
    integer(int64) :: state
    integer :: i, j, swapped
    values = [(i, i = 1, n)]
    state = int(max(seed, 1), int64)
    do i = n, 2, -1
       state = mod(48271_int64 * state, 2147483647_int64)
       j = int(mod(state, int(i, int64)), kind(j)) + 1
       swapped = values(i)
       values(i) = values(j)
       values(j) = swapped
    end do
  end function pseudo_random_permutation

end module benchmark_measurement
