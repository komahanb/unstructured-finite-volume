!=====================================================================!
! PROTOTYPE . SEMANTIC STORAGE AT SCALE
!
! Measures the size of one kernel graph object and extrapolates the
! two candidate designs to the scales the architecture must reach.
! The numbers are printed, not asserted: the assertion is only that
! design A crosses a stated memory limit and design B does not.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program scale

  use graph_fractal, only : graph, branch

  implicit none

  real, parameter :: GB = 1.0e9
  real, parameter :: MEMORY_LIMIT_GB = 64.0        ! one compute node's memory capacity

  integer :: gbytes, bbytes, ibytes
  real    :: a_set_1e9, b_set_1e9, a_nnz_1e10, b_nnz_1e10
  integer :: failures = 0

  block
    type(graph)        :: g
    type(branch) :: b
    integer            :: i
    gbytes = storage_size(g) / 8
    bbytes = storage_size(b) / 8
    ibytes = storage_size(i) / 8
  end block

  write(*,'(1x,a)') "semantic storage at scale"
  print '(1x,a,i0,a,i0,a,i0,a)', ' one graph = ', gbytes, &
       & ' bytes (branch ', bbytes, ', integer ', ibytes, ')'
  print *, ''

  ! A set of 10^9 members.
  !   A: one element graph and one sequence cell per member.
  !   B: one graph for the set, plus one integer for the count.
  a_set_1e9 = 2.0 * 1.0e9 * gbytes / GB
  b_set_1e9 = (gbytes + ibytes) / GB

  ! A relation with 10^10 tuples of arity 2.
  !   A: one tuple-holder graph and two component cells per tuple.
  !   B: one graph for the relation, plus CSR - nnz targets and n+1
  !      offsets, both directions, as flat integer arrays.
  a_nnz_1e10 = 3.0 * 1.0e10 * gbytes / GB
  b_nnz_1e10 = (2.0 * 1.0e10 + 2.0 * 1.0e9) * ibytes / GB

  print '(1x,a)', ' |S| = 10^9'
  print '(1x,a,f12.1,a)', '   A  extensional graph      : ', a_set_1e9,  ' GB of graph objects'
  print '(1x,a,f12.6,a)', '   B  semantic graph + extent: ', b_set_1e9,  ' GB'
  print *, ''
  print '(1x,a)', ' |R| = 10^10 tuples, arity 2'
  print '(1x,a,f12.1,a)', '   A  extensional graph      : ', a_nnz_1e10, ' GB of graph objects'
  print '(1x,a,f12.1,a)', '   B  semantic graph + CSR   : ', b_nnz_1e10, ' GB of flat integer arrays'
  print *, ''

  call check('design A exceeds the stated memory limit for a 10^9 set', &
       & a_set_1e9 .gt. MEMORY_LIMIT_GB)
  call check('design B does not, by nine orders of magnitude', &
       & b_set_1e9 .lt. 1.0e-6 * MEMORY_LIMIT_GB)
  call check('design A exceeds it for 10^10 tuples', &
       & a_nnz_1e10 .gt. MEMORY_LIMIT_GB)
  call check('design B stores tuples as flat arrays, no pointer per nonzero', &
       & b_nnz_1e10 .lt. a_nnz_1e10 / 10.0)

  if (failures .eq. 0) then
     print *, ''
     print *, ' ALL PROPOSITIONS HOLD'
  else
     print *, ''
     print *, ' FAILURES :', failures
     error stop 'scale: a proposition failed'
  end if

contains

  subroutine check(label, satisfied)

    character(len=*), intent(in) :: label
    logical         , intent(in) :: satisfied

    if (satisfied) then
       print *, ' PASS : ', label
    else
       print *, ' FAIL : ', label
       failures = failures + 1
    end if

  end subroutine check

end program scale
