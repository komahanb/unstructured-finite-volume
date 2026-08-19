!=====================================================================!
! The refusals that must die at the dense direct seat, one per
! invocation:
!
!      tolzero       a zero singular tolerance
!      sizemismatch  a solution array whose length disagrees with
!                    the right-hand side
!      singular      dependent rows - no pivot survives elimination
!      nonsquare     a rectangular array at the adapter door
!      probewidth    a probed width carrying a fractional number
!                    per member of the operation's domain
!
! Every case must error stop; a case that returns is a failure of
! the suite.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program refusal

  use iso_fortran_env     , only : dp => REAL64
  use class_graph_stencil , only : stencil_operator
  use class_graph_dense_direct, only : dense_direct, &
       & solve_dense_matrix_with_dense_direct, probed_dense_matrix

  implicit none

  type(stencil_operator) :: statement
  type(dense_direct)     :: solver

  real(dp), allocatable :: xa(:), aprobe(:,:)
  real(dp) :: x2(2), x3(3), achieved
  real(dp) :: rect(2,3)

  character(len=32) :: which

  call get_command_argument(1, which)

  select case (trim(which))

  case ('tolzero')

     statement = stencil_operator([1], [1], [2.0_dp], [0.0_dp], 'one')
     call solver % attach(statement, statement % pattern, &
          & statement % pattern % vertex_set(), &
          & statement % pattern % num_vertices())
     solver % singular_tolerance = 0.0_dp
     x2(1:1) = 0.0_dp
     call solver % solve([6.0_dp], x2(1:1), achieved)

  case ('sizemismatch')

     statement = stencil_operator([1, 1, 2, 2], [1, 2, 1, 2], &
          & [2.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [0.0_dp, 0.0_dp], 'two')
     call solver % attach(statement, statement % pattern, &
          & statement % pattern % vertex_set(), &
          & statement % pattern % num_vertices())
     x3 = 0.0_dp
     call solver % solve([4.0_dp, 7.0_dp], x3, achieved)

  case ('singular')

     ! row 2 = 2 * row 1: no pivot survives elimination
     statement = stencil_operator([1, 1, 2, 2], [1, 2, 1, 2], &
          & [1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp], [0.0_dp, 0.0_dp], 'flat')
     call solver % attach(statement, statement % pattern, &
          & statement % pattern % vertex_set(), &
          & statement % pattern % num_vertices())
     x2 = 0.0_dp
     call solver % solve([1.0_dp, 2.0_dp], x2, achieved)

  case ('nonsquare')

     rect = 1.0_dp
     call solve_dense_matrix_with_dense_direct(rect, [1.0_dp, 1.0_dp], &
          & 1.0e-14_dp, xa, achieved)

  case ('probewidth')

     ! three numbers over a two-member domain is not a whole
     ! number per member
     statement = stencil_operator([1, 1, 2, 2], [1, 2, 1, 2], &
          & [2.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [0.0_dp, 0.0_dp], 'two')
     call probed_dense_matrix(statement, statement % pattern, 3, aprobe)

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(which)

end program refusal
