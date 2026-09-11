!=====================================================================!
! Invalid-input cases for the dense direct solver, one per
! invocation. The case is selected by the first command-line
! argument; each must terminate in error stop with the message
! run.sh expects, and a case that returns normally is reported as
! a failure by run.sh.
!
!      zero_tolerance       a zero singular tolerance
!      size_mismatch  a solution array whose length disagrees with
!                    the right-hand side
!      singular      dependent rows - no pivot survives elimination
!      nonsquare     a rectangular array laid on a stencil
!      nonintegral_width      a compiled-stencil width representing a fractional
!                    number per member of the operation's domain
!      incompatible_components     compiling a stencil at a width whose values per
!                    member violate the operation's argument contract
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program refusal

  use iso_fortran_env     , only : dp => REAL64
  use operation_stencil , only : stencil
  use operation_dense_direct, only : dense_direct

  implicit none

  type(stencil) :: statement, compiled
  type(dense_direct)     :: solver

  real(dp) :: x2(2), x3(3), achieved
  real(dp) :: rectangular_matrix(2,3)

  character(len=32) :: case_name

  call get_command_argument(1, case_name)

  select case (trim(case_name))

  case ('zero_tolerance')

     statement = stencil([1], [1], [2.0_dp], [0.0_dp], 'one')
     call solver % state(statement, statement % pattern, &
          & statement % pattern % vertex_set(), &
          & statement % pattern % num_vertices())
     solver % singular_tolerance = 0.0_dp
     x2(1:1) = 0.0_dp
     call solver % solve([6.0_dp], x2(1:1), achieved)

  case ('size_mismatch')

     statement = stencil([1, 1, 2, 2], [1, 2, 1, 2], &
          & [2.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [0.0_dp, 0.0_dp], 'two')
     call solver % state(statement, statement % pattern, &
          & statement % pattern % vertex_set(), &
          & statement % pattern % num_vertices())
     x3 = 0.0_dp
     call solver % solve([4.0_dp, 7.0_dp], x3, achieved)

  case ('singular')

     ! row 2 = 2 * row 1, so elimination produces no usable pivot
     statement = stencil([1, 1, 2, 2], [1, 2, 1, 2], &
          & [1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp], [0.0_dp, 0.0_dp], 'singular')
     call solver % state(statement, statement % pattern, &
          & statement % pattern % vertex_set(), &
          & statement % pattern % num_vertices())
     x2 = 0.0_dp
     call solver % solve([1.0_dp, 2.0_dp], x2, achieved)

  case ('nonsquare')

     rectangular_matrix = 1.0_dp
     statement = stencil(rectangular_matrix, 'rect')

  case ('nonintegral_width')

     ! three numbers over a two-member domain is not a whole
     ! number per member
     statement = stencil([1, 1, 2, 2], [1, 2, 1, 2], &
          & [2.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [0.0_dp, 0.0_dp], 'two')
     compiled = stencil(statement, statement % pattern, 3)

  case ('incompatible_components')

     ! four numbers over a two-member domain is two per member, and
     ! the stencil's argument contract admits one
     statement = stencil([1, 1, 2, 2], [1, 2, 1, 2], &
          & [2.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [0.0_dp, 0.0_dp], 'two')
     compiled = stencil(statement, statement % pattern, 4)

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case returned normally: ', trim(case_name)

end program refusal
