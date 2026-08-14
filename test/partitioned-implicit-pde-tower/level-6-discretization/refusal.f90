!=====================================================================!
! PARTITIONED IMPLICIT PDE TOWER . LEVEL 6 . THE REFUSAL
!
! One case, and it is the same law every tower in this repository
! has had to learn at its own radius:
!
!   foreign-carrier .... a state field of exactly the RIGHT SIZE,
!                        living on a carrier that is not this
!                        graph's vertex set, must be refused
!
! The foreign carrier holds six members, as V(G) does. Nothing but
! IDENTITY separates them, and a size check would wave it through
! to compute a confident wrong answer over the wrong topology.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program partitioned_pde_level_6_refusal

  use iso_fortran_env  , only : dp => REAL64
  use partitioned_pde_assert, only : NV, Q_EXACT
  use graph_carrier    , only : counted_set
  use graph_grammar    , only : graph_field
  use class_graph      , only : stored_graph
  use class_graph_field, only : field
  use shifted_laplacian_fixture, only : shifted_laplacian

  implicit none

  type(stored_graph)              :: g
  type(counted_set)               :: foreign
  type(shifted_laplacian)         :: shifted
  type(field)                     :: q
  class(graph_field), allocatable :: out
  character(len=32)               :: which

  if (command_argument_count() .lt. 1) then
     error stop 'usage: refusal <case>'
  end if
  call get_command_argument(1, which)

  g = stored_graph(NV, tails=[1,2,3,4,5], heads=[2,3,4,5,6])

  select case (trim(which))

  case ('foreign-carrier')
     ! Six members, like V(G) - and not V(G).
     foreign = counted_set('a carrier of the right size', NV)
     q = field('state on a foreign carrier', foreign)
     call q % set_real_vector(Q_EXACT)
     call shifted % apply(g, [q], out)
     write(*,*) 'a foreign carrier of the right size was accepted'

  case default
     error stop 'usage: refusal <case>'

  end select

end program partitioned_pde_level_6_refusal
