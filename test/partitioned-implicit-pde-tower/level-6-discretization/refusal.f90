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
  use graph_fractal        , only : graph
  use map_set_representation, only : counted_set_representation
  use map_set        , only : set_map
  use field_calculus, only : field
  use view_directed_stored      , only : stored_directed_graph
  use field_stored, only : stored_field
  use shifted_laplacian_fixture, only : shifted_laplacian

  implicit none

  type(stored_directed_graph)              :: g
  type(graph)                 :: foreign
  type(set_map)                   :: sets
  type(shifted_laplacian)         :: shifted
  type(stored_field)                     :: q
  class(field), allocatable :: out
  character(len=32)               :: which

  if (command_argument_count() .lt. 1) then
     error stop 'usage: refusal <case>'
  end if
  call get_command_argument(1, which)

  g = stored_directed_graph(NV, tails=[1,2,3,4,5], heads=[2,3,4,5,6])
  call sets % bind(g % vertex_set(), &
       & counted_set_representation(g % num_vertices()))
  call sets % bind(g % edge_set(), &
       & counted_set_representation(g % num_edges()))

  select case (trim(which))

  case ('foreign-carrier')
     ! Six members, like V(G) - and not V(G). The whole point is
     ! that identity, not cardinality, decides.
     call foreign % declare()
     call sets % bind(foreign, counted_set_representation(NV))
     q = stored_field('state on a foreign carrier', foreign, NV)
     call q % set_real_vector(Q_EXACT)
     call shifted % apply(g, [q], out)
     write(*,*) 'a foreign carrier of the right size was accepted'

  case default
     error stop 'usage: refusal <case>'

  end select

end program partitioned_pde_level_6_refusal
