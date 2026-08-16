!=====================================================================!
! PARTITIONED IMPLICIT PDE TOWER . LEVEL 8 . THE REFUSAL
!
! The composite is DECOMPOSITION-CONTEXT-BOUND: its G1 and G2 were
! cut from one particular chain. Gate B's six-vertex star has the
! same cardinality and a different shape, and it is refused - in
! domain() as well as apply(), so a solver attaching on it dies at
! the earliest honest contract point rather than deep inside a
! matvec with a chain-derived decomposition.
!
!   foreign-host-apply ..... A_part.apply(G_alt, ...)
!   foreign-host-attach .... gmres.attach(A_part, G_alt, V(G_alt))
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program partitioned_pde_level_8_refusal

  use iso_fortran_env  , only : dp => REAL64
  use partitioned_pde_assert, only : NV, Q_EXACT
  use graph_field_calculus, only : graph_field
  use class_graph      , only : stored_graph
  use class_graph_field, only : field
  use class_graph_gmres, only : gmres
  use partitioned_shifted_laplacian_fixture, only : &
       & partitioned_shifted_laplacian

  implicit none

  type(stored_graph)                  :: g, g_alt
  type(partitioned_shifted_laplacian) :: composite
  type(gmres)                         :: solver
  type(field)                         :: q
  class(graph_field), allocatable     :: out
  character(len=32)                   :: which

  if (command_argument_count() .lt. 1) then
     error stop 'usage: refusal <case>'
  end if
  call get_command_argument(1, which)

  g     = stored_graph(NV, tails=[1,2,3,4,5], heads=[2,3,4,5,6])
  g_alt = stored_graph(NV, tails=[1,1,1,1,1], heads=[2,3,4,5,6])

  composite = partitioned_shifted_laplacian(g)

  select case (trim(which))

  case ('foreign-host-apply')
     q = field('state on the star', g_alt % vertex_set(), g_alt % num_vertices())
     call q % set_real_vector(Q_EXACT)
     call composite % apply(g_alt, [q], out)
     write(*,*) 'a chain decomposition acted on the star'

  case ('foreign-host-attach')
     ! attach asks the action for its domain; the guard fires there.
     call solver % attach(composite, g_alt, g_alt % vertex_set(), g_alt % num_vertices())
     write(*,*) 'a chain decomposition was attached to a star host'

  case default
     error stop 'usage: refusal <case>'

  end select

end program partitioned_pde_level_8_refusal
