!=====================================================================!
! The computational-state refusals, each EXPECTED TO DIE for its
! stated reason:
!
!      databottom       asking data of a bottom data seat
!      residualbottom   asking residual of a bottom residual seat
!      unhosted         asking structure of a graph riding nothing
!      twice            a second signing of a declared graph
!      fifth            asking the name of a state that is not one
!                       of the four
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program state_refusal

  use graph_structure, only : relational_graph
  use graph_state    , only : computational_graph, state_name, &
       &                      void_graph, data_graph

  implicit none

  type :: observed_values

     real, allocatable :: entries(:)

  end type observed_values

  type(observed_values)              :: q
  type(computational_graph), target  :: g
  class(*), pointer                  :: seat
  class(relational_graph), pointer   :: host
  character(len=:), allocatable      :: said
  character(len=32)                  :: which

  which = ''
  call get_command_argument(1, which)

  select case (trim(which))

  case ('databottom')
     g = void_graph('empty-handed')
     seat => g % data()

  case ('residualbottom')
     allocate(q % entries, source=[1.0])
     g = data_graph('data only', q)
     seat => g % residual()

  case ('unhosted')
     g = void_graph('adrift')
     host => g % structure()
     write(*,'(1x,a)') host % name()

  case ('twice')
     g = void_graph('signed')
     call g % declare('signed-again')

  case ('fifth')
     said = state_name(5)

  case default
     write(*,'(1x,a)') &
          & "usage: refusal databottom|residualbottom|unhosted|twice|fifth"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program state_refusal
