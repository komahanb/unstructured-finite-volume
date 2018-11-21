!=====================================================================!
!============== Basic implementation of a set tuple ==================!
!=====================================================================!

module class_set

  implicit none

  ! Datatype
  type :: set
     integer, allocatable :: table(:,:)
     integer              :: num_entries
   contains     
     procedure :: add_entry
     procedure :: get_entries
     procedure :: contains
  end type set

  ! Constructor interface for set
  interface set
     module procedure create
  end interface set
  
contains

  !===================================================================!
  ! Constructor implementaion of set
  !===================================================================!

  type(set) function create(max_entries) result(this)

    integer, intent(in) :: max_entries

    allocate(this % table(2,max_entries))
    this % table = 0
    this % num_entries = 0

  end function create

  !===================================================================!
  ! Add an entry into the set
  !===================================================================!

  subroutine add_entry(this, edge)

    class(set), intent(inout) :: this
    integer, intent(in) :: edge(2)

    ! Check if edge is in table

    if (this % contains(edge) .eqv. .false.) then
       this % num_entries = this % num_entries + 1
       this % table(:, this % num_entries) = edge(:)
    end if

  end subroutine add_entry

  !===================================================================!
  ! Check if an entry is contains in the set
  !===================================================================!

  pure type(logical) function contains(this, edge)

    class(set), intent(in) :: this
    integer   , intent(in) :: edge(2)
    integer :: i
    
    ! loop through enties and see if they exist
    contains = .false.

    do i = 1, this % num_entries
       if (edge(1) .eq. this % table (1,i) .and. edge(2) .eq. this % table(2,i)) then
          contains = .true.
       else if (edge(2) .eq. this % table (1,i) .and. edge(1) .eq. this % table(2,i)) then
          contains = .true.                    
       end if
    end do

  end function contains

  !===================================================================!
  ! Get all the entries in the set as an array
  !===================================================================!
  
  pure subroutine get_entries(this, entries)
    
    class(set), intent(in)            :: this
    integer, allocatable, intent(out) :: entries(:,:)

    allocate(entries(2, this % num_entries))
    
    entries(:,:) = this % table(:,1:this % num_entries)

  end subroutine get_entries
  
end module class_set

subroutine test_set

  use class_set, only : set

  implicit none

  type(set) :: faces
  integer   :: i
  integer   :: conn(4,4)
  integer   :: idx1, idx2, ipair, ielem
  integer, parameter :: nelems = 4
  integer   :: num_elem_faces(nelems)
  integer, allocatable :: table(:,:)

  ! Each element has four edges
  num_elem_faces = 4

!!$  7---8---9
!!$  |   |   |
!!$  4---5---6
!!$  |   |   |
!!$  1---2---3

  faces = set(16)    

  conn(:,1) = [1,2,5,4]
  conn(:,2) = [2,3,6,5]
  conn(:,3) = [4,5,8,7]
  conn(:,4) = [5,6,9,8]   

  ! Make ordered pair of nodes
  do ielem = 1, nelems
     do ipair = 1, num_elem_faces(ielem)

        if (ipair .eq. num_elem_faces(ielem)) then           
           idx1 = ipair
           idx2 = 1           
        else           
           idx1 = ipair
           idx2 = ipair+1           
        end if

        ! Add ordered pair of integers to an array
        call faces % add_entry([conn(idx1, ielem), conn(idx2, ielem)])

     end do
  end do

  call faces % get_entries(table)
  do i = 1, 12
     print *, table(:,i)
  end do
  print *, faces % num_entries

end subroutine test_set
