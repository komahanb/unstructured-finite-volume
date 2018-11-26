!=====================================================================!
! Basic implementation of an unordered list tuple
!=====================================================================!

module class_list

  implicit none

  ! Datatype
  type :: list
     integer, allocatable :: table(:,:)
     integer              :: num_entries
     integer              :: num_tuples
   contains
     procedure :: insert
     procedure :: add_entry
     procedure :: get_entries

     ! Destructor
     final :: destroy

  end type list

  ! Constructor interface for list
  interface list
     module procedure create
  end interface list
  
contains

  !===================================================================!
  ! Constructor implementaion of list
  !===================================================================!

  type(list) function create(num_tuples, max_entries) result(this)

    integer, intent(in) :: max_entries
    integer, intent(in) :: num_tuples

    this % num_tuples = num_tuples
    
    allocate(this % table(num_tuples, max_entries))
    this % table = 0
    this % num_entries = 0
    
  end function create

  !===================================================================!
  ! Destructor for list object
  !===================================================================!
  
  pure subroutine destroy(this)

    type(list), intent(inout) :: this

    if(allocated(this % table)) deallocate(this % table)

  end subroutine destroy

  !===================================================================!
  ! Add an entry into the list
  !===================================================================!

  type(logical) function insert(this, tuple)

    class(list), intent(inout) :: this
    integer   , intent(in)    :: tuple(:)

    ! Check if tuple is in table   
    this % num_entries = this % num_entries + 1
    this % table(:, this % num_entries) = tuple(:)
    insert = .true.

  end function insert

  !===================================================================!
  ! Add an entry into the list
  !===================================================================!

  subroutine add_entry(this, tuple)

    class(list), intent(inout) :: this
    integer   , intent(in)    :: tuple(:)

    ! Check if tuple is in table    
    this % num_entries = this % num_entries + 1
    this % table(:, this % num_entries) = tuple(:)

  end subroutine add_entry
  
  !===================================================================!
  ! Get all the entries in the list as an array
  !===================================================================!
  
  pure subroutine get_entries(this, entries)

    class(list), intent(in)            :: this
    integer, allocatable, intent(out) :: entries(:,:)

    allocate(entries(this % num_tuples, this % num_entries))

    entries(:,:) = this % table(:,1:this % num_entries)

  end subroutine get_entries
  
end module class_list

subroutine test_list

  use class_list, only : list

  implicit none

  type(list) :: faces
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

  faces = list(2,16)    

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

end subroutine test_list
