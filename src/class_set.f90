!=====================================================================!
! Basic implementation of an unordered set tuple
!=====================================================================!

module class_set

  implicit none

  ! Datatype
  type :: set
     integer, allocatable :: table(:,:)
     integer              :: num_entries
     integer              :: num_tuples
   contains
     procedure :: insert
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

  type(set) function create(num_tuples, max_entries) result(this)

    integer, intent(in) :: max_entries
    integer, intent(in) :: num_tuples

    this % num_tuples = num_tuples
    
    allocate(this % table(num_tuples, max_entries))
    this % table = 0
    this % num_entries = 0
    
  end function create

  !===================================================================!
  ! Add an entry into the set
  !===================================================================!

  type(logical) function insert(this, tuple)

    class(set), intent(inout) :: this
    integer   , intent(in)    :: tuple(:)

    ! Check if tuple is in table   
    if (this % contains(tuple) .eqv. .false.) then
       this % num_entries = this % num_entries + 1
       this % table(:, this % num_entries) = tuple(:)
       insert = .true.
    else
       insert = .false.
    end if

  end function insert

  !===================================================================!
  ! Add an entry into the set
  !===================================================================!

  subroutine add_entry(this, tuple)

    class(set), intent(inout) :: this
    integer   , intent(in)    :: tuple(:)

    ! Check if tuple is in table

    if (this % contains(tuple) .eqv. .false.) then
       this % num_entries = this % num_entries + 1
       this % table(:, this % num_entries) = tuple(:)
    end if

  end subroutine add_entry

   !===================================================================!
  ! Checks if the first argument is a subset of the second argument
  ! (move elsewhere?)
  !===================================================================!
  
  pure type(logical) function is_subset(small, big)

    integer, intent(in) :: small(:)
    integer, intent(in) :: big(:)

    integer, allocatable :: sub(:)
    integer, allocatable :: set(:)

    integer :: lensub, i

    is_subset = .false.

    if (size(small) .gt. size(big)) return

    ! Create local copy of arrays
    allocate(sub, source = small)
    allocate(set, source = big)

    ! Sort two arrays
    call isort(sub)
    call isort(set)    
    lensub = size(sub)

    ! Check if all entries are equal upto the length of the smallest
    ! array
    is_subset = .true.
    check : do i = 1, lensub
       if (any(set .eq. sub(i)) .eqv. .false.) then
          is_subset = .false.
          exit check
       end if
    end do check

    deallocate(sub,set)    

  end function is_subset
  
  !===================================================================!
  ! Sort an integer array ! move elsewhere?
  !===================================================================!
  
  pure subroutine isort(array)

    integer, intent(inout) :: array(:)
    integer :: temp , j , k
    integer :: n

    n = ubound(array,1)

    do j = 1 , n
       do k = j + 1 , n
          if(array(j) > array(k)) then
             temp     = array(k)
             array(k) = array(j)
             array(j) = temp
          end if
       end do
    end do

  end subroutine isort

  !===================================================================!
  ! Check if an entry is contains in the set
  !===================================================================!

  pure type(logical) function contains(this, tuple)

    class(set), intent(in) :: this
    integer   , intent(in) :: tuple(:)
    integer :: i

    contains = .false.

    ! loop through existing tuples and find if it exists
   check: do i = this % num_entries, 1, -1
       ! Improve logic. This is expensive and unnnecessary
       if (is_subset(tuple, this % table(:,i)) .eqv. .true.) then
          contains = .true.
          return
       end if
    end do check

  end function contains

  !===================================================================!
  ! Get all the entries in the set as an array
  !===================================================================!
  
  pure subroutine get_entries(this, entries)
    
    class(set), intent(in)            :: this
    integer, allocatable, intent(out) :: entries(:,:)
   
    allocate(entries(this % num_tuples, this % num_entries))
    
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

  faces = set(2,16)    

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
