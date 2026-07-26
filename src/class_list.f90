!=====================================================================!
! A basic implementation of an unordered list of tuples.
!=====================================================================!

module class_list

  implicit none

  ! The datatype.
  type :: list
     integer, allocatable :: table(:,:)
     integer              :: num_entries
     integer              :: num_tuples
   contains
     procedure :: insert
     procedure :: add_entry
     procedure :: get_entries

     ! The destructor.
     final :: destroy

  end type list

  ! The constructor interface for a list.
  interface list
     module procedure create
  end interface list
  
contains

  !===================================================================!
  ! The constructor implementation of a list.
  !===================================================================!

  pure type(list) function create(num_tuples, max_entries) result(this)

    integer, intent(in) :: max_entries
    integer, intent(in) :: num_tuples

    this % num_tuples = num_tuples
    
    allocate(this % table(num_tuples, max_entries))
    this % table = 0
    this % num_entries = 0
    
  end function create

  !===================================================================!
  ! The destructor for a list object.
  !===================================================================!
  
  pure subroutine destroy(this)

    type(list), intent(inout) :: this

    if(allocated(this % table)) deallocate(this % table)

  end subroutine destroy

  !===================================================================!
  ! Append a tuple to the list and report success.
  !===================================================================!

  impure type(logical) function insert(this, tuple)

    class(list), intent(inout) :: this
    integer   , intent(in)    :: tuple(:)

    ! Append the tuple to the table without a membership check.
    this % num_entries = this % num_entries + 1
    this % table(:, this % num_entries) = tuple(:)
    insert = .true.

  end function insert

  !===================================================================!
  ! Add an entry into the list.
  !===================================================================!

  pure subroutine add_entry(this, tuple)

    class(list), intent(inout) :: this
    integer   , intent(in)    :: tuple(:)

    ! Append the tuple to the table without a membership check.
    this % num_entries = this % num_entries + 1
    this % table(:, this % num_entries) = tuple(:)

  end subroutine add_entry
  
  !===================================================================!
  ! Get all the entries in the list as an array.
  !===================================================================!
  
  pure subroutine get_entries(this, entries)

    class(list), intent(in)            :: this
    integer, allocatable, intent(out) :: entries(:,:)

    allocate(entries(this % num_tuples, this % num_entries))

    entries(:,:) = this % table(:,1:this % num_entries)

  end subroutine get_entries
  
end module class_list
