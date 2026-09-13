!=====================================================================!
! Counted ownership of immutable storage.
!
! A storage cell S records its version and the set B(S) of live
! bindings, one serial per owner. A counted reference r is bound to S
! by acquisition or by defined assignment, which adds a serial to
! B(S); it is released by finalization or by being assigned over,
! which removes its serial. The contents of S are cleared when B(S)
! becomes empty. The cell object itself is retained and recycled: a
! reference whose binding is stale addresses a cell that still
! exists, the versions differ, and the reference reports that it is
! not live instead of reading released storage.
!
! Release is idempotent in the serial. gfortran 15 assigns a
! containing object by finalizing the destination component in
! place, then calling the component's defined assignment on a
! temporary that still contains the destination's former bytes: the
! finalizer removes the serial, and the later release of the
! temporary finds nothing to remove.
!
! Owners are the references bound through acquire or through defined
! assignment: intrinsic assignment of the reference or of any object
! containing it, function results, arrays, move_alloc. Fortran also
! copies a value bitwise without defined assignment under
! allocate(source=variable), structure constructors, polymorphic
! intrinsic assignment and array constructors of containing objects;
! gfortran 15 finalizes some of those copies. Such a copy is not an
! owner: it has its twin's serial, its finalization releases that
! binding, and the twin is refused at its next access. The measured
! table of mechanisms is in doc/topology-ownership.md.
!
! The version and binding serial counters, the free list of recycled
! cells and every cell's binding list form one registry. Its writers
! (acquire, assign, release) run inside one named OpenMP critical
! region when the library is built with -fopenmp; without it the
! sentinel lines are comments and the serial statements are unchanged.
! An owner's deferred clear runs outside the region, because a clear
! may finalize other counted references and the region is not
! re-entrant; the cleared cell joins the free list in a second entry.
! Reads by an owner (live, num_owners, storage) take no lock: no other
! thread can clear a cell on which this thread has a live binding.
!=====================================================================!

module util_counted_storage

  implicit none

  private

  public :: counted_storage, counted_reference

  !===================================================================!
  ! The cell: a version, the live binding serials, and deferred
  ! clearing of the owner's own contents.
  !===================================================================!

  type, abstract :: counted_storage
     integer, private :: version = 0
     integer, private :: num_owners = 0
     integer, allocatable, private :: bindings(:)
   contains
     procedure(clear_interface), deferred :: clear
  end type counted_storage

  abstract interface
     subroutine clear_interface(this)
       import counted_storage
       class(counted_storage), intent(inout) :: this
     end subroutine clear_interface
  end interface

  type :: cell_link
     class(counted_storage), pointer :: cell => null()
     type(cell_link), pointer :: next => null()
  end type cell_link

  !===================================================================!
  ! The reference: a cell, the version it was bound to, and the
  ! serial of its own binding.
  !===================================================================!

  type :: counted_reference
     class(counted_storage), pointer, private :: cell => null()
     integer, private :: version = 0
     integer, private :: binding = 0
   contains
     procedure :: acquire
     procedure :: live
     procedure :: released
     procedure :: num_owners
     procedure :: storage
     procedure, private :: assign
     generic :: assignment(=) => assign
     final :: release
  end type counted_reference

  type(cell_link), pointer :: released_cells => null()
  integer :: last_version = 0
  integer :: last_binding = 0

contains

  !===================================================================!
  ! Bind this reference to a new cell of the template's dynamic
  ! type, as its one owner. A released cell of that type is reused;
  ! otherwise one is allocated. Any previous binding is released.
  !===================================================================!

  subroutine acquire(this, template)

    class(counted_reference), intent(inout) :: this
    class(counted_storage) , intent(in)    :: template

    type(cell_link), pointer :: link, previous
    class(counted_storage), pointer :: cleared

    !$omp critical(counted_storage_registry)
    call remove_binding(this, cleared)
    !$omp end critical(counted_storage_registry)
    call recycle(cleared)

    !$omp critical(counted_storage_registry)
    previous => null()
    link => released_cells
    do while (associated(link))
       if (same_type_as(link % cell, template)) exit
       previous => link
       link => link % next
    end do

    if (associated(link)) then
       this % cell => link % cell
       if (associated(previous)) then
          previous % next => link % next
       else
          released_cells => link % next
       end if
       deallocate(link)
    else
       allocate(this % cell, mold=template)
    end if

    last_version = last_version + 1
    this % version = last_version
    this % cell % version = last_version
    this % cell % num_owners = 0
    if (.not. allocated(this % cell % bindings)) allocate(this % cell % bindings(4))
    this % binding = new_binding(this % cell)
    !$omp end critical(counted_storage_registry)

  end subroutine acquire

  ! Register one more owner on a live cell and return its serial.
  ! Called inside the registry's critical region.
  integer function new_binding(cell) result(binding)

    class(counted_storage), intent(inout) :: cell

    integer, allocatable :: larger(:)
    character(len=250) :: message

    if (cell % num_owners == size(cell % bindings)) then
       if (cell % num_owners > huge(cell % num_owners) / 2) then
          write(message,'(a,i0)') 'util_counted_storage: doubling the bindings array would &
               &exceed the integer range; num_owners = ', cell % num_owners
          error stop trim(message)
       end if
       allocate(larger(2 * cell % num_owners))
       larger(1:cell % num_owners) = cell % bindings(1:cell % num_owners)
       call move_alloc(larger, cell % bindings)
    end if
    last_binding = last_binding + 1
    cell % num_owners = cell % num_owners + 1
    cell % bindings(cell % num_owners) = last_binding
    binding = last_binding

  end function new_binding

  pure logical function live(this)

    class(counted_reference), intent(in) :: this

    live = associated(this % cell)
    if (live) live = this % cell % version == this % version .and. this % version > 0

  end function live

  ! Bound once and no longer live: the cell was cleared by its last
  ! owner, or this value is a bitwise copy whose twin released the
  ! binding. A default reference has never been bound and is not
  ! released.
  pure logical function released(this)

    class(counted_reference), intent(in) :: this

    released = this % version > 0 .and. .not. live(this)

  end function released

  ! Zero for a reference that is not live.
  pure integer function num_owners(this)

    class(counted_reference), intent(in) :: this

    num_owners = 0
    if (live(this)) num_owners = this % cell % num_owners

  end function num_owners

  ! The owner fills and reads the cell through this pointer; the
  ! reference does not interpret its contents.
  function storage(this) result(cell)

    class(counted_reference), intent(in) :: this
    class(counted_storage), pointer :: cell

    cell => null()
    if (live(this)) cell => this % cell

  end function storage

  !===================================================================!
  ! Defined assignment: the destination joins the owners of the
  ! source's cell with a serial of its own, after leaving its own
  ! cell. Elemental, so that arrays and containing objects of any
  ! rank assign through it. The source is read before the
  ! destination is released, so assigning a value to itself or to
  ! another owner of the same cell keeps the cell live.
  !===================================================================!

  impure elemental subroutine assign(lhs, rhs)

    class(counted_reference), intent(inout) :: lhs
    type(counted_reference) , intent(in)    :: rhs

    class(counted_storage), pointer :: cleared

    call bind_over(lhs, rhs, cleared)
    call recycle(cleared)

  end subroutine assign

  ! The registry part of assignment: one more binding on the source's
  ! cell, then the destination's own binding removed.
  subroutine bind_over(lhs, rhs, cleared)

    class(counted_reference), intent(inout) :: lhs
    type(counted_reference) , intent(in)    :: rhs
    class(counted_storage), pointer, intent(out) :: cleared

    class(counted_storage), pointer :: cell
    integer :: version, binding

    !$omp critical(counted_storage_registry)
    cell => null()
    version = 0
    binding = 0
    if (live(rhs)) then
       cell => rhs % cell
       version = rhs % version
       binding = new_binding(cell)
    end if
    call remove_binding(lhs, cleared)
    lhs % cell => cell
    lhs % version = version
    lhs % binding = binding
    !$omp end critical(counted_storage_registry)

  end subroutine bind_over

  !===================================================================!
  ! Finalization: remove this binding's serial; the last owner clears
  ! the cell and retains it for reuse. A reference whose serial is
  ! not registered, because it was already released or because it is
  ! a bitwise copy of a released owner, changes nothing.
  !===================================================================!

  impure elemental subroutine release(this)

    type(counted_reference), intent(inout) :: this

    class(counted_storage), pointer :: cleared

    call release_binding(this, cleared)
    call recycle(cleared)

  end subroutine release

  subroutine release_binding(this, cleared)

    type(counted_reference), intent(inout) :: this
    class(counted_storage), pointer, intent(out) :: cleared

    !$omp critical(counted_storage_registry)
    call remove_binding(this, cleared)
    !$omp end critical(counted_storage_registry)

  end subroutine release_binding

  ! Remove this reference's serial from its cell and unbind the
  ! reference. When the serial was the last one the cell's version
  ! becomes zero, so no reference is live on it, and the cell is
  ! returned for clearing; otherwise null. Called inside the
  ! registry's critical region.
  subroutine remove_binding(this, cleared)

    type(counted_reference), intent(inout) :: this
    class(counted_storage), pointer, intent(out) :: cleared

    integer :: k, n

    cleared => null()
    if (live(this)) then
       n = this % cell % num_owners
       do k = 1, n
          if (this % cell % bindings(k) == this % binding) exit
       end do
       if (k <= n) then
          this % cell % bindings(k) = this % cell % bindings(n)
          this % cell % num_owners = n - 1
          if (n == 1) then
             this % cell % version = 0
             cleared => this % cell
          end if
       end if
    end if
    this % cell => null()
    this % version = 0
    this % binding = 0

  end subroutine remove_binding

  ! Clear a cell that has no owner and place it on the free list. The
  ! owner's clear runs outside the registry's critical region.
  subroutine recycle(cell)

    class(counted_storage), pointer, intent(in) :: cell

    type(cell_link), pointer :: link

    if (.not. associated(cell)) return
    call cell % clear()
    allocate(link)
    link % cell => cell
    !$omp critical(counted_storage_registry)
    link % next => released_cells
    released_cells => link
    !$omp end critical(counted_storage_registry)

  end subroutine recycle

end module util_counted_storage
