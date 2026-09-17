!=====================================================================!
! THE EXCHANGE: the distribution of the unknowns over the images and
! the operators it induces on vectors of the whole length. Every
! image has every vector at the whole length; the entries an image
! owns are its own, the entries it reads through the coupling but does
! not own are its halo, and every other entry is zero or stale and
! never read. The exchange is
!
!      update(x)     x(halo) := x(halo) of the owning image
!      gather(x)     x(i)    := x(i) of the owner of i, for every i
!      total(s)      the sum of s over the images
!      maximum(s)    the largest s over the images
!
! The ownership is given per unknown; the halo is derived from a
! coupling graph as the tails of the edges entering owned vertices
! that are not owned. Under one image nothing is exchanged.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_exchange

  use util_precision, only : dp
  use view_directed , only : directed_graph

  implicit none

  private
  public :: exchange

  type :: exchange

     integer :: image = 1
     integer :: images = 1
     integer, allocatable :: owner(:)
     logical, allocatable :: is_owned(:)
     integer, allocatable :: owned(:)
     ! the halo, grouped by the owning image: the entries owned by
     ! image p are halo(halo_first(p):halo_first(p + 1) - 1)
     integer, allocatable :: halo(:)
     integer, allocatable :: halo_first(:)

   contains

     procedure :: halo_from
     procedure :: update
     procedure :: gather
     procedure :: total
     procedure :: maximum
     procedure :: num_owned

  end type exchange

  interface exchange
     module procedure create
  end interface exchange

contains

  !===================================================================!
  ! The exchange of the ownership given, one owning image per
  ! unknown, with an empty halo until a coupling states one. An owner
  ! outside 1..num_images() is invalid input.
  !===================================================================!

  function create(owner) result(this)

    integer, intent(in) :: owner(:)
    type(exchange) :: this

    integer :: i
    character(len=250) :: message

    this % image  = this_image()
    this % images = num_images()
    if (any(owner < 1) .or. any(owner > this % images)) then
       write(message,'(a,i0,a,i0)') 'operation_exchange: every unknown is owned by one of the images &
            &1..', this % images, '; owners range to ', maxval(owner)
       error stop trim(message)
    end if
    this % owner    = owner
    this % is_owned = owner == this % image
    this % owned    = pack([(i, i = 1, size(owner))], this % is_owned)
    allocate(this % halo(0), this % halo_first(this % images + 1), source=0)
    this % halo_first = 1

  end function create

  pure integer function num_owned(this)
    class(exchange), intent(in) :: this
    num_owned = size(this % owned)
  end function num_owned

  !===================================================================!
  ! THE HALO from a coupling: the tails of the edges entering owned
  ! vertices that this image does not own, each once, grouped by
  ! owner. A coupling of another vertex count than the ownership is
  ! invalid input.
  !===================================================================!

  subroutine halo_from(this, coupling)

    class(exchange)      , intent(inout) :: this
    class(directed_graph), intent(in)    :: coupling

    logical, allocatable :: in_halo(:)
    integer, allocatable :: count(:)
    integer :: e, t, h, n, p, at
    character(len=250) :: message

    n = size(this % owner)
    if (coupling % num_vertices() /= n) then
       write(message,'(a,i0,a,i0)') 'operation_exchange: the coupling has one vertex per unknown; &
            &vertices = ', coupling % num_vertices(), ', unknowns = ', n
       error stop trim(message)
    end if

    allocate(in_halo(n), source=.false.)
    do e = 1, coupling % num_edges()
       if (.not. coupling % edge_has_head(e)) cycle
       h = coupling % edge_head(e)
       t = coupling % edge_tail(e)
       if (this % is_owned(h) .and. .not. this % is_owned(t)) in_halo(t) = .true.
    end do

    allocate(count(this % images), source=0)
    do t = 1, n
       if (in_halo(t)) count(this % owner(t)) = count(this % owner(t)) + 1
    end do
    this % halo_first(1) = 1
    do p = 1, this % images
       this % halo_first(p + 1) = this % halo_first(p) + count(p)
    end do
    if (allocated(this % halo)) deallocate(this % halo)
    allocate(this % halo(this % halo_first(this % images + 1) - 1))
    count = this % halo_first(1:this % images)
    do t = 1, n
       if (.not. in_halo(t)) cycle
       p  = this % owner(t)
       at = count(p)
       this % halo(at) = t
       count(p) = at + 1
    end do

  end subroutine halo_from

  !===================================================================!
  ! THE HALO UPDATE AND THE GATHER, both by one sum over the images
  ! of the vector with every entry not owned set to zero: each entry
  ! is owned by one image, so the sum is the owner's value everywhere.
  ! A reduction of the whole length costs more than the halo alone
  ! would, and needs no coarray variable: the coarray runtime in use
  ! (OpenCoarrays 2.9.2 under gfortran 14) registers a module coarray
  ! before the runtime is initialised and mis-sizes a vector-subscript
  ! access through a component, while its collectives are sound.
  !===================================================================!

  subroutine update(this, x)

    class(exchange), intent(in)    :: this
    real(dp)       , intent(inout) :: x(:)

    call this % gather(x)

  end subroutine update

  subroutine gather(this, x)

    class(exchange), intent(in)    :: this
    real(dp)       , intent(inout) :: x(:)

    if (this % images == 1) return
    where (.not. this % is_owned) x = 0.0_dp
    call co_sum(x)

  end subroutine gather

  real(dp) function total(this, s)
    class(exchange), intent(in) :: this
    real(dp)       , intent(in) :: s
    total = s
    if (this % images > 1) call co_sum(total)
  end function total

  real(dp) function maximum(this, s)
    class(exchange), intent(in) :: this
    real(dp)       , intent(in) :: s
    maximum = s
    if (this % images > 1) call co_max(maximum)
  end function maximum

end module operation_exchange
