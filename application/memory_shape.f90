! What the representation costs, part by part.
!
! One part per run, because within a single process the allocator's
! arena is already warm and a difference in resident memory reports
! nothing. Each run builds exactly one thing and reports the peak the
! kernel recorded, so the parts are told apart by comparing runs
! against the run that builds nothing.
!
!   memory_shape <none|vertices|edges|field|block> <instants>
program memory_shape

  use util_precision  , only : dp
  use view_directed_stored  , only : stored_directed_graph
  use field_stored          , only : stored_field
  use operation_family_bdf  , only : bdf_family
  use physics_vanderpol     , only : van_der_pol
  use gti_march             , only : block_from, partition
  use gti_block             , only : block_residual
  use gti_expansion         , only : expansion, family_holder
  use operation_grid        , only : uniform_grid

  implicit none

  integer, parameter :: degrees = 3, order = 2

  type(stored_directed_graph) :: gr
  type(stored_field)          :: over
  type(block_residual)        :: rows
type(expansion) :: tower
type(family_holder) :: holder(1)
integer, allocatable :: at(:)
  type(bdf_family)            :: scheme

  character(len=32) :: what, given
  integer , allocatable :: tails(:), heads(:)
  real(dp), allocatable :: dt(:), t(:), held(:)
  integer :: instants, n, m, h, band, i, j, e

  call get_command_argument(1, what)
  call get_command_argument(2, given)
  read(given,*) instants

  scheme = bdf_family(order)
  h      = scheme % history_depth(degrees - 1)
  n      = (instants - h) * degrees
  band   = h * degrees

  e = 0
  do i = 1, n
     do j = max(1, i - band), i
        e = e + 1
     end do
  end do
  m = e

  allocate(tails(m), heads(m))
  e = 0
  do i = 1, n
     do j = max(1, i - band), i
        e = e + 1
        tails(e) = j
        heads(e) = i
     end do
  end do

  select case (trim(what))
  case ('none')
     continue
  case ('vertices')
     gr = stored_directed_graph(n, tails=[integer ::], heads=[integer ::])
  case ('edges')
     gr = stored_directed_graph(n, tails=tails, heads=heads)
  case ('field')
     gr   = stored_directed_graph(n, tails=tails, heads=heads)
     over = stored_field('x', gr % vertex_set(), n)
     call over % set_real_vector(spread(1.0_dp, 1, n))
  case ('block')
     call partition(3.0_dp, instants, dt, t)
     allocate(held(h * degrees), source=0.0_dp)
     allocate(holder(1) % scheme, source=scheme)
     call tower % build(van_der_pol(degrees - 1), holder, [instants], uniform_grid(3.0_dp), 0, 0.0_dp)
     call block_from(tower, 1, scheme, van_der_pol(degrees - 1), held, rows, at)
  case default
     error stop 'memory_shape: the part is none, vertices, edges, field or block'
  end select

  write(*,'(a,i8,i9,i9,f12.3)') trim(what), instants, n, m, peak()

contains

  !-------------------------------------------------------------------!
  ! The peak resident size the kernel recorded, in megabytes.
  !-------------------------------------------------------------------!

  real(dp) function peak() result(mb)

    integer :: u, status, kb
    character(len=80) :: line

    mb = 0.0_dp
    open(newunit=u, file='/proc/self/status', action='read')
    do
       read(u,'(a)',iostat=status) line
       if (status /= 0) exit
       if (line(1:6) == 'VmHWM:') then
          read(line(7:),*) kb
          mb = real(kb, dp) / 1000.0_dp
          exit
       end if
    end do
    close(u)

  end function peak

end program memory_shape
