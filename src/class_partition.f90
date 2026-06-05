!=====================================================================!
! Domain partition + halo metadata for distributed-memory solves.
!
! Given the (already partitioned) dof/connectivity graph, organise the
! result into the bookkeeping a distributed solver needs:
!
!   part_of(cell)   -> owning part (image) 1..nparts
!   owned(k)        -> the cells part k owns        (its matrix rows)
!   ghosts(k)       -> the off-part cells part k needs (the columns its
!                      owned rows reference but does not own) - the halo
!
! ghosts(k) is exactly the set of cells in other parts adjacent (through a
! graph edge / interior face) to a cell owned by k. The owner of a ghost g
! is part_of(g), so a halo exchange pulls g from image part_of(g).
!
! This is pure serial integer bookkeeping (NO coarrays) computed identically
! on every image from the replicated graph, so all images agree on who owns
! and needs what without communicating. It is fully unit-testable in serial.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_partition

  use iso_fortran_env, only : dp => REAL64
  use class_graph    , only : graph

  implicit none

  private
  public :: partition

  type :: partition

     integer :: nparts = 1
     integer :: ncells = 0
     integer :: ncut   = 0                  ! edges crossing parts (cut quality)

     integer, allocatable :: part_of(:)     ! (ncells)   owning part of each cell

     ! owned cells per part, stored csr-style: part k owns
     ! own_list(own_ptr(k) : own_ptr(k+1)-1)
     integer, allocatable :: own_ptr(:)     ! (nparts+1)
     integer, allocatable :: own_list(:)    ! (ncells)

     ! ghost (halo) cells per part, same csr layout
     integer, allocatable :: gh_ptr(:)      ! (nparts+1)
     integer, allocatable :: gh_list(:)     ! (ncut-ish, deduped)

   contains

     procedure :: owned
     procedure :: ghosts
     procedure :: n_owned
     procedure :: n_ghosts
     procedure :: balance
     procedure :: print

  end type partition

  interface partition
     module procedure create_partition
  end interface partition

contains

  !===================================================================!
  ! Build the partition metadata from an already-partitioned graph
  ! (call g % partition(nparts) first - this reads g % vertices % part).
  !===================================================================!

  type(partition) function create_partition(g, nparts) result(p)

    type(graph), intent(in) :: g
    integer    , intent(in) :: nparts

    integer, allocatable :: ptr(:), cnt(:), mark(:)
    integer :: nc, ne, c, e, k, t, h, pt, ph, pos

    nc = g % num_vertices
    ne = g % num_edges
    p % nparts = nparts
    p % ncells = nc

    ! ---- part_of from the stamped vertices ----
    allocate(p % part_of(nc))
    do c = 1, nc
       p % part_of(c) = g % vertices(c) % part
    end do

    ! ---- owned cells per part (counting sort by part) ----
    allocate(p % own_ptr(nparts+1)); p % own_ptr = 0
    do c = 1, nc
       p % own_ptr(p % part_of(c)+1) = p % own_ptr(p % part_of(c)+1) + 1
    end do
    p % own_ptr(1) = 1
    do k = 1, nparts
       p % own_ptr(k+1) = p % own_ptr(k+1) + p % own_ptr(k)
    end do
    allocate(p % own_list(nc))
    allocate(ptr(nparts)); ptr = p % own_ptr(1:nparts)
    do c = 1, nc
       k = p % part_of(c)
       p % own_list(ptr(k)) = c
       ptr(k) = ptr(k) + 1
    end do
    deallocate(ptr)

    ! ---- edge cut ----
    p % ncut = 0
    do e = 1, ne
       if (p % part_of(g % edges(e) % tail) .ne. &
            & p % part_of(g % edges(e) % head)) p % ncut = p % ncut + 1
    end do

    ! ---- ghost cells per part ----
    ! ghost(k) = cells not in k but adjacent to a k-owned cell. Two passes:
    ! count (size the csr) then fill, deduped by stamping mark(cell)=k (k is
    ! monotone, so an old stamp from part k-1 reads as unmarked for part k).
    allocate(cnt(nparts)); cnt = 0
    allocate(mark(nc));    mark = 0
    do k = 1, nparts
       do e = 1, ne
          t = g % edges(e) % tail;  h = g % edges(e) % head
          pt = p % part_of(t);      ph = p % part_of(h)
          if (pt .eq. k .and. ph .ne. k) then
             if (mark(h) .ne. k) then; mark(h) = k; cnt(k) = cnt(k) + 1; end if
          end if
          if (ph .eq. k .and. pt .ne. k) then
             if (mark(t) .ne. k) then; mark(t) = k; cnt(k) = cnt(k) + 1; end if
          end if
       end do
    end do

    allocate(p % gh_ptr(nparts+1))
    p % gh_ptr(1) = 1
    do k = 1, nparts
       p % gh_ptr(k+1) = p % gh_ptr(k) + cnt(k)
    end do
    allocate(p % gh_list(p % gh_ptr(nparts+1)-1))

    mark = 0
    do k = 1, nparts
       pos = p % gh_ptr(k)
       do e = 1, ne
          t = g % edges(e) % tail;  h = g % edges(e) % head
          pt = p % part_of(t);      ph = p % part_of(h)
          if (pt .eq. k .and. ph .ne. k) then
             if (mark(h) .ne. k) then; mark(h) = k; p % gh_list(pos) = h; pos = pos + 1; end if
          end if
          if (ph .eq. k .and. pt .ne. k) then
             if (mark(t) .ne. k) then; mark(t) = k; p % gh_list(pos) = t; pos = pos + 1; end if
          end if
       end do
    end do

  end function create_partition

  !===================================================================!
  ! Cells owned by part k (its matrix rows)
  !===================================================================!

  pure function owned(this, k) result(cells)
    class(partition), intent(in) :: this
    integer         , intent(in) :: k
    integer, allocatable :: cells(:)
    cells = this % own_list(this % own_ptr(k) : this % own_ptr(k+1)-1)
  end function owned

  !===================================================================!
  ! Halo (ghost) cells part k needs from other parts
  !===================================================================!

  pure function ghosts(this, k) result(cells)
    class(partition), intent(in) :: this
    integer         , intent(in) :: k
    integer, allocatable :: cells(:)
    cells = this % gh_list(this % gh_ptr(k) : this % gh_ptr(k+1)-1)
  end function ghosts

  pure integer function n_owned(this, k)
    class(partition), intent(in) :: this
    integer         , intent(in) :: k
    n_owned = this % own_ptr(k+1) - this % own_ptr(k)
  end function n_owned

  pure integer function n_ghosts(this, k)
    class(partition), intent(in) :: this
    integer         , intent(in) :: k
    n_ghosts = this % gh_ptr(k+1) - this % gh_ptr(k)
  end function n_ghosts

  !===================================================================!
  ! Load-balance ratio max(owned)/min(owned) over the parts (1.0 = perfect)
  !===================================================================!

  pure real(dp) function balance(this)
    class(partition), intent(in) :: this
    integer :: k, lo, hi, m
    lo = huge(1); hi = 0
    do k = 1, this % nparts
       m = this % n_owned(k)
       lo = min(lo, m); hi = max(hi, m)
    end do
    if (lo .le. 0) then
       balance = huge(1.0_dp)
    else
       balance = real(hi, dp)/real(lo, dp)
    end if
  end function balance

  !===================================================================!
  ! One-line summary
  !===================================================================!

  subroutine print(this)
    class(partition), intent(in) :: this
    integer :: k
    write(*,'(1x,a,i0,a,i0,a,i0,a,f6.3)') "partition: ", this % nparts, &
         & " parts, ", this % ncells, " cells, edge cut ", this % ncut, &
         & ", balance ", this % balance()
    do k = 1, this % nparts
       write(*,'(3x,a,i0,a,i0,a,i0,a)') "part ", k, ": ", this % n_owned(k), &
            & " owned, ", this % n_ghosts(k), " ghost"
    end do
  end subroutine print

end module class_partition
