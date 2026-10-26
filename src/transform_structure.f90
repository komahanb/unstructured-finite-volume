!=====================================================================!
! Module that states the transform prime: the operation between
! graphs, graph -> graph, with data mapped along. At this level it has
! one symbol, an admissibility predicate: may this transform act on
! that graph. The concrete transforms, partition, assemble, coarsen
! and refine, are each checked by a round-trip law:
!
!      exact        assemble(partition(G)) = G      both ways
!      one-sided    coarsen(refine(G)) = G          one way only
!
!          partition                       refine
!       G ----------> parts             G ------->  G'
!       ^               |               ^           |
!       |   assemble    |               |  coarsen  |
!       +---------------+               +-----------+
!
! The interfaces are impure for a language constraint that no longer
! binds: F2018 C1594, on copying inside a pure subprogram a set that
! contains pointers. An impure interface permits a pure override, and
! every implementation is pure.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module transform_structure

  use util_precision, only : dp
  use view_directed , only : directed_graph

  implicit none

  private

  public :: transform
  public :: through_blocks

  !-------------------------------------------------------------------!
  ! Type that every transform between graphs extends
  !-------------------------------------------------------------------!

  type, abstract :: transform

   contains

     ! Deferred procedures
     procedure(transform_on_graph_interface), deferred :: defined_on_graph

  end type transform

  abstract interface

     !----------------------------------------------------------------!
     ! Interface for the admissibility predicate: whether the transform
     ! may act on the graph
     !----------------------------------------------------------------!

     logical function transform_on_graph_interface(this, input_graph)
       import :: transform, directed_graph
       class(transform)     , intent(in) :: this
       class(directed_graph), intent(in) :: input_graph
     end function transform_on_graph_interface

  end interface

contains

  !===================================================================!
  ! The one map between a level and its blocks, read both ways.
  !
  ! Read forward, it is the restriction R: each block the sum of its
  ! members, or their mean where average is requested. Read backward,
  ! it is R^T: every member takes its block's value, the injected
  ! prolongation. Values are num_components wide per member:
  !
  !   block_of = [1, 1, 2, 2, 2]
  !
  !   fine     [ a ] [ b ] [ c ] [ d ] [ e ]
  !              \     /     \     |     /
  !   R           \   /       \    |    /
  !   coarse      [a + b]     [c + d + e]      (/ 2 and / 3 with average)
  !
  !   coarse      [ A ]       [ B ]
  !              +--+--+     +--+--+-----+
  !   R^T        |     |     |     |     |
  !   fine     [ A ] [ A ] [ B ] [ B ] [ B ]
  !
  ! coarsen(refine(G)) = G is this map read twice, which is why it is
  ! written once.
  !===================================================================!

  pure subroutine through_blocks(block_of, num_blocks, num_components, fine, coarse, &
       & transposed, average)

    integer              , intent(in)           :: block_of(:)
    integer              , intent(in)           :: num_blocks, num_components
    real(dp), allocatable, intent(inout)        :: fine(:)
    real(dp), allocatable, intent(inout)        :: coarse(:)
    logical              , intent(in)           :: transposed
    logical              , intent(in), optional :: average

    integer, allocatable :: tally(:)
    integer              :: v, b, c, nv

    nv = size(block_of)
    if (any(block_of < 1) .or. any(block_of > num_blocks)) then
       error stop 'through_blocks: every entry of block_of must lie in 1..num_blocks, &
            &but at least one does not'
    end if

    if (transposed) then
       ! R^T: every member takes its block's value
       if (size(coarse) /= num_blocks * num_components) then
          error stop 'through_blocks: coarse must hold num_blocks*num_components values, &
               &but its size does not match'
       end if
       if (allocated(fine)) deallocate(fine)
       allocate(fine(nv * num_components))
       do v = 1, nv
          b = block_of(v)
          do c = 1, num_components
             fine((v - 1) * num_components + c) = coarse((b - 1) * num_components + c)
          end do
       end do
    else
       ! R: each block the sum of its members, counted
       if (size(fine) /= nv * num_components) then
          error stop 'through_blocks: fine must hold nv*num_components values, &
               &but its size does not match'
       end if
       if (allocated(coarse)) deallocate(coarse)
       allocate(coarse(num_blocks * num_components), tally(num_blocks))
       coarse = 0.0_dp
       tally  = 0
       do v = 1, nv
          b = block_of(v)
          tally(b) = tally(b) + 1
          do c = 1, num_components
             coarse((b - 1) * num_components + c) = &
                  & coarse((b - 1) * num_components + c) + fine((v - 1) * num_components + c)
          end do
       end do
       ! the mean, where requested, over the blocks that have members
       if (present(average)) then
          if (average) then
             do b = 1, num_blocks
                if (tally(b) > 0) then
                   coarse((b - 1) * num_components + 1 : b * num_components) = &
                        & coarse((b - 1) * num_components + 1 : b * num_components) / real(tally(b), dp)
                end if
             end do
          end if
       end if
    end if

  end subroutine through_blocks

end module transform_structure
