!=====================================================================!
! The transform prime: the operation between graphs, graph -> graph
! with data mapped along. Two symbols at this level, both
! admissibility predicates: may this transform act on that graph, and
! on that data defined on it. The concrete types - partition,
! assemble, coarsen, refine - are each checked by a round-trip law:
!
!      exact        assemble(partition(G)) = G      both ways
!      one-sided    coarsen(refine(G)) = G          one way only
!
! The interfaces are impure for a historical language constraint that
! no longer binds (F2018 C1594 on copying a set that contains pointers
! graph inside a pure subprogram); an impure interface permits a pure
! override, and every implementation is pure.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module transform_structure

  use util_precision, only : dp
  use view_directed , only : directed_graph
  use field_calculus, only : field

  implicit none

  private

  public :: transform
  public :: through_blocks

  type, abstract :: transform

   contains

     procedure(transform_on_graph_interface), deferred :: defined_on_graph
     procedure(transform_on_data_interface) , deferred :: defined_on_data

  end type transform

  abstract interface

     logical function transform_on_graph_interface(this, input_graph)
       import :: transform, directed_graph
       class(transform), intent(in) :: this
       class(directed_graph), intent(in) :: input_graph
     end function transform_on_graph_interface

     logical function transform_on_data_interface(this, input_graph, input_data)
       import :: transform, directed_graph, field
       class(transform), intent(in) :: this
       class(directed_graph), intent(in) :: input_graph
       class(field), intent(in) :: input_data
     end function transform_on_data_interface

  end interface

contains

  !===================================================================!
  ! The one map between a level and its blocks, read both ways. Read
  ! forward it is restriction R: each block the sum of its members,
  ! or their mean where average is requested. Read backward it is R^T,
  ! every member taking its block's value - the injected prolongation.
  ! Values are num_components wide per member. coarsen(refine(G)) = G
  ! is this map read twice, which is why it is written once.
  !===================================================================!

  pure subroutine through_blocks(block_of, num_blocks, num_components, fine, coarse, &
       & transposed, average)

    integer              , intent(in)    :: block_of(:)
    integer              , intent(in)    :: num_blocks, num_components
    real(dp), allocatable, intent(inout) :: fine(:)
    real(dp), allocatable, intent(inout) :: coarse(:)
    logical              , intent(in)    :: transposed
    logical              , intent(in), optional :: average

    integer, allocatable :: tally(:)
    integer :: v, b, c, nv

    nv = size(block_of)
    if (any(block_of < 1) .or. any(block_of > num_blocks)) then
       error stop 'through_blocks: every entry of block_of must lie in 1..num_blocks, &
            &but at least one does not'
    end if

    if (transposed) then
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
