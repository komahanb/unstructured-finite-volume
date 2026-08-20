!=====================================================================!
! The operation prime: the verb within a graph, (graph, fields) ->
! field. Three symbols. name says what it is; domain says which
! member set the result lives on and how many entries it has; apply
! does the work.
!
! apply writes its result into the output argument and never adds to
! what was there. The argument is intent(inout) only so a caller
! already holding a buffer of the right shape can lend it and save
! an allocation; lending changes the cost of the call, not its
! meaning.
!
! A concrete operation is handed the fields it reads when it is
! constructed - a coefficient, a measure, a geometry field arrives as
! an argument the compiler checks - so apply fetches nothing by name.
!
! Every operation also reports max_degree, the highest order of
! exact partial action it computes (0 unless overridden), and
! partial_action, one mixed partial - differentiated once per entry
! of slots(:), contracted against the matching direction field, on
! its own domain. The default stops the program, since an operation
! with max_degree 0 declares no partials.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_action

  use view_directed , only : directed_graph
  use graph_fractal       , only : graph
  use field_calculus, only : field

  implicit none

  private

  public :: operation

  type, abstract :: operation

   contains

     procedure(operation_name_interface)  , deferred :: name
     procedure(operation_domain_interface), deferred :: domain
     procedure(operation_apply_interface) , deferred :: apply

     procedure :: max_degree     => operation_max_degree
     procedure :: partial_action => operation_partial_action

  end type operation

  abstract interface

     pure function operation_name_interface(this) result(name)
       import :: operation
       class(operation), intent(in) :: this
       character(len=:), allocatable :: name
     end function operation_name_interface

     !---------------------------------------------------------------!
     ! Where the result lives: WHICH set, and HOW MANY entries it
     ! has. The count travels beside the identity because every
     ! caller wants exactly those two things - to check the domain
     ! matches, and to size a field.
     !---------------------------------------------------------------!

     subroutine operation_domain_interface(this, input_graph, domain, &
          & num_entries)
       import :: operation, directed_graph, graph
       class(operation), intent(in)  :: this
       class(directed_graph)          , intent(in)  :: input_graph
       type(graph)       , intent(out) :: domain
       integer               , intent(out) :: num_entries
     end subroutine operation_domain_interface

     subroutine operation_apply_interface(this, input_graph, input_data, output)
       import :: operation, directed_graph, field
       class(operation), intent(in) :: this
       class(directed_graph), intent(in) :: input_graph
       class(field), intent(in), optional :: input_data(:)
       class(field), allocatable, intent(inout) :: output
     end subroutine operation_apply_interface

  end interface

contains

  !===================================================================!
  ! The default: no exact partial action of any order.
  !===================================================================!

  pure function operation_max_degree(this) result(degree)

    class(operation), intent(in) :: this
    integer :: degree

    associate (u1 => this); end associate

    degree = 0

  end function operation_max_degree

  !===================================================================!
  ! The default refuses every request, because max_degree is 0. A
  ! concrete type that declares a positive max_degree overrides both
  ! bindings; the order requested must not exceed its max_degree.
  !===================================================================!

  subroutine operation_partial_action(this, input_graph, input_data, &
       & slots, directions, output)

    class(operation), intent(in)             :: this
    class(directed_graph), intent(in)        :: input_graph
    class(field), intent(in)                 :: input_data(:)
    integer, intent(in)                      :: slots(:)
    class(field), intent(in)                 :: directions(:)
    class(field), allocatable, intent(inout) :: output

    associate (u1 => this, u2 => input_graph, u3 => input_data, &
         & u4 => slots, u5 => directions); end associate
    if (allocated(output)) deallocate(output)

    error stop 'operation: the requested order is within max_degree'

  end subroutine operation_partial_action

end module operation_action
