!=====================================================================!
! Things that act on a graph.
!
! An operation reads a graph and some values on it and answers more
! values. Four specializations are declared here because each binds
! one thing and defers the rest: a reduction answers on a one-member
! domain, a broadcast reads from one, a discretization is built from
! another operation by binding it to a graph, and a linearization is
! built from another by freezing a state.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_graph_operation

  use iso_fortran_env, only : dp => REAL64
  use structure_graph, only : graph
  use data_graph_field, only : graph_field, graph_functional

  implicit none

  private

  public :: graph_operation
  public :: graph_reduction
  public :: graph_broadcast
  public :: graph_discretization
  public :: graph_linearization

  !===================================================================!
  ! GRAPH_OPERATION. The verb within a graph: data in, data out.
  !
  !      input_graph, input_data(:)  --- apply --->  output
  !
  ! Three symbols. name says what it is. domain says where the answer
  ! lives: a member set of the input graph. apply does the work,
  ! under the lent-buffer rule of the header.
  !
  ! A concrete operation is handed the fields it reads when it is
  ! constructed - a coefficient, a measure, a geometry field arrives
  ! as an argument the compiler checks, and apply fetches nothing by
  ! name. This is what keeps the call structure visible from
  ! construction to result.
  !
  ! The concretions of this role, arriving on level 1, split by
  ! order: first-order kernels act on fields (a differential
  ! operator, a balance), and the one higher-order citizen acts on
  ! operators (the walk, a traversal whose kernel is not yet bound).
  ! The maps that touch the one-entry domain - the reduction and the
  ! broadcast - stand beside this role rather than under it, one
  ! deliberate deviation, recorded where they are declared.
  !===================================================================!

  type, abstract :: graph_operation

   contains

     procedure(operation_name_interface)  , deferred :: name
     procedure(operation_domain_interface), deferred :: domain
     procedure(operation_apply_interface) , deferred :: apply

  end type graph_operation

  !===================================================================!
  ! GRAPH_REDUCTION. Many values become one: field -> functional.
  !
  ! The whole four-step dance, and why it is four steps:
  !
  !   part 1  [q q q]  --accumulate-->  (sum 6, count 3) --+
  !                                                        +--> combine
  !   part 2  [q q]    --accumulate-->  (sum 14, count 2) -+       |
  !                                                                v
  !                                                   (sum 20, count 5)
  !                                                                |
  !                                                            finalize
  !                                                                |
  !                                                                v
  !                                                           J = 4.0
  !
  ! Means of 2 and 7 do not average to 4.5. They average to 4,
  ! because the sum and the count travel together and the division
  ! happens once, at the very end. A reduction that finished early
  ! on each part would get this wrong.
  !
  ! The measure argument is what turns a bare sum into an integral -
  ! weight each cell by its volume, or each face by its area, and
  ! the answer stops depending on how finely the mesh was cut. The
  ! measure seat is also the inner product's second field:
  !
  !      sum       J = sum q_i
  !      integral  J = sum q_i V_i          <- measure is the volume
  !      average   J = sum q_i V_i / sum V_i
  !      norm      J = sqrt( sum q_i^2 V_i )
  !      product   J = sum u_i v_i          <- measure is the field v
  !===================================================================!

  type, abstract, extends(graph_operation) :: graph_reduction

   contains

     ! The four steps, for a caller that owns the parallel dance.
     procedure(reduction_initialize_interface), deferred :: initialize
     procedure(reduction_accumulate_interface), deferred :: accumulate
     procedure(reduction_combine_interface)   , deferred :: combine
     procedure(reduction_finalize_interface)  , deferred :: finalize

     ! All four in one call, for a caller holding the whole thing.
     ! The inherited apply is this verb's operation face: the field
     ! reduced over its own domain, the measure riding as the second
     ! input field, the functional leaving as the output - for a
     ! functional IS a field, one entry wide.
     procedure(reduction_reduce_interface), deferred :: reduce

  end type graph_reduction

  !===================================================================!
  ! GRAPH_BROADCAST. One value becomes many: functional -> field.
  ! The transpose of the reduction, and the round trip is the law:
  !
  !      average( copy( c ) )  = c
  !      sum(     share( c ) ) = c
  !
  ! One step, not four. A reduction walks parts and must combine
  ! their partial answers; a broadcast writes the same fill on every
  ! part, so there is nothing to communicate and nothing to stage.
  !===================================================================!

  type, abstract, extends(graph_operation) :: graph_broadcast

   contains

     procedure(broadcast_interface), deferred :: broadcast

  end type graph_broadcast

  !===================================================================!
  ! DISCRETIZATION_OPERATOR. An operation built from an operation by
  ! binding it to a graph's arithmetic - the act that turns pde and
  ! ode into algebra. A scheme is a MOTIF stamped along the domain
  ! graph:
  !
  !      backward euler    o<--o          reach 1, weights [1, -1]
  !      bdf-k             o<--o<-..<--o  reach k, the k-step table
  !      two-point flux    one mesh edge, weights -+ 1/delta
  !      fitted, degree p  the p-ring, weights solved by the fit
  !
  ! What every member owes by contract: its dependency PATTERN, as a
  ! graph. The support is to a field what this pattern is to a
  ! derived operator - values sit on members, arithmetic flows on
  ! pairs. The minimizers one level up interrogate the pattern - the
  ! diagonal, the colouring, the triangularity, the Galerkin road -
  ! so it is exposed by law, never by inspection.
  !===================================================================!

  type, abstract, extends(graph_operation) :: graph_discretization

   contains

     procedure(discretization_pattern_interface), deferred :: dependencies

  end type graph_discretization

  !===================================================================!
  ! LINEARIZATION_OPERATOR. An operation built from an operation by
  ! freezing a state: the tangent of S at q, wearing the operation
  ! face, so a governed minimizer sees an ordinary linear question,
  !
  !      J v  at  q        the derivative of S along v, at the
  !                        standing state
  !
  ! One deferred verb beyond the face: freeze, which moves the
  ! standing state (and may cache the base residual) between
  ! newton's steps. The difference road divides two residuals; the
  ! exact road, when a statement speaks its own tangent, joins as a
  ! second concretion with no change here.
  !===================================================================!

  type, abstract, extends(graph_operation) :: graph_linearization

   contains

     procedure(linearization_freeze_interface), deferred :: freeze

  end type graph_linearization

  abstract interface

     pure function operation_name_interface(this) result(name)
       import :: graph_operation
       class(graph_operation), intent(in) :: this
       character(len=:), allocatable :: name
     end function operation_name_interface

     subroutine operation_domain_interface(this, input_graph, domain)
       import :: graph_operation, graph
       class(graph_operation), intent(in) :: this
       class(graph), intent(in) :: input_graph
       class(graph), allocatable, intent(out) :: domain
     end subroutine operation_domain_interface

     subroutine operation_apply_interface(this, input_graph, input_data, output)
       import :: graph_operation, graph, graph_field
       class(graph_operation), intent(in) :: this
       class(graph), intent(in) :: input_graph
       class(graph_field), intent(in), optional :: input_data(:)
       class(graph_field), allocatable, intent(inout) :: output
     end subroutine operation_apply_interface

     pure subroutine reduction_initialize_interface(this, state)
       import :: graph_reduction, graph_functional
       class(graph_reduction), intent(in) :: this
       class(graph_functional), allocatable, intent(inout) :: state
     end subroutine reduction_initialize_interface

     pure subroutine reduction_accumulate_interface(this, field, state, measure)
       import :: graph_reduction, graph_field, graph_functional
       class(graph_reduction), intent(in) :: this
       class(graph_field), intent(in) :: field
       class(graph_functional), intent(inout) :: state
       class(graph_field), intent(in), optional :: measure
     end subroutine reduction_accumulate_interface

     pure subroutine reduction_combine_interface(this, left, right, combined)
       import :: graph_reduction, graph_functional
       class(graph_reduction), intent(in) :: this
       class(graph_functional), intent(in) :: left
       class(graph_functional), intent(in) :: right
       class(graph_functional), allocatable, intent(inout) :: combined
     end subroutine reduction_combine_interface

     pure subroutine reduction_finalize_interface(this, state, functional)
       import :: graph_reduction, graph_functional
       class(graph_reduction), intent(in) :: this
       class(graph_functional), intent(in) :: state
       class(graph_functional), allocatable, intent(inout) :: functional
     end subroutine reduction_finalize_interface

     subroutine reduction_reduce_interface(this, field, functional, measure)
       import :: graph_reduction, graph_field, graph_functional
       class(graph_reduction), intent(in) :: this
       class(graph_field), intent(in) :: field
       class(graph_functional), allocatable, intent(inout) :: functional
       class(graph_field), intent(in), optional :: measure
     end subroutine reduction_reduce_interface

     pure subroutine broadcast_interface(this, functional, field)
       import :: graph_broadcast, graph_functional, graph_field
       class(graph_broadcast) , intent(in)    :: this
       class(graph_functional), intent(in)    :: functional
       class(graph_field)     , intent(inout) :: field
     end subroutine broadcast_interface

     subroutine discretization_pattern_interface(this, pattern)
       import :: graph_discretization, graph
       class(graph_discretization), intent(in) :: this
       class(graph), allocatable, intent(out)     :: pattern
     end subroutine discretization_pattern_interface

     subroutine linearization_freeze_interface(this, at, base)
       import :: graph_linearization, dp
       class(graph_linearization), intent(inout) :: this
       real(dp), intent(in)           :: at(:)
       real(dp), intent(in), optional :: base(:)
     end subroutine linearization_freeze_interface

  end interface

end module operation_graph_operation
