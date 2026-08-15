!=====================================================================!
! NOT YET A LEVEL . THE COMPUTATIONAL GRAPH
!
! The naming pass reserves this vocabulary; the releveling report
! proposes its seat in the tower. What lives here is the epistemic
! pair (COMPUTATIONAL-GRAPH.md):
!
!      G = ( Q, R )
!
! Q the data - ALL computational data participating in a problem,
! never presumed to be one field - and R the residual/operator that
! governs it. Each seat is realized or it is BOTTOM: unknown,
! unrealized, no occupant. Two seats, four states, one type:
!
!      void graph        (bottom, bottom)
!      data graph        (Q     , bottom)
!      operator graph    (bottom, R     )
!      realized graph    (Q     , R     )
!
! Exactly one holds, always - the classifier below is the law, and
! the suite holds it to exactly-one.
!
!                       BOTTOM IS NOT EMPTY
!
! An allocated seat whose occupant has zero entries is REALIZED:
! the knowledge is present and its content happens to be small. An
! unallocated seat is bottom: nothing is known. And void speaks of
! knowledge, never topology - a void computational graph may ride
! on a perfectly whole structural graph GAMMA = (S, P). Structural
! emptiness and epistemic unknownness are different absences.
!
!                      REALIZED IS NOT SOLVED
!
! A realized graph asserts occupancy, nothing more. R(Q) = 0 is a
! separate property with its own words - satisfied, consistent,
! converged - asserted elsewhere, never here. A realized graph may
! hold a converged solution or a deliberately inconsistent pair,
! and this module cannot tell the difference. That is by design.
!
!                     THE SCOPE, KEPT MINIMAL
!
! The type, the four state constants, the classification queries,
! the seat accessors - and NOTHING ELSE. No solving, no fitting,
! no minimization, no adjoints, no discovery. The seats hold their
! occupants as unlimited polymorphics precisely so this module
! fixes no storage representation of Q and decides nothing about
! R's relationship to the legacy graph_operation: vocabulary now,
! architecture when it is earned.
!
! The structure is BORROWED, as views borrow: the pair rides above
! a structural host it does not own, so the host must outlive the
! computational graph that names it. The computational graph is the
! fourth citizen to sign the identity roll, after the carriers, the
! relations and the structural graph - one law, one roll, four
! hands.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_state

  use graph_identity , only : token, mint_token
  use graph_structure, only : relational_graph

  implicit none

  private
  public :: computational_graph
  public :: GRAPH_STATE_VOID, GRAPH_STATE_DATA
  public :: GRAPH_STATE_OPERATOR, GRAPH_STATE_REALIZED
  public :: state_name
  public :: void_graph, data_graph, operator_graph, realized_graph

  !===================================================================!
  ! The four states: the 2 x 2 truth table of (has Q, has R). Four
  ! constants, one type - a finite classification is a value, never
  ! a hierarchy.
  !===================================================================!

  integer, parameter :: GRAPH_STATE_VOID     = 0
  integer, parameter :: GRAPH_STATE_DATA     = 1
  integer, parameter :: GRAPH_STATE_OPERATOR = 2
  integer, parameter :: GRAPH_STATE_REALIZED = 3

  !===================================================================!
  ! The computational graph: two epistemic seats, a borrowed
  ! structural host, an identity. An unallocated seat IS bottom -
  ! the language enforces the ontology.
  !===================================================================!

  type :: computational_graph

     type(token)                  , private :: identity
     character(len=:), allocatable, private :: label

     class(relational_graph), pointer, private :: host => null()

     class(*), allocatable, private :: data_seat
     class(*), allocatable, private :: residual_seat

   contains

     procedure :: has_data
     procedure :: has_operator
     procedure :: state

     procedure :: data      => data_of
     procedure :: residual  => residual_of
     procedure :: structure => structure_of

     procedure :: declare
     procedure :: id
     procedure :: same_as
     procedure :: name

  end type computational_graph

contains

  !===================================================================!
  ! The four constructors: one name per state, one type underneath.
  ! Each may name the structural host it rides on; none is obliged
  ! to, for knowledge of structure and knowledge of data are
  ! different knowledge.
  !===================================================================!

  type(computational_graph) function void_graph(name, structure) &
       & result(this)

    character(len=*)       , intent(in)                     :: name
    class(relational_graph), intent(in), target  , optional :: structure

    if (present(structure)) this % host => structure

    call this % declare(name)

  end function void_graph

  type(computational_graph) function data_graph(name, data, structure) &
       & result(this)

    character(len=*)       , intent(in)                     :: name
    class(*)               , intent(in)                     :: data
    class(relational_graph), intent(in), target  , optional :: structure

    allocate(this % data_seat, source=data)
    if (present(structure)) this % host => structure

    call this % declare(name)

  end function data_graph

  type(computational_graph) function operator_graph(name, residual, &
       & structure) result(this)

    character(len=*)       , intent(in)                     :: name
    class(*)               , intent(in)                     :: residual
    class(relational_graph), intent(in), target  , optional :: structure

    allocate(this % residual_seat, source=residual)
    if (present(structure)) this % host => structure

    call this % declare(name)

  end function operator_graph

  type(computational_graph) function realized_graph(name, data, &
       & residual, structure) result(this)

    character(len=*)       , intent(in)                     :: name
    class(*)               , intent(in)                     :: data
    class(*)               , intent(in)                     :: residual
    class(relational_graph), intent(in), target  , optional :: structure

    allocate(this % data_seat    , source=data)
    allocate(this % residual_seat, source=residual)
    if (present(structure)) this % host => structure

    call this % declare(name)

  end function realized_graph

  !===================================================================!
  ! The classification queries. Allocation status IS the epistemic
  ! status: a seat either holds realized knowledge or it is bottom,
  ! and the two questions below are the whole classification.
  !===================================================================!

  pure logical function has_data(this)

    class(computational_graph), intent(in) :: this

    has_data = allocated(this % data_seat)

  end function has_data

  pure logical function has_operator(this)

    class(computational_graph), intent(in) :: this

    has_operator = allocated(this % residual_seat)

  end function has_operator

  !===================================================================!
  ! The classifier: exactly one state, always. The branches exhaust
  ! the truth table and exclude one another by construction - the
  ! suite still holds them to it.
  !===================================================================!

  pure integer function state(this)

    class(computational_graph), intent(in) :: this

    if (allocated(this % data_seat)) then
       if (allocated(this % residual_seat)) then
          state = GRAPH_STATE_REALIZED
       else
          state = GRAPH_STATE_DATA
       end if
    else
       if (allocated(this % residual_seat)) then
          state = GRAPH_STATE_OPERATOR
       else
          state = GRAPH_STATE_VOID
       end if
    end if

  end function state

  !===================================================================!
  ! The canonical names, for diagnostics and for nobody's synonyms:
  ! no empty, no partial, no physics, no solution graph.
  !===================================================================!

  function state_name(s) result(said)

    integer, intent(in)           :: s
    character(len=:), allocatable :: said

    select case (s)
    case (GRAPH_STATE_VOID)
       said = 'void graph'
    case (GRAPH_STATE_DATA)
       said = 'data graph'
    case (GRAPH_STATE_OPERATOR)
       said = 'operator graph'
    case (GRAPH_STATE_REALIZED)
       said = 'realized graph'
    case default
       error stop 'graph_state: there is no fifth state'
    end select

  end function state_name

  !===================================================================!
  ! The seat accessors: references into owned seats, refused loudly
  ! at bottom. Bottom is not a value, and no accessor manufactures
  ! one.
  !===================================================================!

  function data_of(this) result(q)

    class(computational_graph), target, intent(in) :: this
    class(*), pointer                              :: q

    if (.not. allocated(this % data_seat)) then
       error stop 'graph_state: the data seat is bottom - unknown is not a value'
    end if

    q => this % data_seat

  end function data_of

  function residual_of(this) result(r)

    class(computational_graph), target, intent(in) :: this
    class(*), pointer                              :: r

    if (.not. allocated(this % residual_seat)) then
       error stop 'graph_state: the residual seat is bottom - unknown is not a value'
    end if

    r => this % residual_seat

  end function residual_of

  function structure_of(this) result(host)

    class(computational_graph), intent(in) :: this
    class(relational_graph), pointer       :: host

    if (.not. associated(this % host)) then
       error stop 'graph_state: this computational graph rides on no declared structure'
    end if

    host => this % host

  end function structure_of

  !===================================================================!
  ! The identity block, one law with the carriers', the relations'
  ! and the structural graph's: the fourth hand on one roll.
  !===================================================================!

  subroutine declare(this, name)

    class(computational_graph), intent(inout)        :: this
    character(len=*)          , intent(in), optional :: name

    if (this % identity % declared()) then
       error stop 'graph_state: a computational graph never signs twice'
    end if

    this % identity = mint_token()
    if (present(name)) this % label = name

  end subroutine declare

  pure type(token) function id(this)

    class(computational_graph), intent(in) :: this

    id = this % identity

  end function id

  pure logical function same_as(this, other)

    class(computational_graph), intent(in) :: this
    class(computational_graph), intent(in) :: other

    same_as = this % identity % matches(other % identity)

  end function same_as

  function name(this)

    class(computational_graph), intent(in) :: this
    character(len=:), allocatable          :: name

    if (allocated(this % label)) then
       name = this % label
    else
       name = ''
    end if

  end function name

end module graph_state
