!=====================================================================!
! THE PARTITIONED SHIFTED-LAPLACIAN FIXTURE - test-local, and the
! sealing object of this tower. It represents the SAME global map
!
!      A : V(G) -> V(G),      A(q) = 2q - L(q)
!
! but it never evaluates that map globally. Every call decomposes,
! acts locally on each part's own topology, and reassembles:
!
!      GLOBAL q
!         |
!         +-- refresh G1 overlap -> q1 -> A_G1 -> owned y1 --+
!         |                                                  |
!         +-- refresh G2 overlap -> q2 -> A_G2 -> owned y2 --+
!                                                            v
!                                                    assemble + sum
!                                                            v
!                                                       GLOBAL A q
!
! THREE THINGS ARE KEPT DISTINCT, and the type's shape enforces it:
!
!      STRUCTURAL PARTITION   G -> {G1,G2}    ONCE, at construction
!      NUMERICAL OVERLAP      q -> {q1,q2}    EVERY apply
!      OWNED ASSEMBLY         {y1,y2} -> y    EVERY apply
!
! So the type owns STRUCTURE - the global graph, both parts, the
! partitioners that cut them and the assembler that gathers them -
! and owns NO mutable numerical state whatever: no cached q, no
! cached overlap, no previous result. Nothing in it can go stale,
! because there is nothing in it to go stale.
!
! And unlike Gate B's shifted_laplacian, which is GRAPH-PARAMETERIZED
! and will act on any graph it is handed, this composite is
! DECOMPOSITION-CONTEXT-BOUND: G1 and G2 were cut from one
! particular G, so any other graph - however many vertices it
! happens to have - is refused, in domain() as well as apply(), at
! the earliest honest contract point.
!
! The composite never calls the Laplacian on the whole graph. The
! only numerical topology actions in this file are on G1 and G2.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module partitioned_shifted_laplacian_fixture

  use iso_fortran_env  , only : dp => REAL64
  use fractal_graph      , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation
  use graph_set_map      , only : set_map
  use graph_label_map    , only : label_map
  use graph_inclusion_map, only : inclusion_map
  use graph_operation_view, only : graph_operation
  use graph_directed_view, only : directed_graph
  use graph_field_calculus, only : graph_field
  use class_graph      , only : directed_stored_graph
  use class_graph_field, only : field
  use class_graph_partitioner, only : partitioner, PARTITION_LINEAR
  use class_graph_assembler  , only : assembler
  use shifted_laplacian_fixture, only : shifted_laplacian

  use graph_partition_relation, only : partition_relation
  implicit none

  private
  public :: partitioned_shifted_laplacian

  type, extends(graph_operation) :: partitioned_shifted_laplacian

     ! STRUCTURE, cut once and never rebuilt.
     type(directed_stored_graph)        :: whole
     class(directed_graph), allocatable :: g1, g2
     type(partitioner)         :: p1, p2
     type(assembler)           :: asm

     ! ONE RELATION PER PART, because there is one per part. A single
     ! relation here would hold whichever cut ran last and carry the
     ! wrong numbering onto the other piece.
     type(partition_relation) :: r1, r2
     type(shifted_laplacian)   :: local

     ! ...and deliberately NO numerical state: no q, no overlap,
     ! no residual, no previous answer.

   contains
     procedure :: name   => part_name
     procedure :: domain => part_domain
     procedure :: apply  => part_apply
  end type partitioned_shifted_laplacian

  interface partitioned_shifted_laplacian
     module procedure create_partitioned
  end interface partitioned_shifted_laplacian

contains

  !===================================================================!
  ! Construction is where the STRUCTURAL partition happens - once.
  !===================================================================!

  type(partitioned_shifted_laplacian) function create_partitioned(g) &
       & result(this)

    type(directed_stored_graph), intent(in) :: g

    this % whole = g
    this % p1 = partitioner(PARTITION_LINEAR, nparts=2, part=1)
    this % p2 = partitioner(PARTITION_LINEAR, nparts=2, part=2)
    call this % p1 % partition_graph(g, this % g1, this % r1)
    call this % p2 % partition_graph(g, this % g2, this % r2)

  end function create_partitioned

  pure function part_name(this) result(name)
    class(partitioned_shifted_laplacian), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'partitioned shifted laplacian'
  end function part_name

  !===================================================================!
  ! The composite answers on the decomposition it was built from -
  ! and refuses to answer about any other graph. This is the
  ! earliest honest contract point: a solver attaching on a foreign
  ! host dies here rather than deep inside a matvec.
  !===================================================================!

  subroutine part_domain(this, input_graph, domain, nentries)

    class(partitioned_shifted_laplacian), intent(in) :: this
    class(directed_graph)   , intent(in)  :: input_graph
    type(set_graph), intent(out) :: domain
    integer        , intent(out) :: nentries

    call demand_the_recorded_context(this, input_graph)

    domain   = this % whole % vertex_set()
    nentries = this % whole % num_vertices()

  end subroutine part_domain

  !===================================================================!
  ! One matvec: refresh both overlaps from the CURRENT state, act
  ! locally on each part's topology, keep only what each part owns,
  ! and sum. Nothing global is differentiated; nothing is cached.
  !===================================================================!

  subroutine part_apply(this, input_graph, input_data, output)

    class(partitioned_shifted_laplacian), intent(in) :: this
    class(directed_graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(field)                     :: out
    class(graph_field), allocatable :: q1, q2, a1, a2
    type(set_graph)                 :: dom
    real(dp), allocatable           :: total(:)

    !----------------------------------------------------------------!
    ! The interpretation environment of this call. Everything carved
    ! below is consumed below, so it lives and dies here and no map
    ! enters the operation's own type.
    !----------------------------------------------------------------!

    type(set_map)                   :: sets
    type(label_map)                 :: labels
    type(inclusion_map)             :: inclusions

    call demand_the_recorded_context(this, input_graph)

    if (.not. present(input_data)) then
       error stop 'partitioned action: the action needs a state to read'
    end if
    if (size(input_data) /= 1) then
       error stop 'partitioned action: the action reads exactly one state'
    end if
    dom = input_data(1) % domain()
    if (.not. dom % same_as(this % whole % vertex_set())) then
       error stop 'partitioned action: the state must live on the global vertex carrier'
    end if

    !----------------------------------------------------------------!
    ! The transports CARVE. Everything carved here is consumed here,
    ! so the interpretation is local to the call and no map enters
    ! this operation's type - which holds an identity and a count and
    ! nothing else, like every other action in the tower.
    !----------------------------------------------------------------!

    call sets % bind(this % whole % vertex_set(), &
         & counted_set_representation(this % whole % num_vertices()))
    call bind_part(sets, this % g1)
    call bind_part(sets, this % g2)

    ! -- NUMERICAL OVERLAP REFRESH, from the state handed in NOW
    call this % p1 % partition_data(this % r1, this % whole, &
         & input_data(1), this % g1, sets, labels, inclusions, q1)
    call this % p2 % partition_data(this % r2, this % whole, &
         & input_data(1), this % g2, sets, labels, inclusions, q2)

    ! -- LOCAL TOPOLOGY ACTIONS, each on its own part
    call act_locally(this % local, this % g1, q1, a1)
    call act_locally(this % local, this % g2, q2, a2)

    ! -- OWNED ASSEMBLY, and the sum of the two contributions
    allocate(total(this % whole % num_vertices()))
    total = 0.0_dp
    call add_owned(assembler(), this % r1, this % g1, a1, this % whole, &
         & sets, labels, inclusions, total)
    call add_owned(assembler(), this % r2, this % g2, a2, this % whole, &
         & sets, labels, inclusions, total)

    out = field('partitioned action', this % whole % vertex_set(), &
         & this % whole % num_vertices())
    call out % set_real_vector(total)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine part_apply

  !===================================================================!
  ! The context guard: same cardinality is irrelevant, identity is
  ! everything. G1 and G2 were cut from one particular G.
  !===================================================================!

  subroutine demand_the_recorded_context(this, input_graph)

    class(partitioned_shifted_laplacian), intent(in) :: this
    class(directed_graph)                        , intent(in) :: input_graph

    type(set_graph) :: given, recorded

    given    = input_graph % vertex_set()
    recorded = this % whole % vertex_set()

    if (.not. given % same_as(recorded)) then
       error stop 'partitioned action: this decomposition belongs to another graph'
    end if

  end subroutine demand_the_recorded_context

  !===================================================================!
  ! Apply the local action on ONE part. The repack is only Fortran
  ! bookkeeping - an array constructor needs a concrete type - and
  ! the values are unchanged.
  !===================================================================!

  subroutine act_locally(local, part, q_part, answer)

    type(shifted_laplacian), intent(in)  :: local
    class(directed_graph)           , intent(in)  :: part
    class(graph_field)     , intent(in)  :: q_part
    class(graph_field), allocatable, intent(out) :: answer

    type(field)                     :: qf
    class(graph_field), allocatable :: out
    real(dp), allocatable           :: v(:)

    call q_part % get_real_vector(v)
    qf = field('local state', part % vertex_set(), part % num_vertices())
    call qf % set_real_vector(v)

    call local % apply(part, [qf], out)
    allocate(answer, source=out)

  end subroutine act_locally

  !===================================================================!
  ! Assemble ONE part's local answer home and add its contribution.
  ! The assembler keeps owned members only; whether it hands back a
  ! full global field or a subobject of the global carrier, every
  ! value is placed by MEMBER.
  !===================================================================!

  subroutine add_owned(asm, rel, part, answer, whole, sets, labels, &
       & inclusions, total)

    type(assembler)   , intent(in)    :: asm

    ! The part's OWN relation, travelling beside the part. The two
    ! arrive together or the wrong numbering goes home.
    type(partition_relation), intent(in) :: rel

    class(directed_graph)      , intent(in)    :: part
    class(graph_field), intent(in)    :: answer
    type(directed_stored_graph), intent(in)    :: whole
    type(set_map)     , intent(inout) :: sets
    type(label_map)   , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    real(dp)          , intent(inout) :: total(:)

    class(graph_field), allocatable :: home
    type(set_graph)                 :: dom
    real(dp), allocatable           :: v(:)
    integer , allocatable           :: mem(:)
    integer                         :: i

    call asm % assemble_data(rel, part, answer, whole, sets, labels, &
         & inclusions, home)
    dom = home % domain()
    call home % get_real_vector(v)

    !----------------------------------------------------------------!
    ! The select type asked whether the assembled domain was a subset
    ! TYPE, and wrote the same arithmetic twice. It is one branch now,
    ! and not by approximation: for a counted domain members are
    ! 1..n and local_index is the identity, so the subset arm REDUCES
    ! to the full arm exactly.
    !----------------------------------------------------------------!

    call sets % members_of(dom, mem)
    do i = 1, size(mem)
       total(mem(i)) = total(mem(i)) + v(sets % index_in(dom, mem(i)))
    end do

  end subroutine add_owned

  !===================================================================!
  ! A part is a new graph, so its carriers are new domains and must be
  ! described before a field can be seated on them.
  !===================================================================!

  subroutine bind_part(sets, g)

    type(set_map), intent(inout) :: sets
    class(directed_graph) , intent(in)    :: g

    call sets % bind(g % vertex_set(), &
         & counted_set_representation(g % num_vertices()))
    call sets % bind(g % edge_set(), &
         & counted_set_representation(g % num_edges()))

  end subroutine bind_part

end module partitioned_shifted_laplacian_fixture
