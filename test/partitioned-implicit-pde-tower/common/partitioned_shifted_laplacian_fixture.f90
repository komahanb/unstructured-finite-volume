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
  use graph_carrier    , only : counted_set, subset_set, member_set
  use graph_grammar    , only : graph, graph_field, graph_operation
  use class_graph      , only : stored_graph
  use class_graph_field, only : field
  use class_graph_partitioner, only : partitioner, PARTITION_LINEAR
  use class_graph_assembler  , only : assembler
  use shifted_laplacian_fixture, only : shifted_laplacian

  implicit none

  private
  public :: partitioned_shifted_laplacian

  type, extends(graph_operation) :: partitioned_shifted_laplacian

     ! STRUCTURE, cut once and never rebuilt.
     type(stored_graph)        :: whole
     class(graph), allocatable :: g1, g2
     type(partitioner)         :: p1, p2
     type(assembler)           :: asm
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

    type(stored_graph), intent(in) :: g

    this % whole = g
    this % p1 = partitioner(PARTITION_LINEAR, nparts=2, part=1)
    this % p2 = partitioner(PARTITION_LINEAR, nparts=2, part=2)
    call this % p1 % partition_graph(g, this % g1)
    call this % p2 % partition_graph(g, this % g2)
    this % asm = assembler()

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

  subroutine part_domain(this, input_graph, domain)

    class(partitioned_shifted_laplacian), intent(in) :: this
    class(graph), intent(in) :: input_graph
    class(member_set), allocatable, intent(out) :: domain

    call demand_the_recorded_context(this, input_graph)
    allocate(domain, source=this % whole % vertex_set())

  end subroutine part_domain

  !===================================================================!
  ! One matvec: refresh both overlaps from the CURRENT state, act
  ! locally on each part's topology, keep only what each part owns,
  ! and sum. Nothing global is differentiated; nothing is cached.
  !===================================================================!

  subroutine part_apply(this, input_graph, input_data, output)

    class(partitioned_shifted_laplacian), intent(in) :: this
    class(graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(field)                     :: out
    class(graph_field), allocatable :: q1, q2, a1, a2
    class(member_set), allocatable  :: dom
    real(dp), allocatable           :: total(:)

    call demand_the_recorded_context(this, input_graph)

    if (.not. present(input_data)) then
       error stop 'partitioned action: the action needs a state to read'
    end if
    if (size(input_data) /= 1) then
       error stop 'partitioned action: the action reads exactly one state'
    end if
    call input_data(1) % domain(dom)
    if (.not. dom % same_as(this % whole % vertex_set())) then
       error stop 'partitioned action: the state must live on the global vertex carrier'
    end if

    ! -- NUMERICAL OVERLAP REFRESH, from the state handed in NOW
    call this % p1 % partition_data(this % whole, input_data(1), &
         & this % g1, q1)
    call this % p2 % partition_data(this % whole, input_data(1), &
         & this % g2, q2)

    ! -- LOCAL TOPOLOGY ACTIONS, each on its own part
    call act_locally(this % local, this % g1, q1, a1)
    call act_locally(this % local, this % g2, q2, a2)

    ! -- OWNED ASSEMBLY, and the sum of the two contributions
    allocate(total(this % whole % num_vertices()))
    total = 0.0_dp
    call add_owned(this % asm, this % g1, a1, this % whole, total)
    call add_owned(this % asm, this % g2, a2, this % whole, total)

    out = field('partitioned action', this % whole % vertex_set())
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
    class(graph)                        , intent(in) :: input_graph

    type(counted_set) :: given, recorded

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
    class(graph)           , intent(in)  :: part
    class(graph_field)     , intent(in)  :: q_part
    class(graph_field), allocatable, intent(out) :: answer

    type(field)                     :: qf
    class(graph_field), allocatable :: out
    real(dp), allocatable           :: v(:)

    call q_part % get_real_vector(v)
    qf = field('local state', part % vertex_set())
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

  subroutine add_owned(asm, part, answer, whole, total)

    type(assembler)   , intent(in)    :: asm
    class(graph)      , intent(in)    :: part
    class(graph_field), intent(in)    :: answer
    type(stored_graph), intent(in)    :: whole
    real(dp)          , intent(inout) :: total(:)

    class(graph_field), allocatable :: home
    class(member_set), allocatable  :: dom
    real(dp), allocatable           :: v(:)
    integer , allocatable           :: mem(:)
    integer                         :: i

    call asm % assemble_data(part, answer, whole, home)
    call home % domain(dom)
    call home % get_real_vector(v)

    select type (dom)
    type is (subset_set)
       call dom % members(mem)
       do i = 1, size(mem)
          total(mem(i)) = total(mem(i)) + v(dom % local_index(mem(i)))
       end do
    class default
       do i = 1, size(total)
          total(i) = total(i) + v(i)
       end do
    end select

  end subroutine add_owned

end module partitioned_shifted_laplacian_fixture
