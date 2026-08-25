! The relations a coupling carries, and the weights on them.
!
! One block of five instants under a backward difference of order
! two. Its coupling holds
!
!   branch(1)   the five slices, then two more carriers: the set of
!               state components and the set of constraint instances
!   branch(2)   one relation, from the first of those sets to the
!               second, whose tuples are the scheme reach
!
! The two extra carriers are what the relation's domains refer to.
! The level's own check requires only that the carriers begin with
! its members, so carriers beyond them are permitted, and the
! relational check requires that every domain of every relation is
! one of the carriers. The two checks meet on this spine without
! either being weakened.
!
! One relation serves the whole block. A relation for each edge would
! be one object per instant per row, which is the count view_set
! exists to keep at O(1).
!
! The weights are then printed against the relation's own tuples, so
! that the sparsity and the numbers are shown to be one object read
! two ways.
program coupling_relation

  use util_precision  , only : dp
  use graph_fractal         , only : graph
  use view_level            , only : level_storage, level_consistent, &
       & level_num_members
  use view_relational       , only : relational_binding, num_member_sets, &
       & num_relations, relational_valid, relation_at
  use relation_finitary     , only : relation
  use relation_binary       , only : csr_relation
  use map_set               , only : set_map
  use map_set_representation, only : counted_set_representation
  use view_directed_stored  , only : stored_directed_graph
  use field_calculus        , only : field
  use field_stored          , only : stored_field
  use operation_family_bdf  , only : bdf_family
  use operation_weight      , only : scheme_weight

  implicit none

  integer , parameter :: order = 2
  integer , parameter :: num_instants = 5
  integer , parameter :: num_conditions = 2

  type(level_storage)      :: store
  type(relational_binding) :: binding
  type(set_map)            :: sets
  type(csr_relation)       :: reach

  integer, allocatable :: tails(:), heads(:), source_degree(:), determines(:)
  integer, allocatable :: table(:,:)
  real(dp), allocatable :: weight(:), dt(:)
  integer :: slices(num_instants), source_carrier, target_carrier
  integer :: relation_element, coupling, block
  integer :: j, k

  !--------------------------------------------------------------------!
  ! The scheme reach: every velocity row that fits, then every
  ! acceleration row that fits.
  !--------------------------------------------------------------------!

  tails = [((k - j, j = 0, order), k = order + 1, num_instants), &
       &   ((k - j, j = 0, 2 * order), k = 2 * order + 1, num_instants)]
  heads = [((k, j = 0, order), k = order + 1, num_instants), &
       &   ((k, j = 0, 2 * order), k = 2 * order + 1, num_instants)]
  determines = [((1, j = 0, order), k = order + 1, num_instants), &
       &        ((2, j = 0, 2 * order), k = 2 * order + 1, num_instants)]
  allocate(source_degree(size(tails)), source=0)

  !--------------------------------------------------------------------!
  ! The carriers, the relation element, and the coupling over them.
  !--------------------------------------------------------------------!

  do k = 1, num_instants
     slices(k) = store % assemble([integer ::], 0)
  end do
  source_carrier   = store % assemble([integer ::], 0)
  target_carrier   = store % assemble([integer ::], 0)
  relation_element = store % assemble([integer ::], 0)

  coupling = store % couple([slices, source_carrier, target_carrier], &
       & [relation_element])
  block    = store % assemble(slices, coupling)

  !--------------------------------------------------------------------!
  ! The two sets the relation runs between: one component per slice,
  ! and one constraint instance per slice per condition.
  !--------------------------------------------------------------------!

  call describe(source_carrier, num_instants)
  call describe(target_carrier, num_instants * num_conditions)

  allocate(table(2, size(tails)))
  table(1,:) = tails
  table(2,:) = (heads - 1) * num_conditions + determines

  reach = built_reach(table)

  do k = 1, num_instants
     call bind_carrier(slices(k))
  end do
  call bind_carrier(source_carrier)
  call bind_carrier(target_carrier)
  call bind_reach()

  !--------------------------------------------------------------------!
  ! What the two checks say.
  !--------------------------------------------------------------------!

  write(*,'(a)')      ' the block and its coupling'
  write(*,'(a,i3)')   '   members of the block        ', level_num_members(store % node(block))
  write(*,'(a,l3)')   '   block is consistent         ', level_consistent(store % node(block))
  write(*,'(a,i3)')   '   carriers of the coupling    ', num_member_sets(store % node(coupling))
  write(*,'(a,i3)')   '   relations of the coupling   ', num_relations(store % node(coupling))
  write(*,'(a,l3)')   '   coupling is relationally valid', &
       & relational_valid(store % node(coupling), binding)

  call show_tuples()

contains

  !-------------------------------------------------------------------!
  ! Record the extent of one carrier, which is what lets a relation
  ! number its own rows.
  !-------------------------------------------------------------------!

  subroutine describe(at, n)

    integer, intent(in) :: at, n

    type(graph), pointer :: g

    g => store % node(at)
    call sets % bind(g, counted_set_representation(n))

  end subroutine describe

  !-------------------------------------------------------------------!
  ! The relation over the two carriers. Built here so that the two
  ! graphs are named through pointers, which is what keeps the
  ! signature referring to them and not to a copy.
  !-------------------------------------------------------------------!

  function built_reach(tuples) result(r)

    integer, intent(in) :: tuples(:,:)
    type(csr_relation) :: r

    type(graph), pointer :: from, into

    from => store % node(source_carrier)
    into => store % node(target_carrier)
    r = csr_relation('scheme reach', from, into, tuples, sets)

  end function built_reach

  subroutine bind_carrier(at)

    integer, intent(in) :: at

    type(graph), pointer :: g

    g => store % node(at)
    call binding % bind_set(g, g)

  end subroutine bind_carrier

  subroutine bind_reach()

    type(graph), pointer :: g

    g => store % node(relation_element)
    call binding % bind_relation(g, reach)

  end subroutine bind_reach

  !-------------------------------------------------------------------!
  ! The weights on the same edges, printed against the relation's own
  ! tuples. The relation numbers its rows itself, so each tuple is
  ! matched to its edge rather than assumed to be in the same order.
  !-------------------------------------------------------------------!

  subroutine show_tuples()

    class(relation), pointer :: r
    integer, allocatable :: held(:,:)
    integer :: i, e, target_index

    dt = [0.0_dp, 0.30_dp, 0.20_dp, 0.40_dp, 0.25_dp]
    call edge_weights(dt, weight)

    r => relation_at(store % node(coupling), binding, 1)

    write(*,'(a)')    ' '
    write(*,'(a,i3)') ' tuples the relation holds     ', r % num_tuples()
    write(*,'(a)')    '   component   constraint   instant  determines      weight'

    call r % tuples(held)

    do i = 1, size(held, 2)
       do e = 1, size(tails)
          target_index = (heads(e) - 1) * num_conditions + determines(e)
          if (tails(e) == held(1, i) .and. target_index == held(2, i)) then
             write(*,'(i12,i13,i10,i12,f12.5)') held(1, i), held(2, i), &
                  & heads(e), determines(e), weight(e)
             exit
          end if
       end do
    end do

  end subroutine show_tuples

  subroutine edge_weights(steps, w)

    real(dp), intent(in) :: steps(:)
    real(dp), allocatable, intent(out) :: w(:)

    type(stored_directed_graph) :: edges
    type(stored_field) :: step_field, degrees, conditions
    type(scheme_weight) :: weights
    class(field), allocatable :: out

    edges = stored_directed_graph(num_instants, tails=tails, heads=heads)

    step_field = stored_field('dt', edges % vertex_set(), num_instants)
    degrees    = stored_field('source degree', edges % edge_set(), size(tails))
    conditions = stored_field('determines', edges % edge_set(), size(tails))
    call step_field % set_real_vector(steps)
    call degrees    % set_integer_vector(source_degree)
    call conditions % set_integer_vector(determines)

    weights = scheme_weight(bdf_family(order))
    call weights % apply(edges, [step_field, degrees, conditions], out)
    call out % real_vector(w)

  end subroutine edge_weights

end program coupling_relation
