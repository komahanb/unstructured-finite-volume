!=====================================================================!
! ADJOINT TOWER . LEVEL 4 . THE REFUSAL
!
! One case, and it is the rung's whole point:
!
!   cyclic-order .... the implicitly coupled state system has no
!                     topological order, and the walk REFUSES rather
!                     than inventing one
!
! The coupling handed to the walk is derived exactly as the test
! derives it - C_Q = J_Q o J_Q^T from the one dependency source -
! so what dies here is the real structure of the specimen, not a
! contrived cycle built to provoke an error message.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program adjoint_level_4_refusal

  use adjoint_assert, only : VAR_P, VAR_U, VAR_V
  use adjoint_assert, only : TGT_R1, TGT_R2, TGT_F
  use graph_fractal        , only : set_graph => graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map
  use map_label      , only : label_map
  use map_inclusion  , only : inclusion_map, declared_subobject
  use relation_finitary, only : stored_relation
  use relation_algebra, only : compose_binary
  use relation_binary , only : csr_relation, transposed_view, &
       &                             transpose_of, inclusion_of
  use relation_algorithms, only : topological_order
  use graph_fractal        , only : graph, known_branch, null_branch
  use view_relational, only : relational_binding, &
       & num_member_sets, member_set_at, num_relations, relation_at, &
       & holds_set

  implicit none


  type(set_graph)              :: v, t
  type(set_graph)               :: q_dom, y_dom
  type(stored_relation)          :: dep
  type(csr_relation), target     :: inc_y, inc_q, jq
  type(transposed_view)          :: inc_q_t, jq_t
  type(csr_relation)             :: coupling
  type(graph)             , target :: g
  type(graph)             , target :: scell(1), selem(1)
  type(graph)             , target :: rcell(1), relem(1)
  type(relational_binding)         :: bnd
  integer                          :: kcell
  integer, allocatable           :: order(:)
  integer                        :: table(2, 9)
  character(len=32)              :: which
  type(set_map)     :: sets
  type(label_map)     :: labels
  type(inclusion_map)     :: inclusions

  if (command_argument_count() .lt. 1) then
     error stop 'usage: refusal <case>'
  end if
  call get_command_argument(1, which)

  select case (trim(which))

  case ('cyclic-order')
     call v % declare()
     call sets % bind(v, counted_set_representation(3))
     call t % declare()
     call sets % bind(t, counted_set_representation(3))
     call q_dom % declare()
     call sets       % bind(q_dom, listed_set_representation([VAR_U, VAR_V]))
     call inclusions % include_in(q_dom, v)
     call y_dom % declare()
     call sets       % bind(y_dom, listed_set_representation([TGT_R1, TGT_R2]))
     call inclusions % include_in(y_dom, t)

     table(:, 1) = [TGT_R1, VAR_P]
     table(:, 2) = [TGT_R1, VAR_U]
     table(:, 3) = [TGT_R1, VAR_V]
     table(:, 4) = [TGT_R2, VAR_P]
     table(:, 5) = [TGT_R2, VAR_U]
     table(:, 6) = [TGT_R2, VAR_V]
     table(:, 7) = [TGT_F , VAR_P]
     table(:, 8) = [TGT_F , VAR_U]
     table(:, 9) = [TGT_F , VAR_V]
     dep = stored_relation('dependency', [t, v], table, sets)

     inc_y    = inclusion_of(y_dom, t, sets, labels)
     inc_q    = inclusion_of(q_dom, v, sets, labels)
     inc_q_t  = transpose_of(inc_q)
     jq       = compose_binary(compose_binary(inc_y, dep, sets), inc_q_t, sets)
     jq_t     = transpose_of(jq)
     coupling = compose_binary(jq_t, jq, sets)

     ! 'state coupling': (S, P) as one sequence on each branch.
     call g % declare()
     do kcell = 1, 1
        call scell(kcell) % declare()
        call selem(kcell) % declare()
     end do
     do kcell = 1, 1
        call rcell(kcell) % declare()
        call relem(kcell) % declare()
     end do

     call bnd % bind_set(selem(1), q_dom)
     call bnd % bind_relation(relem(1), coupling)

     do kcell = 1, 1
        scell(kcell) % branch(1) = known_branch(selem(kcell))
        if (kcell .lt. 1) scell(kcell) % branch(2) = &
             & known_branch(scell(kcell + 1))
     end do
     do kcell = 1, 1
        rcell(kcell) % branch(1) = known_branch(relem(kcell))
        if (kcell .lt. 1) rcell(kcell) % branch(2) = &
             & known_branch(rcell(kcell + 1))
     end do

     g % branch(1) = known_branch(scell(1))
     g % branch(2) = known_branch(rcell(1))

     call topological_order(coupling, sets, order)
     write(*,*) 'an implicit system was given an execution order', order

  case default
     error stop 'usage: refusal <case>'

  end select

end program adjoint_level_4_refusal
