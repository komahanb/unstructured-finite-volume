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
  use graph_carrier , only : counted_set, subset_set
  use graph_relation, only : stored_relation
  use graph_relation_algebra, only : compose_binary
  use graph_binary_relation , only : csr_relation, transposed_view, &
       &                             transpose_of, inclusion_of
  use graph_profile  , only : directed_adjacency_view
  use graph_algorithms, only : topological_order
  use fractal_graph        , only : graph, known_branch, null_branch
  use graph_relational_view, only : relational_binding, &
       & num_member_sets, member_set_at, num_relations, relation_at, &
       & holds_set

  implicit none


  type(counted_set)              :: v, t
  type(subset_set)               :: q_dom, y_dom
  type(stored_relation)          :: dep
  type(csr_relation), target     :: inc_y, inc_q, jq
  type(transposed_view)          :: inc_q_t, jq_t
  type(csr_relation)             :: coupling
  type(graph)             , target :: g
  type(graph)             , target :: scell(1), selem(1)
  type(graph)             , target :: rcell(1), relem(1)
  type(relational_binding)         :: bnd
  integer                          :: kcell
  type(directed_adjacency_view)  :: view
  integer, allocatable           :: order(:)
  integer                        :: table(2, 9)
  character(len=32)              :: which

  if (command_argument_count() .lt. 1) then
     error stop 'usage: refusal <case>'
  end if
  call get_command_argument(1, which)

  select case (trim(which))

  case ('cyclic-order')
     v = counted_set('variables', 3)
     t = counted_set('targets'  , 3)
     q_dom = subset_set('state'   , v, [VAR_U, VAR_V])
     y_dom = subset_set('residual', t, [TGT_R1, TGT_R2])

     table(:, 1) = [TGT_R1, VAR_P]
     table(:, 2) = [TGT_R1, VAR_U]
     table(:, 3) = [TGT_R1, VAR_V]
     table(:, 4) = [TGT_R2, VAR_P]
     table(:, 5) = [TGT_R2, VAR_U]
     table(:, 6) = [TGT_R2, VAR_V]
     table(:, 7) = [TGT_F , VAR_P]
     table(:, 8) = [TGT_F , VAR_U]
     table(:, 9) = [TGT_F , VAR_V]
     dep = stored_relation('dependency', [t, v], table)

     inc_y    = inclusion_of(y_dom)
     inc_q    = inclusion_of(q_dom)
     inc_q_t  = transpose_of(inc_q)
     jq       = compose_binary(compose_binary(inc_y, dep), inc_q_t)
     jq_t     = transpose_of(jq)
     coupling = compose_binary(jq_t, jq)

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
     view = directed_adjacency_view(g, bnd, coupling)

     call topological_order(view, order)
     write(*,*) 'an implicit system was given an execution order', order

  case default
     error stop 'usage: refusal <case>'

  end select

end program adjoint_level_4_refusal
