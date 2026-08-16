!=====================================================================!
! The interpretation and algorithm refusals, each EXPECTED TO DIE
! for its stated reason:
!
!      notowned     the selector names a relation the graph lacks
!      notbinary    a ternary relation offered as adjacency
!      notsquare    a binary relation over two different domains
!      cycle        a lawful cyclic view whose topological order
!                   is rightly refused - the VIEW is valid; only
!                   the ordering is undefined
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program algorithms_refusal

  use graph_carrier        , only : counted_set
  use graph_relation       , only : stored_relation, slot
  use graph_binary_relation, only : csr_relation
  use graph_profile        , only : directed_adjacency_view
  use graph_algorithms     , only : topological_order
  use fractal_graph        , only : graph, known_branch, null_branch
  use graph_relational_view, only : relational_binding, &
       & num_member_sets, member_set_at, num_relations, relation_at, &
       & holds_set

  implicit none


  type(counted_set)              :: a, b, c
  type(csr_relation)             :: adj, stray, lopsided, ring
  type(stored_relation)          :: fat
  type(graph)             , target :: g
  ! Sized for the widest case; each case builds its own.
  type(graph)             , target :: scell(2), selem(2)
  type(graph)             , target :: rcell(1), relem(1)
  type(relational_binding)         :: bnd
  integer                          :: kcell
  type(directed_adjacency_view)  :: view
  integer, allocatable           :: order(:)
  character(len=32)              :: which

  which = ''
  call get_command_argument(1, which)

  a = counted_set('a-things', 3)
  b = counted_set('b-things', 2)

  select case (trim(which))

  case ('notowned')
     adj   = csr_relation('inside' , a, a, reshape([1, 2], [2, 1]))
     stray = csr_relation('outside', a, a, reshape([2, 3], [2, 1]))
     ! 'bad': (S, P) as one sequence on each branch.
     call g % declare()
     do kcell = 1, 1
        call scell(kcell) % declare()
        call selem(kcell) % declare()
     end do
     do kcell = 1, 1
        call rcell(kcell) % declare()
        call relem(kcell) % declare()
     end do

     call bnd % bind_set(selem(1), a)
     call bnd % bind_relation(relem(1), adj)

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
     view = directed_adjacency_view(g, bnd, stray)

  case ('notbinary')
     c   = counted_set('c-things', 2)
     fat = stored_relation('fat', [slot(a), slot(a), slot(c)], &
          & reshape([1, 2, 1], [3, 1]))
     ! 'bad': (S, P) as one sequence on each branch.
     call g % declare()
     do kcell = 1, 2
        call scell(kcell) % declare()
        call selem(kcell) % declare()
     end do
     do kcell = 1, 1
        call rcell(kcell) % declare()
        call relem(kcell) % declare()
     end do

     call bnd % bind_set(selem(1), a)
     call bnd % bind_set(selem(2), c)
     call bnd % bind_relation(relem(1), fat)

     do kcell = 1, 2
        scell(kcell) % branch(1) = known_branch(selem(kcell))
        if (kcell .lt. 2) scell(kcell) % branch(2) = &
             & known_branch(scell(kcell + 1))
     end do
     do kcell = 1, 1
        rcell(kcell) % branch(1) = known_branch(relem(kcell))
        if (kcell .lt. 1) rcell(kcell) % branch(2) = &
             & known_branch(rcell(kcell + 1))
     end do

     g % branch(1) = known_branch(scell(1))
     g % branch(2) = known_branch(rcell(1))
     view = directed_adjacency_view(g, bnd, fat)

  case ('notsquare')
     lopsided = csr_relation('lopsided', a, b, reshape([1, 1], [2, 1]))
     ! 'bad': (S, P) as one sequence on each branch.
     call g % declare()
     do kcell = 1, 2
        call scell(kcell) % declare()
        call selem(kcell) % declare()
     end do
     do kcell = 1, 1
        call rcell(kcell) % declare()
        call relem(kcell) % declare()
     end do

     call bnd % bind_set(selem(1), a)
     call bnd % bind_set(selem(2), b)
     call bnd % bind_relation(relem(1), lopsided)

     do kcell = 1, 2
        scell(kcell) % branch(1) = known_branch(selem(kcell))
        if (kcell .lt. 2) scell(kcell) % branch(2) = &
             & known_branch(scell(kcell + 1))
     end do
     do kcell = 1, 1
        rcell(kcell) % branch(1) = known_branch(relem(kcell))
        if (kcell .lt. 1) rcell(kcell) % branch(2) = &
             & known_branch(rcell(kcell + 1))
     end do

     g % branch(1) = known_branch(scell(1))
     g % branch(2) = known_branch(rcell(1))
     view = directed_adjacency_view(g, bnd, lopsided)

  case ('cycle')
     ring = csr_relation('ring', a, a, &
          & reshape([1,2,  2,3,  3,1], [2, 3]))
     ! 'cyclic': (S, P) as one sequence on each branch.
     call g % declare()
     do kcell = 1, 1
        call scell(kcell) % declare()
        call selem(kcell) % declare()
     end do
     do kcell = 1, 1
        call rcell(kcell) % declare()
        call relem(kcell) % declare()
     end do

     call bnd % bind_set(selem(1), a)
     call bnd % bind_relation(relem(1), ring)

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
     ! The view is lawful: a cycle is a directed graph.
     view = directed_adjacency_view(g, bnd, ring)
     ! Only the ordering is undefined.
     call topological_order(view, order)

  case default
     write(*,'(1x,a)') "usage: refusal notowned|notbinary|notsquare|cycle"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program algorithms_refusal
