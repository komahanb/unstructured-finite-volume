!=====================================================================!
! The profile refusals, each EXPECTED TO DIE for its stated reason:
!
!      tailless     an edge with no tuple in T
!      twotailed    an edge with two tuples in T
!      twoheaded    an edge with two tuples in H
!      ternary      a relation of the wrong arity offered as T
!      mismatched   H over domains T does not share
!      selfsame     one domain playing both edges and vertices
!
! Each case binds the member sets and the two relations it means; the
! representation is then wired once, below, and the profile reads it.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program ordinary_refusal

  use graph_set_representation, only : counted_set_representation
  use graph_set_map           , only : set_map
  use graph_relation       , only : stored_relation
  use graph_binary_relation, only : csr_relation
  use graph_profile        , only : directed_incidence_view
  use fractal_graph        , only : graph, known_branch
  use graph_relational_view, only : relational_binding

  implicit none

  type(graph)         :: verts, edges, other, roles
  type(set_map)         :: sets
  type(csr_relation)        :: t, h
  type(stored_relation)     :: fat
  type(directed_incidence_view) :: view
  character(len=32)         :: which

  type(graph)             , target :: g
  type(graph)             , target :: scell(3), selem(3), rcell(2), relem(2)
  type(relational_binding)         :: bnd
  integer                          :: ns, k

  which = ''
  call get_command_argument(1, which)

  call g % declare()
  do k = 1, 3
     call scell(k) % declare(); call selem(k) % declare()
  end do
  do k = 1, 2
     call rcell(k) % declare(); call relem(k) % declare()
  end do

  call verts % declare()
  call sets % bind(verts, counted_set_representation(3))
  call edges % declare()
  call sets % bind(edges, counted_set_representation(2))

  select case (trim(which))

  case ('tailless')
     t = csr_relation('tail', edges, verts, reshape([1, 1], [2, 1]), sets)
     h = csr_relation('head', edges, verts, reshape([1, 2], [2, 1]), sets)
     ns = 2
     call bnd % bind_set(selem(1), verts)
     call bnd % bind_set(selem(2), edges)
     call bnd % bind_relation(relem(1), t)
     call bnd % bind_relation(relem(2), h)

  case ('twotailed')
     t = csr_relation('tail', edges, verts, &
          & reshape([1,1,  1,2,  2,2], [2, 3]), sets)
     h = csr_relation('head', edges, verts, reshape([1, 2], [2, 1]), sets)
     ns = 2
     call bnd % bind_set(selem(1), verts)
     call bnd % bind_set(selem(2), edges)
     call bnd % bind_relation(relem(1), t)
     call bnd % bind_relation(relem(2), h)

  case ('twoheaded')
     t = csr_relation('tail', edges, verts, &
          & reshape([1,1,  2,2], [2, 2]), sets)
     h = csr_relation('head', edges, verts, &
          & reshape([1,2,  1,3,  2,3], [2, 3]), sets)
     ns = 2
     call bnd % bind_set(selem(1), verts)
     call bnd % bind_set(selem(2), edges)
     call bnd % bind_relation(relem(1), t)
     call bnd % bind_relation(relem(2), h)

  case ('ternary')
     call roles % declare()
     call sets % bind(roles, counted_set_representation(2))
     fat = stored_relation('endpoint', &
          & [edges, verts, roles], &
          & reshape([1,1,1,  2,2,1], [3, 2]), sets)
     h = csr_relation('head', edges, verts, reshape([1, 2], [2, 1]), sets)
     ns = 3
     call bnd % bind_set(selem(1), verts)
     call bnd % bind_set(selem(2), edges)
     call bnd % bind_set(selem(3), roles)
     call bnd % bind_relation(relem(1), fat)
     call bnd % bind_relation(relem(2), h)

  case ('mismatched')
     call other % declare()
     call sets % bind(other, counted_set_representation(3))
     t = csr_relation('tail', edges, verts, &
          & reshape([1,1,  2,2], [2, 2]), sets)
     h = csr_relation('head', edges, other, reshape([1, 2], [2, 1]), sets)
     ns = 3
     call bnd % bind_set(selem(1), verts)
     call bnd % bind_set(selem(2), edges)
     call bnd % bind_set(selem(3), other)
     call bnd % bind_relation(relem(1), t)
     call bnd % bind_relation(relem(2), h)

  case ('selfsame')
     t = csr_relation('tail', edges, edges, reshape([1, 2], [2, 1]), sets)
     h = csr_relation('head', edges, edges, reshape([2, 1], [2, 1]), sets)
     ns = 1
     call bnd % bind_set(selem(1), edges)
     call bnd % bind_relation(relem(1), t)
     call bnd % bind_relation(relem(2), h)

  case default
     write(*,'(1x,a)') &
          & "usage: refusal tailless|twotailed|twoheaded|ternary|mismatched|selfsame"
     error stop 'no case chosen'

  end select

  do k = 1, ns
     scell(k) % branch(1) = known_branch(selem(k))
     if (k .lt. ns) scell(k) % branch(2) = known_branch(scell(k + 1))
  end do

  do k = 1, 2
     rcell(k) % branch(1) = known_branch(relem(k))
     if (k .lt. 2) rcell(k) % branch(2) = known_branch(rcell(k + 1))
  end do

  g % branch(1) = known_branch(scell(1))
  g % branch(2) = known_branch(rcell(1))

  view = directed_incidence_view(g, bnd, sets, tail_at=1, head_at=2)

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program ordinary_refusal
