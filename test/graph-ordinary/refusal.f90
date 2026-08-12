!=====================================================================!
! The profile refusals, each EXPECTED TO DIE for its stated reason:
!
!      tailless     an edge with no tuple in T
!      twotailed    an edge with two tuples in T
!      twoheaded    an edge with two tuples in H
!      ternary      a relation of the wrong arity offered as T
!      mismatched   H over domains T does not share
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program ordinary_refusal

  use graph_carrier        , only : counted_set
  use graph_relation       , only : stored_relation, slot
  use graph_binary_relation, only : csr_relation
  use graph_structure      , only : relational_graph, held_set, &
       &                            held_relation
  use graph_profile        , only : ordinary_graph_view

  implicit none

  type(counted_set)              :: verts, edges, other, roles
  type(csr_relation)             :: t, h
  type(stored_relation)          :: fat
  type(relational_graph), target :: g
  type(ordinary_graph_view)      :: view
  character(len=32)              :: which

  which = ''
  call get_command_argument(1, which)

  verts = counted_set('vertices', 3)
  edges = counted_set('edges'   , 2)

  select case (trim(which))

  case ('tailless')
     t = csr_relation('tail', edges, verts, reshape([1, 1], [2, 1]))
     h = csr_relation('head', edges, verts, reshape([1, 2], [2, 1]))
     g = relational_graph('bad', [held_set(verts), held_set(edges)], &
          & [held_relation(t), held_relation(h)])
     view = ordinary_graph_view(g, tail_at=1, head_at=2)

  case ('twotailed')
     t = csr_relation('tail', edges, verts, &
          & reshape([1,1,  1,2,  2,2], [2, 3]))
     h = csr_relation('head', edges, verts, reshape([1, 2], [2, 1]))
     g = relational_graph('bad', [held_set(verts), held_set(edges)], &
          & [held_relation(t), held_relation(h)])
     view = ordinary_graph_view(g, tail_at=1, head_at=2)

  case ('twoheaded')
     t = csr_relation('tail', edges, verts, &
          & reshape([1,1,  2,2], [2, 2]))
     h = csr_relation('head', edges, verts, &
          & reshape([1,2,  1,3,  2,3], [2, 3]))
     g = relational_graph('bad', [held_set(verts), held_set(edges)], &
          & [held_relation(t), held_relation(h)])
     view = ordinary_graph_view(g, tail_at=1, head_at=2)

  case ('ternary')
     roles = counted_set('roles', 2)
     fat = stored_relation('endpoint', &
          & [slot(edges), slot(verts), slot(roles)], &
          & reshape([1,1,1,  2,2,1], [3, 2]))
     h = csr_relation('head', edges, verts, reshape([1, 2], [2, 1]))
     g = relational_graph('bad', &
          & [held_set(verts), held_set(edges), held_set(roles)], &
          & [held_relation(fat), held_relation(h)])
     view = ordinary_graph_view(g, tail_at=1, head_at=2)

  case ('mismatched')
     other = counted_set('elsewhere', 3)
     t = csr_relation('tail', edges, verts, &
          & reshape([1,1,  2,2], [2, 2]))
     h = csr_relation('head', edges, other, reshape([1, 2], [2, 1]))
     g = relational_graph('bad', &
          & [held_set(verts), held_set(edges), held_set(other)], &
          & [held_relation(t), held_relation(h)])
     view = ordinary_graph_view(g, tail_at=1, head_at=2)

  case default
     write(*,'(1x,a)') &
          & "usage: refusal tailless|twotailed|twoheaded|ternary|mismatched"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program ordinary_refusal
