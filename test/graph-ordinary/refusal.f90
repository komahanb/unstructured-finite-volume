!=====================================================================!
! The interpretation refusals, each EXPECTED TO DIE for its stated reason:
!
!      tailless     an edge with no tuple in T
!      twotailed    an edge with two tuples in T
!      twoheaded    an edge with two tuples in H
!      ternary      a relation of the wrong arity offered as T
!      mismatched   H over domains T does not share
!      selfsame     one domain playing both edges and vertices
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program ordinary_refusal

  use graph_set        , only : index_set
  use graph_relation       , only : stored_relation, declared_domain
  use graph_binary_relation, only : csr_relation
  use graph_structure      , only : related_graph, declared_set, &
       &                            declared_relation
  use graph_interpretation        , only : ordinary_graph_view

  implicit none

  type(index_set)              :: verts, edges, other, roles
  type(csr_relation)             :: t, h
  type(stored_relation)          :: fat
  type(related_graph), target :: g
  type(ordinary_graph_view)      :: view
  character(len=32)              :: which

  which = ''
  call get_command_argument(1, which)

  verts = index_set('vertices', 3)
  edges = index_set('edges'   , 2)

  select case (trim(which))

  case ('tailless')
     t = csr_relation('tail', edges, verts, reshape([1, 1], [2, 1]))
     h = csr_relation('head', edges, verts, reshape([1, 2], [2, 1]))
     g = related_graph('bad', [declared_set(verts), declared_set(edges)], &
          & [declared_relation(t), declared_relation(h)])
     view = ordinary_graph_view(g, tail_at=1, head_at=2)

  case ('twotailed')
     t = csr_relation('tail', edges, verts, &
          & reshape([1,1,  1,2,  2,2], [2, 3]))
     h = csr_relation('head', edges, verts, reshape([1, 2], [2, 1]))
     g = related_graph('bad', [declared_set(verts), declared_set(edges)], &
          & [declared_relation(t), declared_relation(h)])
     view = ordinary_graph_view(g, tail_at=1, head_at=2)

  case ('twoheaded')
     t = csr_relation('tail', edges, verts, &
          & reshape([1,1,  2,2], [2, 2]))
     h = csr_relation('head', edges, verts, &
          & reshape([1,2,  1,3,  2,3], [2, 3]))
     g = related_graph('bad', [declared_set(verts), declared_set(edges)], &
          & [declared_relation(t), declared_relation(h)])
     view = ordinary_graph_view(g, tail_at=1, head_at=2)

  case ('ternary')
     roles = index_set('roles', 2)
     fat = stored_relation('endpoint', &
          & [declared_domain(edges), declared_domain(verts), declared_domain(roles)], &
          & reshape([1,1,1,  2,2,1], [3, 2]))
     h = csr_relation('head', edges, verts, reshape([1, 2], [2, 1]))
     g = related_graph('bad', &
          & [declared_set(verts), declared_set(edges), declared_set(roles)], &
          & [declared_relation(fat), declared_relation(h)])
     view = ordinary_graph_view(g, tail_at=1, head_at=2)

  case ('mismatched')
     other = index_set('elsewhere', 3)
     t = csr_relation('tail', edges, verts, &
          & reshape([1,1,  2,2], [2, 2]))
     h = csr_relation('head', edges, other, reshape([1, 2], [2, 1]))
     g = related_graph('bad', &
          & [declared_set(verts), declared_set(edges), declared_set(other)], &
          & [declared_relation(t), declared_relation(h)])
     view = ordinary_graph_view(g, tail_at=1, head_at=2)

  case ('selfsame')
     t = csr_relation('tail', edges, edges, reshape([1, 2], [2, 1]))
     h = csr_relation('head', edges, edges, reshape([2, 1], [2, 1]))
     g = related_graph('bad', [declared_set(edges)], &
          & [declared_relation(t), declared_relation(h)])
     view = ordinary_graph_view(g, tail_at=1, head_at=2)

  case default
     write(*,'(1x,a)') &
          & "usage: refusal tailless|twotailed|twoheaded|ternary|mismatched|selfsame"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program ordinary_refusal
