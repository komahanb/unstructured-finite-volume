!=====================================================================!
! The structure refusals, each EXPECTED TO DIE for its stated
! reason:
!
!      unheld      a relation over a domain the graph does not hold
!      twice        a second signing of a declared graph
!      undeclared   a seat holding a domain that never signed
!      dupset       one domain seated twice
!      duprel       one relation seated twice
!      view         a borrowing view offered for ownership
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program structure_refusal

  use graph_set        , only : index_set, unrelated_graph, declared_set
  use graph_binary_relation, only : csr_relation, transposed_view, &
       &                            transpose_of
  use graph_structure      , only : related_graph, declared_relation, &
       &                            declare_relations

  implicit none

  type(index_set)          :: ops, unheld, raw
  type(csr_relation), target :: dep
  type(transposed_view)      :: flipped
  type(related_graph)     :: g
  type(unrelated_graph)   :: ug
  type(declared_relation)        :: none(0)
  character(len=32)          :: which

  which = ''
  call get_command_argument(1, which)

  ops     = index_set('operations', 3)
  unheld = index_set('a domain this graph does not hold', 3)

  select case (trim(which))

  case ('unheld-domain')
     dep = csr_relation('feeds', ops, unheld, reshape([1, 1], [2, 1]))
     g = related_graph('bad', [declared_set(ops)], [declared_relation(dep)])

  case ('twice')
     ug = unrelated_graph('good', [declared_set(ops)])
     call ug % declare('good-again')

  case ('empty-relation-family')
     ! Declared domains, and no relation declared over them. That is
     ! an unrelated_graph, and related_graph is not that type.
     g = related_graph('bad', [declared_set(ops)], none)

  case ('map-empty-family')
     ! The map delegates, so the constructor's |R| > 0 law fires.
     ug = unrelated_graph('over', [declared_set(ops)])
     g  = declare_relations(ug, none, 'bad')

  case ('map-unheld-domain')
     dep = csr_relation('feeds', ops, unheld, reshape([1, 1], [2, 1]))
     ug  = unrelated_graph('over', [declared_set(ops)])
     g   = declare_relations(ug, [declared_relation(dep)], 'bad')

  case ('map-view')
     dep     = csr_relation('feeds', ops, ops, reshape([1, 2], [2, 1]))
     flipped = transpose_of(dep)
     ug      = unrelated_graph('over', [declared_set(ops)])
     g       = declare_relations(ug, [declared_relation(flipped)], 'bad')

  case ('undeclared')
     g = related_graph('bad', [declared_set(ops), declared_set(raw)], none)

  case ('dupset')
     g = related_graph('bad', [declared_set(ops), declared_set(ops)], none)

  case ('duprel')
     dep = csr_relation('feeds', ops, ops, reshape([1, 2], [2, 1]))
     g = related_graph('bad', [declared_set(ops)], &
          & [declared_relation(dep), declared_relation(dep)])

  case ('view')
     dep = csr_relation('feeds', ops, ops, reshape([1, 2], [2, 1]))
     flipped = transpose_of(dep)
     g = related_graph('bad', [declared_set(ops)], &
          & [declared_relation(flipped)])

  case default
     write(*,'(1x,a)') &
          & "usage: refusal unheld-domain|empty-relation-family|map-empty-family|" // &
          & "map-unheld-domain|map-view|twice|undeclared|dupset|duprel|view"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program structure_refusal
