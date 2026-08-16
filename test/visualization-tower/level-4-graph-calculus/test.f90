!=====================================================================!
! VISUALIZATION TOWER . LEVEL 4 . STRUCTURAL INTERPRETATION
!
! The level answers one question: CAN THE EXISTING STRUCTURE BE
! RENDERED MEANINGFULLY.
!
! Five representations, all of them generated:
!
!      A   the forward chain     X0 --D1--> X1 --D2--> X2 --D3--> X3
!      B   the reverse chain     X3 --D3^T--> ... --D1^T--> X0
!      C   D1, D2, D3            one sparsity picture each
!      D   D3:1                  the whole chain's skeleton
!      E   D3:1^T                and the same, reversed
!
! and one architectural question asked of the nucleus itself.
!
!                       THE PICTURE IS DOWNSTREAM
!
! Every glyph below comes from relation % has, every axis from
! carrier % member, every name from the object that carries it. The
! expected pictures written into this test are the ORACLE - what a
! person worked out on paper from the twelve occurrences. They are
! not an input to the drawing, and the renderer has never seen them.
!
! The forbidden direction would be to write the picture first and
! then invent the relation that reproduces it. That is why the
! oracle lives here, in the test, and not one line of it is
! reachable from common/structural_renderer_fixture.f90.
!
!                     RENDERED FROM THE TRANSPOSE, NOT THE STRING
!
! Representation E is drawn twice: once from D3:1^T, and once from
! D1^T o D2^T o D3^T, composed independently. The two pictures must
! agree grid for grid while differing in their name line - which is
! the composition/transpose law of Level 2, seen through a
! representation instead of through same_extension.
!
!                        THE HOSTILE CARRIER
!
! One carrier here declares its members { 30, 10, 20 } in that order,
! and the renderer must print them in that order. Nothing else in the
! specimen could catch a renderer that sorted, or that assumed a
! member is its own position, because every other carrier counts from
! one and would look correct either way.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program visualization_level_4

  use visualization_assert , only : report, verdict
  use visualization_assert , only : ND31
  use visualization_assert , only : X0_A, X0_B, X0_C, X0_D
  use visualization_assert , only : X1_P, X1_Q, X1_R
  use visualization_assert , only : X2_U, X2_V, X2_W
  use visualization_assert , only : X3_M, X3_N
  use graph_carrier        , only : counted_set, subset_set, member_set
  use graph_relation       , only : relation
  use graph_binary_relation, only : csr_relation, binary_relation
  use fractal_graph        , only : graph, known_branch, null_branch
  use graph_relational_view, only : relational_binding, &
       & num_member_sets, member_set_at, num_relations, relation_at, &
       & holds_set
  use visualization_carriers_fixture , only : structural_carriers, label_for
  use visualization_relations_fixture, only : occurrences_of_a1
  use visualization_relations_fixture, only : occurrences_of_a2
  use visualization_relations_fixture, only : occurrences_of_a3
  use visualization_algebra_fixture  , only : derive_dependency
  use visualization_algebra_fixture  , only : derive_composition
  use visualization_algebra_fixture  , only : materialized_transpose
  use visualization_algebra_fixture  , only : same_extension
  use structural_renderer_fixture    , only : picture, stage, glyph_at
  use structural_renderer_fixture    , only : sparsity_picture
  use structural_renderer_fixture    , only : chain_picture
  use structural_renderer_fixture    , only : dependency_listing

  implicit none

  type(counted_set)          :: x0, x1, x2, x3, e1, e2, e3
  type(csr_relation)         :: t1, h1, t2, h2, t3, h3
  type(csr_relation), target :: d1, d2, d3, d21, d31
  type(csr_relation), target :: d1t, d2t, d3t
  type(csr_relation)         :: middle, drev, d31t
  type(graph)             , target :: g
  type(graph)             , target :: scell(7), selem(7)
  type(graph)             , target :: rcell(6), relem(6)
  type(relational_binding)         :: bnd
  integer                          :: kcell
  integer                    :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "visualization tower . level 4 . structural interpretation"
  write(*,'(1x,a)') "============================================="

  call structural_carriers(x0, x1, x2, x3, e1, e2, e3)
  call occurrences_of_a1(e1, x0, x1, t1, h1)
  call occurrences_of_a2(e2, x1, x2, t2, h2)
  call occurrences_of_a3(e3, x2, x3, t3, h3)

  ! 'the operator chain A3 o A2 o A1': (S, P) as one sequence on each branch.
  call g % declare()
  do kcell = 1, 7
     call scell(kcell) % declare()
     call selem(kcell) % declare()
  end do
  do kcell = 1, 6
     call rcell(kcell) % declare()
     call relem(kcell) % declare()
  end do

  call bnd % bind_set(selem(1), x0)
  call bnd % bind_set(selem(2), x1)
  call bnd % bind_set(selem(3), x2)
  call bnd % bind_set(selem(4), x3)
  call bnd % bind_set(selem(5), e1)
  call bnd % bind_set(selem(6), e2)
  call bnd % bind_set(selem(7), e3)
  call bnd % bind_relation(relem(1), t1)
  call bnd % bind_relation(relem(2), h1)
  call bnd % bind_relation(relem(3), t2)
  call bnd % bind_relation(relem(4), h2)
  call bnd % bind_relation(relem(5), t3)
  call bnd % bind_relation(relem(6), h3)

  do kcell = 1, 7
     scell(kcell) % branch(1) = known_branch(selem(kcell))
     if (kcell .lt. 7) scell(kcell) % branch(2) = &
          & known_branch(scell(kcell + 1))
  end do
  do kcell = 1, 6
     rcell(kcell) % branch(1) = known_branch(relem(kcell))
     if (kcell .lt. 6) rcell(kcell) % branch(2) = &
          & known_branch(rcell(kcell + 1))
  end do

  g % branch(1) = known_branch(scell(1))
  g % branch(2) = known_branch(rcell(1))

  d1 = derive_dependency('D1', t1, h1)
  d2 = derive_dependency('D2', t2, h2)
  d3 = derive_dependency('D3', t3, h3)

  d21 = derive_composition('D2:1', d1,  d2)
  d31 = derive_composition('D3:1', d21, d3)

  d1t = materialized_transpose('D1^T', d1)
  d2t = materialized_transpose('D2^T', d2)
  d3t = materialized_transpose('D3^T', d3)

  middle = derive_composition('D2^T o D3^T', d3t, d2t)
  drev   = derive_composition('D1^T o D2^T o D3^T', middle, d1t)
  d31t   = materialized_transpose('D3:1^T', d31)

  call say_the_structure()

  call check_representation_a_forward_chain(nfail)
  call check_representation_b_reverse_chain(nfail)
  call check_representation_c_individual_sparsity(nfail)
  call check_representation_d_composed_sparsity(nfail)
  call check_representation_e_reverse_sparsity(nfail)
  call check_the_renderer_invents_nothing(nfail)
  call check_declaration_order_rules(nfail)
  call check_the_ordinary_graph_question(nfail)

  call verdict(nfail, "level 4")

contains

  !===================================================================!
  ! The generated summary, printed for a person to look at. Every
  ! line of it was derived a moment ago from twelve occurrences.
  !===================================================================!

  subroutine say_the_structure()

    type(picture) :: pic

    write(*,'(1x,a)') "---------------------------------------------"

    pic = chain_picture('FORWARD CHAIN', [stage(d1), stage(d2), stage(d3)])
    call pic % emit()
    write(*,*)

    pic = sparsity_picture(d1);   call pic % emit(); write(*,*)
    pic = sparsity_picture(d2);   call pic % emit(); write(*,*)
    pic = sparsity_picture(d3);   call pic % emit(); write(*,*)
    pic = sparsity_picture(d21);  call pic % emit(); write(*,*)
    pic = sparsity_picture(d31);  call pic % emit(); write(*,*)

    pic = chain_picture('STRUCTURAL REVERSE', &
         &              [stage(d3t), stage(d2t), stage(d1t)])
    call pic % emit()
    write(*,*)

    pic = sparsity_picture(d31t); call pic % emit(); write(*,*)
    pic = sparsity_picture(drev); call pic % emit(); write(*,*)

    pic = dependency_listing(d31); call pic % emit()

    write(*,'(1x,a)') "---------------------------------------------"

  end subroutine say_the_structure

  !===================================================================!
  ! REPRESENTATION A. The forward chain, with every carrier name and
  ! every arrow read off the relations.
  !===================================================================!

  subroutine check_representation_a_forward_chain(nfail)

    integer, intent(inout) :: nfail

    type(picture) :: pic

    pic = chain_picture('FORWARD CHAIN', [stage(d1), stage(d2), stage(d3)])

    call report(pic % at(2) .eq. 'X0 --D1--> X1 --D2--> X2 --D3--> X3', &
         & "FORWARD : X0 --D1--> X1 --D2--> X2 --D3--> X3", nfail)

    pic = dependency_listing(d1)
    call report(pic % at(2) .eq. 'a -> p' .and. pic % at(3) .eq. 'b -> p q' .and. &
         &      pic % at(4) .eq. 'c -> q' .and. pic % at(5) .eq. 'd -> r', &
         & "and the same structure as a deterministic listing: " // &
         & "a->p, b->p q, c->q, d->r", nfail)

    pic = dependency_listing(d31)
    call report(pic % at(2) .eq. 'a -> m' .and. pic % at(3) .eq. 'b -> m n' .and. &
         &      pic % at(4) .eq. 'c -> m n' .and. pic % at(5) .eq. 'd -> n', &
         & "the chain's whole reach, listed: a->m, b->m n, c->m n, d->n", &
         & nfail)

  end subroutine check_representation_a_forward_chain

  !===================================================================!
  ! REPRESENTATION B. The reverse chain - the SAME routine, handed
  ! the transposed relations. No string was reversed to make it.
  !===================================================================!

  subroutine check_representation_b_reverse_chain(nfail)

    integer, intent(inout) :: nfail

    type(picture) :: pic

    pic = chain_picture('STRUCTURAL REVERSE', &
         &              [stage(d3t), stage(d2t), stage(d1t)])

    call report(pic % at(2) .eq. 'X3 --D3^T--> X2 --D2^T--> X1 --D1^T--> X0', &
         & "REVERSE : X3 --D3^T--> X2 --D2^T--> X1 --D1^T--> X0 - " // &
         & "drawn from the transposed relations, not from the " // &
         & "forward line read backwards", nfail)

    pic = dependency_listing(d1t)
    call report(pic % at(2) .eq. 'p -> a b' .and. pic % at(3) .eq. 'q -> b c' .and. &
         &      pic % at(4) .eq. 'r -> d', &
         & "and D1^T listed: p is reached from a and b, q from b and " // &
         & "c, r from d", nfail)

  end subroutine check_representation_b_reverse_chain

  !===================================================================!
  ! REPRESENTATION C. One picture per operator, rows the codomain,
  ! columns the domain, glyph the relation's own membership.
  !===================================================================!

  subroutine check_representation_c_individual_sparsity(nfail)

    integer, intent(inout) :: nfail

    type(picture) :: pic

    pic = sparsity_picture(d1)
    call report(pic % rows() .eq. 5 .and. &
         &      pic % at(1) .eq. 'D1' .and. &
         &      pic % at(2) .eq. '        a b c d' .and. &
         &      pic % at(3) .eq. 'p       # # . .' .and. &
         &      pic % at(4) .eq. 'q       . # # .' .and. &
         &      pic % at(5) .eq. 'r       . . . #', &
         & "D1 renders exactly the sparsity of { a->p, b->p, b->q, " // &
         & "c->q, d->r }", nfail)

    pic = sparsity_picture(d2)
    call report(pic % at(1) .eq. 'D2' .and. &
         &      pic % at(2) .eq. '        p q r' .and. &
         &      pic % at(3) .eq. 'u       # # .' .and. &
         &      pic % at(4) .eq. 'v       . # .' .and. &
         &      pic % at(5) .eq. 'w       . . #', &
         & "D2 renders exactly the sparsity of { p->u, q->u, q->v, " // &
         & "r->w }", nfail)

    pic = sparsity_picture(d3)
    call report(pic % at(1) .eq. 'D3' .and. &
         &      pic % at(2) .eq. '        u v w' .and. &
         &      pic % at(3) .eq. 'm       # . .' .and. &
         &      pic % at(4) .eq. 'n       . # #', &
         & "D3 renders exactly the sparsity of { u->m, v->n, w->n }", &
         & nfail)

    ! The pictures are rectangular because the relations are. D1 is
    ! three rows by four columns, and nothing about the nucleus
    ! wanted it square.
    block
      type(picture) :: one, three
      one   = sparsity_picture(d1)
      three = sparsity_picture(d3)
      call report(one % rows() .eq. 2 + 3 .and. three % rows() .eq. 2 + 2, &
           & "and each picture is |codomain| rows by |domain| columns - " // &
           & "RECTANGULAR, as the typed dependency is", nfail)
    end block

  end subroutine check_representation_c_individual_sparsity

  !===================================================================!
  ! REPRESENTATION D. The composed skeleton - one of this tower's
  ! two central pictures.
  !===================================================================!

  subroutine check_representation_d_composed_sparsity(nfail)

    integer, intent(inout) :: nfail

    type(picture) :: pic

    pic = sparsity_picture(d21)
    call report(pic % at(2) .eq. '        a b c d' .and. &
         &      pic % at(3) .eq. 'u       # # # .' .and. &
         &      pic % at(4) .eq. 'v       . # # .' .and. &
         &      pic % at(5) .eq. 'w       . . . #', &
         & "D2:1 renders the two-operator skeleton - and b->u fills " // &
         & "ONE cell though two walks reach it", nfail)

    pic = sparsity_picture(d31)
    call report(pic % rows() .eq. 4 .and. &
         &      pic % at(1) .eq. 'D3:1' .and. &
         &      pic % at(2) .eq. '        a b c d' .and. &
         &      pic % at(3) .eq. 'm       # # # .' .and. &
         &      pic % at(4) .eq. 'n       . # # #', &
         & "D3:1 renders THE STRUCTURAL SKELETON OF THE COMPLETE " // &
         & "OPERATOR CHAIN", nfail)

    call report(filled_cells(d31) .eq. ND31, &
         & "six filled cells, six tuples - the picture counts what " // &
         & "the relation holds", nfail)

  end subroutine check_representation_d_composed_sparsity

  !===================================================================!
  ! REPRESENTATION E. The reverse skeleton, drawn twice from two
  ! independently composed relations, and required to agree.
  !===================================================================!

  subroutine check_representation_e_reverse_sparsity(nfail)

    integer, intent(inout) :: nfail

    type(picture) :: from_transpose, from_reverse_chain
    integer       :: k
    logical       :: agree

    from_transpose = sparsity_picture(d31t)

    call report(from_transpose % rows() .eq. 6 .and. &
         &      from_transpose % at(1) .eq. 'D3:1^T' .and. &
         &      from_transpose % at(2) .eq. '        m n' .and. &
         &      from_transpose % at(3) .eq. 'a       # .' .and. &
         &      from_transpose % at(4) .eq. 'b       # #' .and. &
         &      from_transpose % at(5) .eq. 'c       # #' .and. &
         &      from_transpose % at(6) .eq. 'd       . #', &
         & "D3:1^T renders the reverse skeleton: a from m; b and c " // &
         & "from both; d from n", nfail)

    from_reverse_chain = sparsity_picture(drev)

    agree = (from_transpose % rows() .eq. from_reverse_chain % rows())
    do k = 2, from_transpose % rows()
       agree = agree .and. &
            &  (from_transpose % at(k) .eq. from_reverse_chain % at(k))
    end do

    call report(agree, &
         & "AND D1^T o D2^T o D3^T RENDERS THE IDENTICAL GRID - " // &
         & "visualization respects structural transpose and " // &
         & "composition", nfail)

    call report(from_transpose % at(1) .ne. from_reverse_chain % at(1) .and. &
         &      same_extension(d31t, drev), &
         & "two relations, two names, one structure - the agreement " // &
         & "is about content, not identity", nfail)

  end subroutine check_representation_e_reverse_sparsity

  !===================================================================!
  ! The renderer invents nothing and forgets nothing: every cell of
  ! every picture is the relation's own answer, and the listing and
  ! the grid tell the same story because they ask the same question.
  !===================================================================!

  subroutine check_the_renderer_invents_nothing(nfail)

    integer, intent(inout) :: nfail

    call report(cells_match_membership(d1) .and. cells_match_membership(d2) .and. &
         &      cells_match_membership(d3) .and. cells_match_membership(d21) .and. &
         &      cells_match_membership(d31) .and. cells_match_membership(d31t), &
         & "every cell of every picture equals relation % has on that " // &
         & "pair - THE PICTURE IS DOWNSTREAM OF THE MATHEMATICS", nfail)

    call report(grid_agrees_with_listing(d31) .and. &
         &      grid_agrees_with_listing(d1t), &
         & "and the grid and the listing agree cell for cell: two " // &
         & "representations, one structure", nfail)

    ! The chain the renderer drew is a chain it checked: each leg's
    ! codomain is the next leg's domain, by identity.
    call report(legs_meet(d1, d2) .and. legs_meet(d2, d3) .and. &
         &      legs_meet(d3t, d2t) .and. legs_meet(d2t, d1t), &
         & "and the drawn chains genuinely compose - the renderer " // &
         & "refuses arrows between relations that do not meet", nfail)

  end subroutine check_the_renderer_invents_nothing

  !===================================================================!
  ! THE HOSTILE ORDER. A carrier that declares { 30, 10, 20 } must be
  ! drawn 30, 10, 20. Every other carrier in the specimen counts from
  ! one and would look right even to a renderer that sorted.
  !===================================================================!

  subroutine check_declaration_order_rules(nfail)

    integer, intent(inout) :: nfail

    type(counted_set)  :: ambient
    type(subset_set)   :: shuffled
    type(csr_relation) :: probe
    type(picture)      :: pic

    ambient  = counted_set('an ambient roster', 30)
    shuffled = subset_set('declared out of order', ambient, [30, 10, 20])

    probe = csr_relation('P', shuffled, x3, &
         &               reshape([30, X3_M, 20, X3_N], [2, 2]))

    pic = sparsity_picture(probe)

    call report(pic % at(2) .eq. '        30 10 20', &
         & "a carrier declaring { 30, 10, 20 } is drawn 30 10 20 - " // &
         & "NOT SORTED", nfail)

    call report(pic % at(3) .eq. 'm       #  .  .' .and. &
         &      pic % at(4) .eq. 'n       .  .  #', &
         & "and the glyphs follow the same order: 30 reaches m, 20 " // &
         & "reaches n, and the columns say so", nfail)

    call report(shuffled % member(1) .eq. 30 .and. &
         &      shuffled % local_index(30) .eq. 1 .and. &
         &      label_for(shuffled, 30) .eq. '30', &
         & "order came from member(i) and local_index, and the label " // &
         & "from a carrier the fixture has never heard of", nfail)

  end subroutine check_declaration_order_rules

  !===================================================================!
  ! THE GATE-A ARCHITECTURAL EXPERIMENT.
  !
  ! Question one: can the rectangular typed relation D1 : X0 -> X1 be
  ! rendered directly from the relational nucleus? Everything above
  ! this line is the answer, and it is yes.
  !
  ! Question two: would forcing D1 into ordinary_graph_view collapse
  ! X0 and X1, or otherwise change its meaning? The two profiles'
  ! contracts are inspected rather than provoked:
  !
  !   ordinary_graph_view  requires  T <= E x V  and  H <= E x V
  !                                  with ONE V, and refuses with
  !                                  'the head relation must share
  !                                  the tail''s domains' otherwise
  !
  !   directed_adjacency_view  requires  A <= V x V, and refuses with
  !                                  'a directed adjacency runs over
  !                                  one domain' otherwise
  !
  ! Both demand a SINGLE vertex carrier. The specimen's premises are
  ! checked here; the conclusion is Gate A's to state.
  !===================================================================!

  subroutine check_the_ordinary_graph_question(nfail)

    integer, intent(inout) :: nfail

    class(member_set), allocatable :: tail_end, head_end, from, to
    type(picture)                  :: pic

    ! ---- ANSWERED: the nucleus alone drew it.
    pic = sparsity_picture(d1)
    call report(pic % at(3) .eq. 'p       # # . .', &
         & "D1 : X0 -> X1 WAS RENDERED FROM member_set, relation and " // &
         & "graph-owned structure alone - no ordinary graph was " // &
         & "required at this radius", nfail)

    ! ---- The premise the ordinary profile's schema cannot satisfy.
    tail_end = t1 % target()
    head_end = h1 % target()
    call report(.not. tail_end % same_as(head_end), &
         & "T1 lands in X0 and H1 in X1, so no SINGLE vertex carrier " // &
         & "V exists for ordinary_graph_view to read them over", nfail)

    from = d1 % source()
    to   = d1 % target()
    call report(.not. from % same_as(to), &
         & "and D1's two ends are two domains, which is exactly what " // &
         & "directed_adjacency_view refuses", nfail)

    ! ---- What a union would cost, stated exactly.
    !
    !      DIRECTION IS NOT WHAT IS LOST. An ordinary directed graph
    !      preserves ordered-pair direction perfectly well - the tuple
    !      order of a same-domain relation IS the direction, as
    !      directed_adjacency_view documents.
    !
    !      What a union loses is the intrinsic distinct TYPED source
    !      and codomain: D : X_i -> X_j against D^T : X_j -> X_i. Under
    !      a union both become relations over one V, and the first
    !      thing same_extension tests - domain identity - has nothing
    !      left to compare.
    from = d31 % source()
    to   = d31 % target()
    call report(.not. from % same_as(to) .and. .not. same_extension(d31, d31t), &
         & "D3:1 and D3:1^T are told apart by their TYPED CARRIERS, " // &
         & "which a union would erase - direction would survive as " // &
         & "tuple order; the distinct source and codomain would not", nfail)

    call report(num_member_sets(g) .eq. 7 .and. .not. holds_set(g, bnd, union_like()), &
         & "so nothing was collapsed to make a picture: seven typed " // &
         & "carriers went in, and seven came out", nfail)

  end subroutine check_the_ordinary_graph_question

  !===================================================================!
  ! Helpers.
  !===================================================================!

  type(counted_set) function union_like()

    union_like = counted_set('V = X0 u X1 u X2 u X3', 12)

  end function union_like

  logical function legs_meet(first, second)

    class(binary_relation), intent(in) :: first, second

    class(member_set), allocatable :: here, there

    here  = first % target()
    there = second % source()
    legs_meet = here % same_as(there)

  end function legs_meet

  integer function filled_cells(r)

    class(relation), intent(in) :: r

    class(member_set), allocatable :: cols, rows
    integer                        :: i, j

    cols = r % domain(1)
    rows = r % domain(2)

    filled_cells = 0
    do j = 1, cols % size()
       do i = 1, rows % size()
          if (glyph_at(r, cols % member(j), rows % member(i)) .eq. '#') then
             filled_cells = filled_cells + 1
          end if
       end do
    end do

  end function filled_cells

  !-------------------------------------------------------------------!
  ! Every cell of the drawn grid, read back out of the picture and
  ! held against the relation's own membership test.
  !-------------------------------------------------------------------!

  logical function cells_match_membership(r)

    class(relation), intent(in) :: r

    class(member_set), allocatable :: cols, rows
    type(picture)                  :: pic
    character(len=:), allocatable  :: line
    integer                        :: i, j, stub, wide, at
    logical                        :: drawn

    cols = r % domain(1)
    rows = r % domain(2)
    pic  = sparsity_picture(r)

    stub = max(label_width(rows), 7) + 2
    wide = label_width(cols) + 1

    cells_match_membership = .true.
    do i = 1, rows % size()
       line = pic % at(2 + i)
       do j = 1, cols % size()
          at    = stub + (j - 1) * wide
          drawn = (line(at:at) .eq. '#')
          cells_match_membership = cells_match_membership .and. &
               & (drawn .eqv. r % has([cols % member(j), rows % member(i)]))
       end do
    end do

  end function cells_match_membership

  !-------------------------------------------------------------------!
  ! The listing says exactly what the grid says: for each member of
  ! the domain, the labels named in the listing are the members whose
  ! cell is filled.
  !-------------------------------------------------------------------!

  logical function grid_agrees_with_listing(r)

    class(relation), intent(in) :: r

    class(member_set), allocatable :: cols, rows
    type(picture)                  :: listing
    character(len=:), allocatable  :: said, wanted
    integer                        :: i, j

    cols    = r % domain(1)
    rows    = r % domain(2)
    listing = dependency_listing(r)

    grid_agrees_with_listing = .true.
    do j = 1, cols % size()
       wanted = label_for(cols, cols % member(j)) // ' ->'
       do i = 1, rows % size()
          if (glyph_at(r, cols % member(j), rows % member(i)) .eq. '#') then
             wanted = wanted // ' ' // label_for(rows, rows % member(i))
          end if
       end do
       said = listing % at(1 + j)
       grid_agrees_with_listing = grid_agrees_with_listing .and. &
            & (said .eq. wanted)
    end do

  end function grid_agrees_with_listing

  integer function label_width(carrier)

    class(member_set), intent(in) :: carrier

    integer :: k

    label_width = 1
    do k = 1, carrier % size()
       label_width = max(label_width, &
            &            len(label_for(carrier, carrier % member(k))))
    end do

  end function label_width

end program visualization_level_4
