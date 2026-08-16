!=====================================================================!
! VISUALIZATION TOWER . LEVEL 6 . DISCRETIZATION
!
! The level answers one question: DOES A PRODUCTION DISCRETIZATION
! OPERATOR EXPOSE THE SAME STRUCTURAL SKELETON ALREADY DERIVED
! RELATIONALLY.
!
! This is the first level to name production machinery, and the first
! executable consumer that
!
!      discretization_operator % dependencies()
!
! has ever had - the Level-6 census found ZERO callers of it anywhere
! in the repository.
!
!                    THE MEASUREMENT, NOT THE VERDICT
!
! Level 6 is an architectural measurement. A demonstrated
! specialization is a valid result, and the level is written so that
! "production is narrower than the tower" and "production is wrong"
! cannot be confused with each other.
!
! Nothing is applied. No stencil and no step is ever evaluated; both
! are CONSTRUCTED and INTERROGATED, and the import gate refuses an
! apply call in any tower source so that claim is mechanical.
!
!                        THE THREE MEASUREMENTS
!
! 1. THE STENCIL WITNESS. A production stencil is built with D2's
!    Boolean occupancy in production's own coordinates. Its
!    dependencies() renders an identical grid to the relational
!    D2 : X1 -> X2. Then the typed question is asked, and answered
!    differently:
!
!        same pixels     =>  YES
!        same object     =>  NO
!
!    because production's answer stands on ONE vertex carrier where
!    the relation stands on two, and X1 is not X2 even though both
!    hold three members.
!
! 2. THE RECTANGULAR WITNESS. D1 : X0 -> X1 runs 4 -> 3. Production's
!    constructor takes a single vertex count, so |V| would have to be
!    4 and 3 at once. Given 4 it can hold every edge - and it then
!    reports a FOURTH row that the typed codomain does not have. That
!    phantom row is the whole finding, and it is read off safely: no
!    out-of-range index, no manufactured union carrier.
!
! 3. THE STEP WITNESS. The same family verb, asked of the other
!    concrete citizen, answers on a different AXIS: the stencil on the
!    independent variable. BDF2's residual at the newest instant reads
!    three instants, so its stencil is a FAN-IN - 1->3, 2->3, 3->3 -
!    and not the succession 1->2->3, which is a chronology and a
!    different relation entirely.
!
!    Two steps wrapping actions with DIFFERENT state sparsity return
!    the SAME stencil, which is executable evidence that a step's
!    dependencies() describes its own axis and not the wrapped
!    action's algebra.
!
!    ONE CONTRACT, ONE MEANING:
!
!        dependencies() = the stencil on the axis this concrete type
!                         represents
!
!    dependent for the stencil operator, independent for the step.
!
!                     THE THEOREM THE LEVEL EXISTS FOR
!
!      V(R1) = V(R2)   does NOT imply   R1 = R2
!
! when the representation V forgets carrier identity. Hence every
! picture below carries its SIGNATURE line: the grid alone is exactly
! the part that cannot tell the two apart.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program visualization_level_6

  use visualization_assert , only : report, verdict
  use visualization_assert , only : ND1, ND2
  use visualization_assert , only : NX0, NX1, NX2
  use visualization_assert , only : X0_A, X0_B, X0_C, X0_D
  use visualization_assert , only : X1_P, X1_Q, X1_R
  use visualization_assert , only : X2_U, X2_V, X2_W
  use fractal_graph        , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation
  use graph_set_map        , only : set_map
  use graph_label_map      , only : label_map
  use graph_relation       , only : relation
  use graph_binary_relation, only : csr_relation
  use graph_grammar        , only : graph
  use class_graph_stencil  , only : stencil_operator
  use class_graph_step     , only : step_operator
  use visualization_carriers_fixture , only : structural_carriers
  use visualization_relations_fixture, only : occurrences_of_a1
  use visualization_relations_fixture, only : occurrences_of_a2
  use visualization_relations_fixture, only : occurrences_of_a3
  use visualization_algebra_fixture  , only : derive_dependency
  use structural_renderer_fixture    , only : picture, sparsity_picture
  use production_discretization_fixture, only : d2_coordinate_stencil
  use production_discretization_fixture, only : d1_coordinate_stencil
  use production_discretization_fixture, only : diagonal_stencil, bdf2_around
  use production_pattern_renderer_fixture, only : pattern_picture
  use production_pattern_renderer_fixture, only : production_has
  use production_pattern_renderer_fixture, only : same_coordinate_pattern
  use production_pattern_renderer_fixture, only : coordinate_shapes_fit
  use production_pattern_renderer_fixture, only : signature_of_pattern
  use production_pattern_renderer_fixture, only : signature_of_relation

  implicit none

  type(set_graph)      :: x0, x1, x2, x3, e1, e2, e3
  type(set_map)     :: sets
  type(label_map)     :: labels
  type(csr_relation)     :: t1, h1, t2, h2, t3, h3
  type(csr_relation)     :: d1, d2, d3

  type(stencil_operator) :: sten_d2, sten_d1, sten_diag
  type(step_operator)    :: clock_d2, clock_diag

  class(graph), allocatable :: pat_d2, pat_d1, pat_diag
  class(graph), allocatable :: motif_d2, motif_diag

  integer :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "visualization tower . level 6 . discretization"
  write(*,'(1x,a)') "============================================="

  ! ---- the tower's own structure, untouched since Gate A.
  call structural_carriers(x0, x1, x2, x3, e1, e2, e3, sets, labels)
  call occurrences_of_a1(e1, x0, x1, t1, h1, sets)
  call occurrences_of_a2(e2, x1, x2, t2, h2, sets)
  call occurrences_of_a3(e3, x2, x3, t3, h3, sets)

  d1 = derive_dependency('D1', t1, h1, sets)
  d2 = derive_dependency('D2', t2, h2, sets)
  d3 = derive_dependency('D3', t3, h3, sets)

  ! ---- production. CONSTRUCTED AND INTERROGATED, NEVER APPLIED.
  sten_d2   = d2_coordinate_stencil()
  sten_d1   = d1_coordinate_stencil()
  sten_diag = diagonal_stencil()

  call sten_d2 % dependencies(pat_d2)
  if (.not. sets % describes(pat_d2 % vertex_set())) &
       & call sets % bind(pat_d2 % vertex_set(), &
       &      counted_set_representation(pat_d2 % num_vertices()))
  if (.not. sets % describes(pat_d2 % edge_set())) &
       & call sets % bind(pat_d2 % edge_set(), &
       &      counted_set_representation(pat_d2 % num_edges()))
  ! the production graph names its own two domains
  if (.not. labels % labelled(pat_d2 % vertex_set())) &
       & call labels % bind(pat_d2 % vertex_set(), 'vertices')
  if (.not. labels % labelled(pat_d2 % edge_set())) &
       & call labels % bind(pat_d2 % edge_set(), 'edges')
  call sten_d1 % dependencies(pat_d1)
  if (.not. sets % describes(pat_d1 % vertex_set())) &
       & call sets % bind(pat_d1 % vertex_set(), &
       &      counted_set_representation(pat_d1 % num_vertices()))
  if (.not. sets % describes(pat_d1 % edge_set())) &
       & call sets % bind(pat_d1 % edge_set(), &
       &      counted_set_representation(pat_d1 % num_edges()))
  ! the production graph names its own two domains
  if (.not. labels % labelled(pat_d1 % vertex_set())) &
       & call labels % bind(pat_d1 % vertex_set(), 'vertices')
  if (.not. labels % labelled(pat_d1 % edge_set())) &
       & call labels % bind(pat_d1 % edge_set(), 'edges')
  call sten_diag % dependencies(pat_diag)
  if (.not. sets % describes(pat_diag % vertex_set())) &
       & call sets % bind(pat_diag % vertex_set(), &
       &      counted_set_representation(pat_diag % num_vertices()))
  if (.not. sets % describes(pat_diag % edge_set())) &
       & call sets % bind(pat_diag % edge_set(), &
       &      counted_set_representation(pat_diag % num_edges()))
  ! the production graph names its own two domains
  if (.not. labels % labelled(pat_diag % vertex_set())) &
       & call labels % bind(pat_diag % vertex_set(), 'vertices')
  if (.not. labels % labelled(pat_diag % edge_set())) &
       & call labels % bind(pat_diag % edge_set(), 'edges')

  clock_d2   = bdf2_around(sten_d2)
  clock_diag = bdf2_around(sten_diag)

  call clock_d2 % dependencies(motif_d2)
  if (.not. sets % describes(motif_d2 % vertex_set())) &
       & call sets % bind(motif_d2 % vertex_set(), &
       &      counted_set_representation(motif_d2 % num_vertices()))
  if (.not. sets % describes(motif_d2 % edge_set())) &
       & call sets % bind(motif_d2 % edge_set(), &
       &      counted_set_representation(motif_d2 % num_edges()))
  ! the production graph names its own two domains
  if (.not. labels % labelled(motif_d2 % vertex_set())) &
       & call labels % bind(motif_d2 % vertex_set(), 'vertices')
  if (.not. labels % labelled(motif_d2 % edge_set())) &
       & call labels % bind(motif_d2 % edge_set(), 'edges')
  call clock_diag % dependencies(motif_diag)
  if (.not. sets % describes(motif_diag % vertex_set())) &
       & call sets % bind(motif_diag % vertex_set(), &
       &      counted_set_representation(motif_diag % num_vertices()))
  if (.not. sets % describes(motif_diag % edge_set())) &
       & call sets % bind(motif_diag % edge_set(), &
       &      counted_set_representation(motif_diag % num_edges()))
  ! the production graph names its own two domains
  if (.not. labels % labelled(motif_diag % vertex_set())) &
       & call labels % bind(motif_diag % vertex_set(), 'vertices')
  if (.not. labels % labelled(motif_diag % edge_set())) &
       & call labels % bind(motif_diag % edge_set(), 'edges')

  call say_the_measurement()

  call check_the_stencil_coordinate_pattern(nfail)
  call check_the_typed_identity_differs(nfail)
  call check_the_rectangular_witness(nfail)
  call check_the_step_pattern(nfail)
  call check_state_is_not_time(nfail)
  call check_the_motif_is_independent_of_what_it_wraps(nfail)
  call check_the_visual_equality_theorem(nfail)

  call verdict(nfail, "level 6")

contains

  !===================================================================!
  ! The measurement, printed for a person. Every grid below was
  ! generated a moment ago - the left ones from twelve occurrences,
  ! the right ones from production's own dependencies().
  !===================================================================!

  subroutine say_the_measurement()

    type(picture)     :: pic
    type(set_graph) :: verts

    write(*,'(1x,a)') "---------------------------------------------"

    write(*,'(4x,a)') "RELATIONAL D2"
    write(*,'(4x,a)') "signature: " // signature_of_relation(d2, labels)
    pic = sparsity_picture(d2, sets, labels); call say_grid(pic, 2)

    write(*,'(4x,a)') "STENCIL dependencies()"
    pic = pattern_picture(pat_d2, '', sets, labels); call say_grid(pic, 2)

    verts = carrier_of(pat_d2)
    call say_verdict("coordinate pattern equal", &
         &           same_coordinate_pattern(d2, pat_d2, sets))
    call say_verdict("typed source identity equal", verts % same_as(x1))
    call say_verdict("typed target identity equal", verts % same_as(x2))
    write(*,*)

    write(*,'(4x,a)') "RECTANGULAR D1"
    write(*,'(4x,a)') "signature: " // signature_of_relation(d1, labels)
    write(*,'(4x,a,i0,a,i0,a)') "shape: ", NX1, " rows x ", NX0, " columns"
    pic = sparsity_picture(d1, sets, labels); call say_grid(pic, 2)

    write(*,'(4x,a)') "STENCIL dependencies() for the same occupancy"
    pic = pattern_picture(pat_d1, '', sets, labels); call say_grid(pic, 2)

    call say_verdict("production contract preserves this typed signature", &
         &           coordinate_shapes_fit(d1, pat_d1, sets))
    write(*,*)

    write(*,'(4x,a)') "BDF2 dependencies()"
    pic = pattern_picture(motif_d2, '', sets, labels); call say_grid(pic, 2)

    call say_verdict("wrapped state pattern equal", &
         &           same_production_pattern(motif_d2, pat_d2))
    call say_verdict("same BDF2 motif under second wrapped sparsity", &
         &           same_production_pattern(motif_d2, motif_diag))

    write(*,'(1x,a)') "---------------------------------------------"

  end subroutine say_the_measurement

  subroutine say_grid(pic, from)

    type(picture), intent(in) :: pic
    integer      , intent(in) :: from

    integer :: k

    do k = from, pic % rows()
       write(*,'(4x,a)') pic % at(k)
    end do
    write(*,*)

  end subroutine say_grid

  subroutine say_verdict(label, yes)

    character(len=*), intent(in) :: label
    logical         , intent(in) :: yes

    character(len=48) :: dots

    dots = repeat('.', max(1, 48 - len(label)))
    if (yes) then
       write(*,'(4x,a,1x,a,1x,a)') label, trim(dots), "YES"
    else
       write(*,'(4x,a,1x,a,1x,a)') label, trim(dots), "NO"
    end if

  end subroutine say_verdict

  !===================================================================!
  ! MEASUREMENT ONE, first half. The Boolean occupancy agrees.
  !===================================================================!

  subroutine check_the_stencil_coordinate_pattern(nfail)

    integer, intent(inout) :: nfail

    type(picture) :: pic

    call report(pat_d2 % num_vertices() .eq. NX1 .and. &
         &      pat_d2 % num_edges() .eq. ND2, &
         & "production's dependencies() answers a graph of 3 vertices " // &
         & "and 4 edges - one edge per stencil coefficient", nfail)

    call report(coordinate_shapes_fit(d2, pat_d2, sets), &
         & "and its vertex count matches BOTH of D2's carriers, so " // &
         & "the two can be laid over one another at all", nfail)

    call report(same_coordinate_pattern(d2, pat_d2, sets), &
         & "SAME COORDINATE PATTERN as D2 : X1 -> X2, cell by cell, " // &
         & "columns and rows in declaration order", nfail)

    pic = pattern_picture(pat_d2, 'STENCIL', sets, labels)
    call report(pic % at(3) .eq. '        1 2 3' .and. &
         &      pic % at(4) .eq. '1       # # .' .and. &
         &      pic % at(5) .eq. '2       . # .' .and. &
         &      pic % at(6) .eq. '3       . . #', &
         & "and it renders the identical grid the tower drew at " // &
         & "Level 4, glyph for glyph", nfail)

    ! Read directly off production's own edge ends, not off the picture.
    call report(production_has(pat_d2, 1, 1) .and. &
         &      production_has(pat_d2, 2, 1) .and. &
         &      production_has(pat_d2, 2, 2) .and. &
         &      production_has(pat_d2, 3, 3) .and. &
         &      .not. production_has(pat_d2, 1, 2) .and. &
         &      .not. production_has(pat_d2, 3, 1), &
         & "read off production's own edge ends: column -> row, four " // &
         & "arrows and no others", nfail)

  end subroutine check_the_stencil_coordinate_pattern

  !===================================================================!
  ! MEASUREMENT ONE, second half. THE SAME PIXELS ARE NOT THE SAME
  ! OBJECT.
  !===================================================================!

  subroutine check_the_typed_identity_differs(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: from, to
    type(set_graph)              :: verts

    verts = carrier_of(pat_d2)
    from  = d2 % domain(1)
    to    = d2 % domain(2)

    call report(.not. from % same_as(to) .and. from % same_as(x1) .and. &
         &      to % same_as(x2), &
         & "D2 stands on TWO declared carriers: X1 -> X2, and X1 is " // &
         & "not X2 though both hold three members", nfail)

    call report(.not. verts % same_as(x1) .and. .not. verts % same_as(x2), &
         & "production's pattern stands on a carrier that is NEITHER " // &
         & "X1 NOR X2 - it declared its own", nfail)

    call report(signature_of_pattern(pat_d2, labels) .eq. 'vertices -> vertices' .and. &
         &      signature_of_relation(d2, labels) .eq. 'X1 -> X2', &
         & "ONE carrier in both places against TWO: the signatures " // &
         & "read 'vertices -> vertices' and 'X1 -> X2'", nfail)

    call report(sets % size_of(verts) .eq. NX1 .and. sets % size_of(verts) .eq. NX2 .and. &
         &      same_coordinate_pattern(d2, pat_d2, sets), &
         & "SAME PIXELS, NOT THE SAME TYPED STRUCTURAL OBJECT - equal " // &
         & "occupancy, equal counts, and three distinct carriers", nfail)

  end subroutine check_the_typed_identity_differs

  !===================================================================!
  ! MEASUREMENT TWO. The rectangular witness: 4 -> 3 against a single
  ! vertex count.
  !===================================================================!

  subroutine check_the_rectangular_witness(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: from, to
    type(picture)                  :: pic
    integer                        :: j
    logical                        :: phantom_empty

    from = d1 % domain(1)
    to   = d1 % domain(2)

    call report(sets % size_of(from) .eq. NX0 .and. sets % size_of(to) .eq. NX1 .and. &
         &      NX0 .ne. NX1, &
         & "D1 : X0 -> X1 is RECTANGULAR - four columns, three rows", &
         & nfail)

    call report(pat_d1 % num_vertices() .eq. NX0 .and. &
         &      pat_d1 % num_edges() .eq. ND1, &
         & "production held every one of D1's five arrows - given a " // &
         & "vertex count of 4, the only count that fits the columns", &
         & nfail)

    call report(.not. coordinate_shapes_fit(d1, pat_d1, sets), &
         & "AND THE SHAPES DO NOT FIT: |V| = 4 would have to be 3 as " // &
         & "well, and one number cannot be both", nfail)

    ! The phantom row: vertex 4 exists as a HEAD position because the
    ! carrier is one set, though D1's codomain has no fourth member.
    phantom_empty = .true.
    do j = 1, pat_d1 % num_vertices()
       phantom_empty = phantom_empty .and. .not. production_has(pat_d1, j, 4)
    end do

    pic = pattern_picture(pat_d1, 'STENCIL', sets, labels)
    call report(pic % rows() .eq. 3 + NX0 .and. phantom_empty .and. &
         &      pic % at(7) .eq. '4       . . . .', &
         & "so the picture has a FOURTH ROW that D1's codomain does " // &
         & "not have - empty, and present only because one carrier " // &
         & "serves both axes", nfail)

    call report(pic % at(4) .eq. '1       # # . .' .and. &
         &      pic % at(5) .eq. '2       . # # .' .and. &
         &      pic % at(6) .eq. '3       . . . #', &
         & "its first three rows do carry D1's occupancy - the " // &
         & "arrows survived; the SIGNATURE did not", nfail)

    call report(.not. same_coordinate_pattern(d1, pat_d1, sets), &
         & "and same_coordinate_pattern refuses the comparison rather " // &
         & "than padding D1 to make it fit", nfail)

  end subroutine check_the_rectangular_witness

  !===================================================================!
  ! MEASUREMENT THREE. The same family verb, asked of the other
  ! concrete citizen.
  !===================================================================!

  subroutine check_the_step_pattern(nfail)

    integer, intent(inout) :: nfail

    type(picture) :: pic

    call report(motif_d2 % num_vertices() .eq. 3 .and. &
         &      motif_d2 % num_edges() .eq. 3, &
         & "BDF2's dependencies() answers three vertices and three " // &
         & "edges - reach + 1 instants, and one arrow from each", nfail)

    call report(production_has(motif_d2, 1, 3) .and. &
         &      production_has(motif_d2, 2, 3) .and. &
         &      production_has(motif_d2, 3, 3) .and. &
         &      .not. production_has(motif_d2, 1, 2) .and. &
         &      .not. production_has(motif_d2, 2, 1), &
         & "and its arrows FAN IN on the newest instant: 1->3, 2->3, " // &
         & "3->3 - what the residual actually reads", nfail)

    pic = pattern_picture(motif_d2, 'BDF2', sets, labels)
    call report(pic % at(4) .eq. '1       . . .' .and. &
         &      pic % at(5) .eq. '2       . . .' .and. &
         &      pic % at(6) .eq. '3       # # #', &
         & "rendered, one full row and nothing else - THE STENCIL OF " // &
         & "THE INSTANT BEING SOLVED FOR", nfail)

    ! A stencil is not a chronology. The succession 1 -> 2 -> 3 is a
    ! true relation about which instant follows which, and it is not
    ! what this contract answers.
    call report(.not. production_has(motif_d2, 1, 2) .and. &
         &      production_has(motif_d2, 3, 3), &
         & "A STENCIL IS NOT A CHRONOLOGY: succession would hold " // &
         & "1->2, and the self-arrow that makes the newest instant " // &
         & "an unknown would be missing", nfail)

  end subroutine check_the_step_pattern

  !===================================================================!
  ! STATE DEPENDENCY IS NOT TEMPORAL DEPENDENCY. Both answers have
  ! three vertices here, and that is the trap.
  !===================================================================!

  subroutine check_state_is_not_time(nfail)

    integer, intent(inout) :: nfail

    call report(motif_d2 % num_vertices() .eq. pat_d2 % num_vertices(), &
         & "the step's answer and its wrapped action's answer have " // &
         & "the SAME VERTEX COUNT - three each", nfail)

    call report(.not. same_production_pattern(motif_d2, pat_d2), &
         & "AND DIFFERENT EDGE STRUCTURE - equal cardinality is not " // &
         & "equal semantics", nfail)

    call report(production_has(pat_d2, 2, 1) .and. &
         &      .not. production_has(motif_d2, 2, 1), &
         & "the state stencil couples unknown 2 to unknown 1; the " // &
         & "step's independent axis has no such arrow", nfail)

    call report(production_has(motif_d2, 1, 3) .and. &
         &      .not. production_has(pat_d2, 1, 3), &
         & "the step's residual reads instant 1; the state stencil " // &
         & "has no coefficient coupling unknown 1 to unknown 3", nfail)

  end subroutine check_state_is_not_time

  !===================================================================!
  ! THE INDEPENDENCE PROBE. Two wrapped actions with genuinely
  ! different state sparsity, one motif.
  !===================================================================!

  subroutine check_the_motif_is_independent_of_what_it_wraps(nfail)

    integer, intent(inout) :: nfail

    call report(pat_diag % num_vertices() .eq. pat_d2 % num_vertices() .and. &
         &      .not. same_production_pattern(pat_diag, pat_d2), &
         & "the second action's state sparsity is genuinely different " // &
         & "- diagonal against D2's occupancy, same three vertices", &
         & nfail)

    call report(same_production_pattern(motif_d2, motif_diag), &
         & "YET BDF2 RETURNS THE IDENTICAL MOTIF FOR BOTH - so the " // &
         & "step's dependencies() does not describe what it wraps", &
         & nfail)

    call report(motif_diag % num_vertices() .eq. 3 .and. &
         &      production_has(motif_diag, 1, 3) .and. &
         &      production_has(motif_diag, 3, 3), &
         & "it describes THE STENCIL ON ITS OWN AXIS: reach + 1 " // &
         & "instants fanning in, whatever action rides on them", nfail)

  end subroutine check_the_motif_is_independent_of_what_it_wraps

  !===================================================================!
  ! THE THEOREM.
  !
  !      V(R1) = V(R2)  does not imply  R1 = R2
  !
  ! when V forgets carrier identity. Both halves are checked: the
  ! grids agree, and the objects do not.
  !===================================================================!

  subroutine check_the_visual_equality_theorem(nfail)

    integer, intent(inout) :: nfail

    type(picture)     :: relational, production
    type(set_graph) :: verts
    integer           :: k
    logical           :: grids_agree

    relational = sparsity_picture(d2, sets, labels)
    production = pattern_picture(pat_d2, 'STENCIL', sets, labels)

    ! The relational picture is name + header + rows; the production
    ! one is title + signature + header + rows. Compare the ROWS,
    ! stripped of their labels and their spacing - which is all a
    ! reader sees when no signature is shown.
    grids_agree = .true.
    do k = 1, NX2
       grids_agree = grids_agree .and. &
            & (glyphs_of(relational % at(2 + k)) .eq. &
            &  glyphs_of(production % at(3 + k)))
       grids_agree = grids_agree .and. &
            & (len(glyphs_of(relational % at(2 + k))) .eq. NX1)
    end do

    call report(grids_agree, &
         & "THE TWO GRIDS AGREE glyph for glyph, though one was drawn " // &
         & "from a typed relation and the other from a production " // &
         & "graph", nfail)

    verts = carrier_of(pat_d2)
    call report(.not. verts % same_as(x1) .and. .not. verts % same_as(x2) .and. &
         &      .not. x1 % same_as(x2), &
         & "AND THE OBJECTS DO NOT: three declared carriers stand " // &
         & "where the pictures show one shape", nfail)

    call report(signature_of_relation(d2, labels) .ne. signature_of_pattern(pat_d2, labels), &
         & "which is why every Level-6 picture carries its SIGNATURE " // &
         & "- the grid alone is exactly the part that cannot tell " // &
         & "them apart", nfail)

  end subroutine check_the_visual_equality_theorem

  !===================================================================!
  ! Helpers.
  !===================================================================!

  type(set_graph) function carrier_of(p)

    class(graph), intent(in) :: p

    carrier_of = p % vertex_set()

  end function carrier_of

  !-------------------------------------------------------------------!
  ! Two production patterns, compared extensionally: same vertex
  ! count, and the same arrow at every pair of local coordinates.
  !-------------------------------------------------------------------!

  logical function same_production_pattern(p, q)

    class(graph), intent(in) :: p, q

    integer :: i, j

    same_production_pattern = (p % num_vertices() .eq. q % num_vertices())
    if (.not. same_production_pattern) return

    do j = 1, p % num_vertices()
       do i = 1, p % num_vertices()
          same_production_pattern = same_production_pattern .and. &
               & (production_has(p, j, i) .eqv. production_has(q, j, i))
       end do
    end do

  end function same_production_pattern

  !-------------------------------------------------------------------!
  ! A rendered row's glyphs alone, with its label and spacing thrown
  ! away - what a reader sees when the signature is not shown.
  !-------------------------------------------------------------------!

  function glyphs_of(line) result(marks)

    character(len=*), intent(in) :: line

    character(len=:), allocatable :: marks
    integer                       :: k

    marks = ''
    do k = 1, len(line)
       if (line(k:k) .eq. '#' .or. line(k:k) .eq. '.') then
          marks = marks // line(k:k)
       end if
    end do

  end function glyphs_of

end program visualization_level_6
