!=====================================================================!
! GTI ADAPTIVE GRAPH GROWTH
!
! The first adaptive-composition layer: the time graph grows one
! candidate at a time, transactionally -
!
!      candidate vertex + candidate relation
!        -> append to G_time
!        -> solve the appended relation
!        -> accept keeps the graph
!        -> reject rolls the graph back.
!
! Adaptive composition IS graph growth: G_n = (S_n, R_n) proposes
! S_{n+1} = S_n u {new sample}, R_{n+1} = R_n u {new relation},
! and the only new law is transactional mutation - an accepted
! append persists, a rejected append leaves no trace, and a later
! accepted relation may read an earlier accepted vertex as
! history. The local solve is the unchanged forward driver's
! solve_relation; nothing here estimates an error or chooses a
! step size - the accept decision arrives as a logical from an
! EXTERNAL controller, and a later phase may automate it.
!
! A candidate names its own vertex by the sentinel -1 in its
! relation tuple: "the vertex being appended by this candidate."
! The candidate laws: the relation must reference the appended
! vertex, and the appended vertex must be its unknown - a
! candidate that solves something else is not a growth step.
!
! A non-converged candidate is an answer: it reports solved and
! unconverged, is rolled back whole, and the graph stands exactly
! as before - the forward driver never wrote, and the appended
! seats are gone.
!
! The driver carries nothing: no graph, no controller, no
! estimator, no forms, no state.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_time_adaptive_growth_drivers

  use gti_design_bundles   , only : gti_design_bundle
  use gti_form_interface   , only : gti_differentiable_form
  use gti_time_graphs      , only : gti_time_graph, gti_time_vertex, &
       & gti_time_relation
  use gti_time_forward_drivers, only : gti_time_forward_driver, &
       & gti_time_forward_options, gti_time_forward_step_result

  implicit none

  private
  public :: gti_time_growth_candidate
  public :: gti_time_growth_step_result
  public :: gti_time_adaptive_growth_driver

  !===================================================================!
  ! One proposed growth step: the vertex to append, and the
  ! relation that solves it. The sentinel -1 in the relation's
  ! vertex tuple names the vertex being appended.
  !===================================================================!

  type :: gti_time_growth_candidate

     type(gti_time_vertex)   :: vertex
     type(gti_time_relation) :: relation

  end type gti_time_growth_candidate

  !===================================================================!
  ! What one growth step reports: whether it was kept, whether
  ! and how its solve went, where it appended, and the graph
  ! sizes it left behind.
  !===================================================================!

  type :: gti_time_growth_step_result

     logical :: accepted = .false.
     logical :: solved = .false.
     logical :: converged = .false.
     integer :: appended_vertex = 0
     integer :: appended_relation = 0
     integer :: num_vertices_after = 0
     integer :: num_relations_after = 0
     type(gti_time_forward_step_result) :: forward_step

  end type gti_time_growth_step_result

  !===================================================================!
  ! The stateless growth verb-set. The types keep their public
  ! singular names; Fortran denies a type its host module's name,
  ! so the module speaks in the plural.
  !===================================================================!

  type :: gti_time_adaptive_growth_driver

   contains

     procedure :: append_candidate
     procedure :: discard_last_candidate
     procedure :: try_candidate

  end type gti_time_adaptive_growth_driver

contains

  !===================================================================!
  ! Append one candidate: the vertex becomes the last vertex, the
  ! relation - its sentinel resolved to the new index - becomes
  ! the last relation, the candidate laws are proven, and the
  ! grown graph is validated whole.
  !===================================================================!

  subroutine append_candidate(this, graph, candidate, appended_vertex, &
       & appended_relation)

    class(gti_time_adaptive_growth_driver), intent(in)    :: this
    type(gti_time_graph)                  , intent(inout) :: graph
    type(gti_time_growth_candidate)       , intent(in)    :: candidate
    integer                               , intent(out)   :: appended_vertex
    integer                               , intent(out)   :: appended_relation

    type(gti_time_relation) :: relation
    integer :: k
    logical :: references_appended

    if (graph % num_vertices() < 1) then
       error stop 'gti_time_adaptive_growth_driver: graph has an initial vertex'
    end if

    call append_vertex(graph, candidate % vertex)
    appended_vertex = graph % num_vertices()

    !----------------------------------------------------------------!
    ! Resolve the sentinel: -1 means the vertex just appended. A
    ! candidate whose relation never names its own vertex is not
    ! a growth step.
    !----------------------------------------------------------------!

    relation = candidate % relation

    references_appended = .false.
    if (allocated(relation % sample_vertex)) then
       do k = 1, size(relation % sample_vertex)
          if (relation % sample_vertex(k) == -1) then
             relation % sample_vertex(k) = appended_vertex
             references_appended = .true.
          end if
       end do
    end if

    if (.not. references_appended) then
       error stop 'gti_time_adaptive_growth_driver: candidate relation references appended vertex'
    end if

    call append_relation(graph, relation)
    appended_relation = graph % num_relations()

    if (graph % relation(appended_relation) % unknown_vertex() /= &
         & appended_vertex) then
       error stop 'gti_time_adaptive_growth_driver: candidate unknown is appended vertex'
    end if

    call graph % validate()

  end subroutine append_candidate

  !===================================================================!
  ! Roll back the most recent append: the last relation and the
  ! last vertex leave, and everything earlier - vertices,
  ! relations, q values, solution flags - stands exactly as it
  ! was. The expected indices are the caller's proof it is
  ! discarding what it appended.
  !===================================================================!

  subroutine discard_last_candidate(this, graph, expected_vertex, &
       & expected_relation)

    class(gti_time_adaptive_growth_driver), intent(in)    :: this
    type(gti_time_graph)                  , intent(inout) :: graph
    integer                               , intent(in)    :: expected_vertex
    integer                               , intent(in)    :: expected_relation

    if (expected_vertex /= graph % num_vertices()) then
       error stop 'gti_time_adaptive_growth_driver: appended vertex is last'
    end if

    if (expected_relation /= graph % num_relations()) then
       error stop 'gti_time_adaptive_growth_driver: appended relation is last'
    end if

    call shrink_relations(graph, graph % num_relations() - 1)
    call shrink_vertices(graph, graph % num_vertices() - 1)

    if (graph % num_relations() >= 1) then
       call graph % validate()
    end if

  end subroutine discard_last_candidate

  !===================================================================!
  ! One transactional growth step: append, solve the appended
  ! relation alone, and keep or roll back. Non-convergence rolls
  ! back and reports - an answer, never an error stop. The accept
  ! logical is the external controller's decision, nothing more.
  !===================================================================!

  subroutine try_candidate(this, residual_form, graph, candidate, design, &
       & options, accept, result)

    class(gti_time_adaptive_growth_driver), intent(in)    :: this
    class(gti_differentiable_form)        , intent(in)    :: residual_form
    type(gti_time_graph)                  , intent(inout) :: graph
    type(gti_time_growth_candidate)       , intent(in)    :: candidate
    type(gti_design_bundle)               , intent(in)    :: design
    type(gti_time_forward_options)        , intent(in)    :: options
    logical                               , intent(in)    :: accept
    type(gti_time_growth_step_result)     , intent(inout) :: result

    type(gti_time_forward_driver) :: forward
    integer :: appended_vertex, appended_relation

    result % accepted  = .false.
    result % solved    = .false.
    result % converged = .false.
    result % appended_vertex = 0
    result % appended_relation = 0
    result % forward_step = gti_time_forward_step_result()

    call this % append_candidate(graph, candidate, appended_vertex, &
         & appended_relation)

    result % appended_vertex   = appended_vertex
    result % appended_relation = appended_relation

    call forward % solve_relation(residual_form, graph, appended_relation, &
         & design, options, result % forward_step)

    result % solved    = .true.
    result % converged = result % forward_step % converged

    if (.not. result % converged) then
       call this % discard_last_candidate(graph, appended_vertex, &
            & appended_relation)
    else if (accept) then
       result % accepted = .true.
    else
       call this % discard_last_candidate(graph, appended_vertex, &
            & appended_relation)
    end if

    result % num_vertices_after  = graph % num_vertices()
    result % num_relations_after = graph % num_relations()

  end subroutine try_candidate

  !===================================================================!
  ! The private growth mechanics: grow and shrink the graph's
  ! arrays while preserving every earlier seat exactly.
  !===================================================================!

  subroutine append_vertex(graph, vertex)

    type(gti_time_graph) , intent(inout) :: graph
    type(gti_time_vertex), intent(in)    :: vertex

    type(gti_time_vertex), allocatable :: grown(:)
    integer :: n

    n = graph % num_vertices()
    allocate(grown(n + 1))
    if (n > 0) grown(1:n) = graph % vertex(1:n)
    grown(n + 1) = vertex

    call move_alloc(grown, graph % vertex)

  end subroutine append_vertex

  subroutine append_relation(graph, relation)

    type(gti_time_graph)   , intent(inout) :: graph
    type(gti_time_relation), intent(in)    :: relation

    type(gti_time_relation), allocatable :: grown(:)
    integer :: n

    n = graph % num_relations()
    allocate(grown(n + 1))
    if (n > 0) grown(1:n) = graph % relation(1:n)
    grown(n + 1) = relation

    call move_alloc(grown, graph % relation)

  end subroutine append_relation

  subroutine shrink_vertices(graph, n)

    type(gti_time_graph), intent(inout) :: graph
    integer             , intent(in)    :: n

    type(gti_time_vertex), allocatable :: kept(:)

    allocate(kept(n))
    if (n > 0) kept(1:n) = graph % vertex(1:n)

    call move_alloc(kept, graph % vertex)

  end subroutine shrink_vertices

  subroutine shrink_relations(graph, n)

    type(gti_time_graph), intent(inout) :: graph
    integer             , intent(in)    :: n

    type(gti_time_relation), allocatable :: kept(:)

    allocate(kept(n))
    if (n > 0) kept(1:n) = graph % relation(1:n)

    call move_alloc(kept, graph % relation)

  end subroutine shrink_relations

end module gti_time_adaptive_growth_drivers
