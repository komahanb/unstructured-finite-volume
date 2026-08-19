!=====================================================================!
! GTI TIME GRAPH REPRESENTATION
!
! The first explicit time graph: the relation view
!
!      G_time = (S_time, R_time)
!
! stored as plain arrays. S_time is the set of time sample
! vertices - each carrying one sample (a q field and its time) and
! a solution flag. R_time is the set of local motif applications -
! each an ordered tuple of sample vertices, one unknown position
! among them, one coefficient motif, and one evaluation time. A
! vertex shared by two relations is how one step's solution
! becomes the next step's history: r1 writes vertex 2, r2 reads it.
!
! This seat stores, validates, and materializes. It does not
! march, does not traverse, does not solve: no Newton, tangent, or
! adjoint driver is named here, and the local calculus consumes
! what build_samples hands it without ever seeing the graph. The
! kernel graph is deliberately absent too - this is a
! representation of the GTI time relation, not a rival ontology,
! and identity questions stay where they live.
!
! The laws validate proves before anything is materialized:
!
!      at least one vertex; every relation has vertices, all in
!      range, its unknown among them, a motif with rules, every
!      rule weighted once per sample, and every referenced vertex
!      providing q.
!
! Writing a solution back replaces VALUES, never identity: the
! unknown vertex's q is minted from its own field through the
! unknown-problem injection, so the domain rides unchanged, the
! sample time and non-q seats survive, and the vertex is marked
! solved.
!
! The graph carries vertices and relations, and nothing else: no
! form, no solver, no design, no schedule, no map.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_time_graphs

  use iso_fortran_env      , only : dp => REAL64
  use gti_value_buffers    , only : gti_value_buffer
  use gti_state_bundles    , only : GTI_STATE_Q
  use gti_time_local_schemes, only : gti_time_sample, gti_time_motif
  use gti_time_local_unknown_problems, only : gti_time_local_unknown_problem

  implicit none

  private
  public :: gti_time_vertex
  public :: gti_time_relation
  public :: gti_time_graph

  !===================================================================!
  ! One time sample vertex: the sample it carries, and whether a
  ! solve has filled it.
  !===================================================================!

  type :: gti_time_vertex

     type(gti_time_sample) :: sample
     logical               :: has_solution = .false.

  end type gti_time_vertex

  !===================================================================!
  ! One local motif application: which vertices, in what order,
  ! which of them is unknown, under which coefficient rows, at
  ! which time.
  !===================================================================!

  type :: gti_time_relation

     type(gti_time_motif) :: motif
     integer, allocatable :: sample_vertex(:)
     integer              :: unknown_sample = 0
     real(dp)             :: evaluation_time = 0.0_dp

   contains

     procedure :: arity          => relation_arity
     procedure :: unknown_vertex => relation_unknown_vertex

  end type gti_time_relation

  !===================================================================!
  ! The time graph: vertices and relations, and nothing else. The
  ! types keep their public singular names; Fortran denies a type
  ! its host module's name, so the module speaks in the plural.
  !===================================================================!

  type :: gti_time_graph

     type(gti_time_vertex)  , allocatable :: vertex(:)
     type(gti_time_relation), allocatable :: relation(:)

   contains

     procedure :: num_vertices
     procedure :: num_relations
     procedure :: validate
     procedure :: build_samples
     procedure :: write_unknown_q

  end type gti_time_graph

contains

  !===================================================================!
  ! How many sample vertices one relation touches.
  !===================================================================!

  pure function relation_arity(this) result(n)

    class(gti_time_relation), intent(in) :: this
    integer :: n

    if (allocated(this % sample_vertex)) then
       n = size(this % sample_vertex)
    else
       n = 0
    end if

  end function relation_arity

  !===================================================================!
  ! The unknown's GLOBAL vertex index: the local unknown position
  ! translated through the relation's vertex tuple. An invalid
  ! position answers 0 - absence, for the validator to refuse.
  !===================================================================!

  pure function relation_unknown_vertex(this) result(vertex_index)

    class(gti_time_relation), intent(in) :: this
    integer :: vertex_index

    vertex_index = 0
    if (this % unknown_sample < 1) return
    if (this % unknown_sample > this % arity()) return

    vertex_index = this % sample_vertex(this % unknown_sample)

  end function relation_unknown_vertex

  pure function num_vertices(this) result(n)

    class(gti_time_graph), intent(in) :: this
    integer :: n

    if (allocated(this % vertex)) then
       n = size(this % vertex)
    else
       n = 0
    end if

  end function num_vertices

  pure function num_relations(this) result(n)

    class(gti_time_graph), intent(in) :: this
    integer :: n

    if (allocated(this % relation)) then
       n = size(this % relation)
    else
       n = 0
    end if

  end function num_relations

  !===================================================================!
  ! The representation laws, refused in order. Validation proves
  ! the graph is well-formed; it evaluates nothing and solves
  ! nothing.
  !===================================================================!

  pure subroutine validate(this)

    class(gti_time_graph), intent(in) :: this

    integer :: r, k, m, arity, vertex_index

    if (this % num_vertices() < 1) then
       error stop 'gti_time_graph: at least one vertex is required'
    end if

    do r = 1, this % num_relations()

       arity = this % relation(r) % arity()
       if (arity < 1) then
          error stop 'gti_time_graph: relation has sample vertices'
       end if

       do k = 1, arity
          vertex_index = this % relation(r) % sample_vertex(k)
          if (vertex_index < 1 .or. vertex_index > this % num_vertices()) then
             error stop 'gti_time_graph: relation sample vertex is in range'
          end if
       end do

       if (this % relation(r) % unknown_sample < 1 .or. &
            & this % relation(r) % unknown_sample > arity) then
          error stop 'gti_time_graph: relation unknown sample is in range'
       end if

       if (this % relation(r) % motif % size() < 1) then
          error stop 'gti_time_graph: relation motif has rules'
       end if

       do m = 1, this % relation(r) % motif % size()
          if (.not. allocated(this % relation(r) % motif % rule(m) % weight)) then
             error stop 'gti_time_graph: motif rule has weights'
          end if
          if (size(this % relation(r) % motif % rule(m) % weight) /= arity) then
             error stop 'gti_time_graph: motif rule weight count matches relation arity'
          end if
       end do

       do k = 1, arity
          vertex_index = this % relation(r) % sample_vertex(k)
          if (.not. this % vertex(vertex_index) % sample % state % &
               & has_component(GTI_STATE_Q)) then
             error stop 'gti_time_graph: referenced vertex provides q'
          end if
       end do

    end do

  end subroutine validate

  !===================================================================!
  ! Materialize one relation's local sample set, in relation-local
  ! order, as deep copies: the local calculus consumes these
  ! without ever seeing the graph, and mutating them leaves the
  ! graph alone.
  !===================================================================!

  subroutine build_samples(this, relation_index, samples)

    class(gti_time_graph)             , intent(in)  :: this
    integer                           , intent(in)  :: relation_index
    type(gti_time_sample), allocatable, intent(out) :: samples(:)

    integer :: k, vertex_index, arity

    if (relation_index < 1 .or. relation_index > this % num_relations()) then
       error stop 'gti_time_graph: relation index is in range'
    end if

    arity = this % relation(relation_index) % arity()
    if (arity < 1) then
       error stop 'gti_time_graph: relation has sample vertices'
    end if

    allocate(samples(arity))

    do k = 1, arity
       vertex_index = this % relation(relation_index) % sample_vertex(k)
       if (vertex_index < 1 .or. vertex_index > this % num_vertices()) then
          error stop 'gti_time_graph: relation sample vertex is in range'
       end if
       samples(k) = this % vertex(vertex_index) % sample
    end do

  end subroutine build_samples

  !===================================================================!
  ! Write a solved q back to the shared vertex: values replaced,
  ! identity preserved - the injection mints from the vertex's own
  ! field - and the vertex marked solved. Sample time and non-q
  ! seats survive. Every relation sharing the vertex sees the new
  ! q on its next materialization; that shared seeing IS the time
  ! coupling.
  !===================================================================!

  subroutine write_unknown_q(this, relation_index, q)

    class(gti_time_graph) , intent(inout) :: this
    integer               , intent(in)    :: relation_index
    type(gti_value_buffer), intent(in)    :: q

    type(gti_time_local_unknown_problem) :: problem
    type(gti_time_sample)                :: one(1)
    type(gti_time_sample), allocatable   :: work(:)

    integer :: vertex_index

    if (relation_index < 1 .or. relation_index > this % num_relations()) then
       error stop 'gti_time_graph: relation index is in range'
    end if

    vertex_index = this % relation(relation_index) % unknown_vertex()
    if (vertex_index < 1) then
       error stop 'gti_time_graph: relation unknown sample is in range'
    end if
    if (vertex_index > this % num_vertices()) then
       error stop 'gti_time_graph: relation sample vertex is in range'
    end if

    one(1) = this % vertex(vertex_index) % sample

    ! a defined descriptor quiets a gfortran maybe-uninitialized
    ! false positive; the injection replaces it wholesale
    allocate(work(0))

    call problem % inject_trial_q(one, 1, q, work)

    this % vertex(vertex_index) % sample       = work(1)
    this % vertex(vertex_index) % has_solution = .true.

  end subroutine write_unknown_q

end module gti_time_graphs
