#include "scalar.fpp"

!=====================================================================!
! Continuous, pointwise, discretization-agnostic physics operators.
!
! A conservation law is written as
!
!     dq/dt + div F(q, grad q) = S(q, grad q)
!
! and a function of interest as  J = integral f(q, grad q) dV. The three
! operators - flux F, source S, objective f - are pointwise functions of
! the state q and its gradient grad q (first derivatives only; the
! divergence theorem / integration by parts removes the second). They
! provide their value and their state/design partials. They know nothing
! about meshes or cells: a finite-volume assembler integrates F over
! faces and S over the volume; a finite-element assembler would weight
! the same operators differently. This is the seam that keeps the
! framework law-agnostic (and discretization-agnostic).
!
! flux, source and objective extend a common base `physics`. Fortran has
! single inheritance, so a concrete law supplies a flux object and a
! source object separately (composed at setup), not one type that is
! both.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module interface_physics

  use iso_fortran_env, only : dp => REAL64
  use interface_graph, only : graph

  implicit none

  private
  public :: physics, source, objective, point_state

  !-------------------------------------------------------------------!
  ! Pointwise evaluation context: the state and its gradient at a point
  ! (a face centroid for the flux, a cell centroid for the source). The
  ! state is type(scalar) so complex-step partials drop in later; the
  ! spatial point x is geometry (real).
  !-------------------------------------------------------------------!

  type :: point_state
     integer                   :: nv = 1       ! number of variables
     type(scalar), allocatable :: q(:)         ! q        (nv)
     type(scalar), allocatable :: gradq(:,:)   ! grad q   (3, nv)
     real(dp)                  :: x(3) = 0.0_dp ! spatial point
  end type point_state

  !-------------------------------------------------------------------!
  ! Common base: a physics IS a graph - the coupling graph of the law.
  !
  !    vertices = the physical variables (the inherited num_vertices
  !               is the variable count - one fact, one home)
  !    edges    = "variable i's equation involves variable j"
  !
  !    scalar law:            (q1)              no edges
  !    two fields, decoupled: (q1)   (q2)       no edges - and that IS
  !                                             the decoupling, stated
  !                                             as structure
  !    a coupled law:         (q1)───(q2)       the edge is the reason
  !                                             the matrix will need
  !                                             cross-variable entries
  !-------------------------------------------------------------------!

  type, abstract, extends(graph) :: physics
   contains
     ! the graph contract, answered by the stored coupling adjacency
     procedure :: neighbours
     procedure :: degree
     ! become the coupling graph (called by every law's constructor)
     procedure :: form_coupling
     procedure :: num_design_vars => physics_num_design_vars  ! default 0
     procedure :: set_design_vars => physics_set_design_vars  ! default no-op
     procedure :: get_design_vars => physics_get_design_vars  ! default no-op
  end type physics

  !-------------------------------------------------------------------!
  ! Volumetric source S(q, grad q): value -> S(nv) + partials.
  !-------------------------------------------------------------------!

  type, extends(physics), abstract :: source
   contains
     procedure(source_value_interface), deferred :: value
     procedure :: dsource_dq      => source_dq_zero           ! (nv,nv)
     procedure :: dsource_dgradq  => source_dgradq_zero       ! (nv,3,nv)
     procedure :: dsource_ddesign => source_ddesign_zero      ! (nv) for design var k
  end type source

  !-------------------------------------------------------------------!
  ! Function-of-interest integrand f(q, grad q): value -> scalar + parts.
  !-------------------------------------------------------------------!

  type, extends(physics), abstract :: objective
   contains
     procedure(objective_value_interface), deferred :: value
     procedure :: dobj_dq      => objective_dq_zero           ! (nv)
     procedure :: dobj_dgradq  => objective_dgradq_zero       ! (3,nv)
     procedure :: dobj_ddesign => objective_ddesign_zero      ! scalar for design var k
  end type objective

  !-------------------------------------------------------------------!
  ! Deferred interfaces
  !-------------------------------------------------------------------!

  abstract interface

     pure function source_value_interface(this, st) result(S)
       import :: source, point_state
       class(source)    , intent(in) :: this
       type(point_state), intent(in) :: st
       type(scalar)                  :: S(this % num_vertices)
     end function source_value_interface

     pure function objective_value_interface(this, st) result(f)
       import :: objective, point_state
       class(objective) , intent(in) :: this
       type(point_state), intent(in) :: st
       type(scalar)                  :: f
     end function objective_value_interface

  end interface

contains

  !===================================================================!
  ! The graph contract: the coupling adjacency is stored by
  ! form_coupling below, so both queries are one-line delegations to
  ! the retained mechanism (the same shape the mesh uses).
  !===================================================================!

  pure function neighbours(this, v) result(nbrs)
    class(physics), intent(in) :: this
    integer       , intent(in) :: v
    integer, allocatable       :: nbrs(:)
    nbrs = this % stored_neighbours(v)
  end function neighbours

  pure integer function degree(this, v)
    class(physics), intent(in) :: this
    integer       , intent(in) :: v
    degree = this % stored_degree(v)
  end function degree

  !===================================================================!
  ! Become the coupling graph, in place: one vertex per variable, one
  ! edge per coupled pair the law declares. No edges means the
  ! variables are independent - which is every law we have today, and
  ! exactly why the assembled same-variable sparsity is right. A
  ! future coupled law hands its pairs here, and this graph becomes
  ! the instruction sheet for the cross-variable pattern.
  !===================================================================!

  pure subroutine form_coupling(this, n, tails, heads)

    class(physics), intent(inout)        :: this
    integer       , intent(in)           :: n
    integer       , intent(in), optional :: tails(:), heads(:)

    integer :: i

    this % num_vertices = n

    this % num_edges = 0
    if (present(tails)) this % num_edges = size(tails)
    if (allocated(this % edges)) deallocate(this % edges)
    allocate(this % edges(this % num_edges))
    do i = 1, this % num_edges
       this % edges(i) % tail = tails(i)
       this % edges(i) % head = heads(i)
    end do

    if (allocated(this % vertices)) deallocate(this % vertices)
    allocate(this % vertices(n))
    do i = 1, n
       this % vertices(i) % number = i
       this % vertices(i) % part   = 1
    end do

    call this % build_adjacency()

  end subroutine form_coupling

  !===================================================================!
  ! Design-variable defaults (no design dependence)
  !===================================================================!

  pure integer function physics_num_design_vars(this) result(n)
    class(physics), intent(in) :: this
    n = 0
  end function physics_num_design_vars

  pure subroutine physics_set_design_vars(this, x)
    class(physics), intent(inout) :: this
    real(dp)      , intent(in)    :: x(:)
  end subroutine physics_set_design_vars

  pure subroutine physics_get_design_vars(this, x)
    class(physics), intent(in)  :: this
    real(dp)      , intent(out) :: x(:)
    x = 0.0_dp
  end subroutine physics_get_design_vars

  !===================================================================!
  ! Default partials: zero (overridden by laws that depend on them)
  !===================================================================!

  pure function source_dq_zero(this, st) result(dS)
    class(source)    , intent(in) :: this
    type(point_state), intent(in) :: st
    type(scalar)                  :: dS(this % num_vertices, this % num_vertices)
    dS = 0.0_dp
  end function source_dq_zero

  pure function source_dgradq_zero(this, st) result(dS)
    class(source)    , intent(in) :: this
    type(point_state), intent(in) :: st
    type(scalar)                  :: dS(this % num_vertices, 3, this % num_vertices)
    dS = 0.0_dp
  end function source_dgradq_zero

  pure function source_ddesign_zero(this, st, k) result(dS)
    class(source)    , intent(in) :: this
    type(point_state), intent(in) :: st
    integer          , intent(in) :: k
    type(scalar)                  :: dS(this % num_vertices)
    dS = 0.0_dp
  end function source_ddesign_zero

  pure function objective_dq_zero(this, st) result(df)
    class(objective) , intent(in) :: this
    type(point_state), intent(in) :: st
    type(scalar)                  :: df(this % num_vertices)
    df = 0.0_dp
  end function objective_dq_zero

  pure function objective_dgradq_zero(this, st) result(df)
    class(objective) , intent(in) :: this
    type(point_state), intent(in) :: st
    type(scalar)                  :: df(3, this % num_vertices)
    df = 0.0_dp
  end function objective_dgradq_zero

  pure function objective_ddesign_zero(this, st, k) result(df)
    class(objective) , intent(in) :: this
    type(point_state), intent(in) :: st
    integer          , intent(in) :: k
    type(scalar)                  :: df
    df = 0.0_dp
  end function objective_ddesign_zero

end module interface_physics
