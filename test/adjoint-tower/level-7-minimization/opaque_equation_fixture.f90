!=====================================================================!
! THE OPAQUE EQUATION FIXTURE - test-local, deliberately outside
! src/. Level 7 is allowed to SUPPLY its numerical equations rather
! than generate them, exactly as every earlier tower supplied a
! residual before constituting one. Two equations stand here, and
! they are deliberately INDEPENDENT of each other:
!
!      primal    R(q, p=2) = [ 2u +  v -  8,
!                              3u + 4v - 22 ]      Q -> Y
!
!      adjoint   A^T lambda - c = [ 2L1 + 3L2 - 1,
!                                   1L1 + 4L2 - 2 ] Y -> Q
!
! That duplication is the rung's point and its temporary sin: Level
! 7 tests whether ONE minimizer family can govern equations whose
! unknown and residual domains are exchanged. It does not yet test
! whether the two equations come from one truth - Level 8 removes
! the independence and generates both from a single constitution.
!
! Note the orientations. The primal takes a state on Q and answers
! on Y; the adjoint takes a covector on Y and answers on Q. Both
! domains have two members, so nothing but IDENTITY separates them,
! and each operation checks its input with same_as before reading a
! single number. A is non-symmetric on purpose: solving the primal
! orientation where the transpose belongs would return [0.4, 0.2]
! instead of [-0.4, 0.6], and no dimension check would notice.
!
! Every read and write goes through the domain's own local_index -
! never a raw member id, never an assumed position.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module opaque_equation_fixture

  use iso_fortran_env  , only : dp => REAL64
  use adjoint_assert   , only : VAR_U, VAR_V, TGT_R1, TGT_R2
  use fractal_graph        , only : set_graph => graph
  use graph_set_map        , only : set_map
  use graph_set_representation, only : set_representation
  use graph_operation_view, only : graph_operation
  use graph_directed_view, only : directed_graph
  use graph_field_calculus, only : graph_field
  use class_graph_field, only : field

  implicit none

  private
  public :: opaque_primal, opaque_adjoint

  !===================================================================!
  ! The primal equation: a state on Q, a residual on Y.
  !===================================================================!

  type, extends(graph_operation) :: opaque_primal
     ! Identity and count. Both domains are counted 1..n here,
     ! so no coordinates beyond the count are needed.
     type(set_graph) :: q_dom, y_dom
     integer         :: n_q_dom = 0, n_y_dom = 0
     ! Both domains are LISTED, so a member is not its position: the
     ! action keeps its own coordinates, copied at construction.
     class(set_representation), allocatable :: c_q, c_y
   contains
     procedure :: name   => primal_name
     procedure :: domain => primal_domain
     procedure :: apply  => primal_apply
  end type opaque_primal

  !===================================================================!
  ! The adjoint equation: a covector on Y, a residual on Q. The
  ! orientation is exchanged, and nothing else is.
  !===================================================================!

  type, extends(graph_operation) :: opaque_adjoint
     ! Identity and count. Both domains are counted 1..n here,
     ! so no coordinates beyond the count are needed.
     type(set_graph) :: q_dom, y_dom
     integer         :: n_q_dom = 0, n_y_dom = 0
     ! Both domains are LISTED, so a member is not its position: the
     ! action keeps its own coordinates, copied at construction.
     class(set_representation), allocatable :: c_q, c_y
   contains
     procedure :: name   => adjoint_name
     procedure :: domain => adjoint_domain
     procedure :: apply  => adjoint_apply
  end type opaque_adjoint

  interface opaque_primal
     module procedure create_primal
  end interface opaque_primal

  interface opaque_adjoint
     module procedure create_adjoint
  end interface opaque_adjoint

contains

  type(opaque_primal) function create_primal(q_dom, y_dom, sets) result(this)
    type(set_graph), intent(in) :: q_dom, y_dom
    type(set_map)  , intent(in) :: sets
    this % q_dom = q_dom ; this % n_q_dom = sets % size_of(q_dom)
    this % y_dom = y_dom ; this % n_y_dom = sets % size_of(y_dom)
    call sets % extent_of(q_dom, this % c_q)
    call sets % extent_of(y_dom, this % c_y)
  end function create_primal

  type(opaque_adjoint) function create_adjoint(y_dom, q_dom, sets) result(this)
    type(set_graph), intent(in) :: y_dom, q_dom
    type(set_map)  , intent(in) :: sets
    this % y_dom = y_dom ; this % n_y_dom = sets % size_of(y_dom)
    this % q_dom = q_dom ; this % n_q_dom = sets % size_of(q_dom)
    call sets % extent_of(y_dom, this % c_y)
    call sets % extent_of(q_dom, this % c_q)
  end function create_adjoint

  pure function primal_name(this) result(name)
    class(opaque_primal), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'opaque primal equation'
  end function primal_name

  pure function adjoint_name(this) result(name)
    class(opaque_adjoint), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'opaque adjoint equation'
  end function adjoint_name

  subroutine primal_domain(this, input_graph, domain, nentries)
    class(opaque_primal), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    type(set_graph), intent(out) :: domain
    integer        , intent(out) :: nentries
    associate (u1 => input_graph); end associate
    domain   = this % y_dom
    nentries = this % n_y_dom     ! answers on Y
  end subroutine primal_domain

  subroutine adjoint_domain(this, input_graph, domain, nentries)
    class(opaque_adjoint), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    type(set_graph), intent(out) :: domain
    integer        , intent(out) :: nentries
    associate (u1 => input_graph); end associate
    domain   = this % q_dom
    nentries = this % n_q_dom     ! answers on Q
  end subroutine adjoint_domain

  !===================================================================!
  ! R(q, p=2) on Y, read from a state on Q. The host graph is
  ! ignored: this equation is a function of its input field alone.
  !===================================================================!

  subroutine primal_apply(this, input_graph, input_data, output)

    class(opaque_primal), intent(in)               :: this
    class(directed_graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(field)                    :: out
    type(set_graph) :: dom
    real(dp), allocatable          :: q(:), r(:)
    real(dp)                       :: u, v
    type(set_map)     :: sets

    associate (u1 => input_graph); end associate

    if (.not. present(input_data)) then
       error stop 'opaque primal: the equation needs a state to judge'
    end if
    dom = input_data(1) % domain()
    if (.not. dom % same_as(this % q_dom)) then
       error stop 'opaque primal: the state must live on the state domain'
    end if
    call input_data(1) % get_real_vector(q)

    u = q(this % c_q % local_index(VAR_U))
    v = q(this % c_q % local_index(VAR_V))

    allocate(r(this % n_y_dom))
    r(this % c_y % local_index(TGT_R1)) = 2.0_dp*u +      v -  8.0_dp
    r(this % c_y % local_index(TGT_R2)) = 3.0_dp*u + 4.0_dp*v - 22.0_dp

    out = field('residual', this % y_dom, this % n_y_dom)
    call out % set_real_vector(r)
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine primal_apply

  !===================================================================!
  ! A^T lambda - c on Q, read from a covector on Y. The rows of this
  ! residual are indexed by STATE slots, because that is what the
  ! transposed operator answers on.
  !===================================================================!

  subroutine adjoint_apply(this, input_graph, input_data, output)

    class(opaque_adjoint), intent(in)              :: this
    class(directed_graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(field)                    :: out
    type(set_graph) :: dom
    real(dp), allocatable          :: lam(:), r(:)
    real(dp)                       :: l1, l2
    type(set_map)     :: sets

    associate (u1 => input_graph); end associate

    if (.not. present(input_data)) then
       error stop 'opaque adjoint: the equation needs a covector to judge'
    end if
    dom = input_data(1) % domain()
    if (.not. dom % same_as(this % y_dom)) then
       error stop 'opaque adjoint: the covector must live on the residual-row domain'
    end if
    call input_data(1) % get_real_vector(lam)

    l1 = lam(this % c_y % local_index(TGT_R1))
    l2 = lam(this % c_y % local_index(TGT_R2))

    allocate(r(this % n_q_dom))
    r(this % c_q % local_index(VAR_U)) = 2.0_dp*l1 + 3.0_dp*l2 - 1.0_dp
    r(this % c_q % local_index(VAR_V)) =      l1 + 4.0_dp*l2 - 2.0_dp

    out = field('adjoint residual', this % q_dom, this % n_q_dom)
    call out % set_real_vector(r)
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine adjoint_apply

  !===================================================================!
  ! Both domains here enumerate 1..n, so a member IS its position and
  ! the action needs no stored representation to say where it stands.
  !===================================================================!


end module opaque_equation_fixture
