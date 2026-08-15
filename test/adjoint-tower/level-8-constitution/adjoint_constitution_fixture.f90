!=====================================================================!
! THE ADJOINT CONSTITUTION FIXTURE - test-local, deliberately
! outside src/: one 2x2 specimen has not earned a production
! derivative API. This module is the ONE place in the tower where
! the law data live:
!
!      R(q,p) = A q - b p        f(q,p) = c^T q + d p
!
!      A: (r1,u)=2  (r1,v)=1     b: (r1,p)=4   c: (f,u)=1
!         (r2,u)=3  (r2,v)=4        (r2,p)=11     (f,v)=2
!                                                d: (f,p)=2
!
! and every coefficient is stored EXACTLY ONCE, keyed by the pair
! of MEMBERS it belongs to - never by position, never by a second
! transposed table. The forward state action and the reverse state
! action read the same coeff_state; the response's forward and
! reverse actions read the same coeff_response_state. There is no
! A^T anywhere in this file.
!
! Nothing here loops over "two rows and two columns". Every action
! walks the STRUCTURAL SUPPORT derived in Gate A - J_Q, J_P, F_Q,
! F_P - asking the relation which incidences exist and asking the
! law what each one is worth:
!
!      forward:   out(row) +=  coeff(row,col) * in(col)
!      reverse:   out(col) +=  coeff(row,col) * in(row)
!
! one lookup, two directions, the reverse accumulating with += as
! it must. Because the support answers by member and every vector
! is indexed through its domain's own local_index, the whole file
! is blind to how any role happens to be enumerated - which is what
! the permutation test at Level 8 exists to prove.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module adjoint_constitution_fixture

  use iso_fortran_env  , only : dp => REAL64
  use adjoint_assert   , only : VAR_P, VAR_U, VAR_V
  use adjoint_assert   , only : TGT_R1, TGT_R2, TGT_F
  use graph_set    , only : set
  use graph_relation   , only : relation
  use graph_grammar    , only : ordinary_graph, graph_field, graph_operation
  use class_graph_field, only : field

  implicit none

  private
  public :: coeff_state, coeff_param, coeff_response_state
  public :: coeff_response_param
  public :: residual_of, response_of
  public :: rq_forward, rq_reverse, rp_forward
  public :: fq_forward, fq_reverse, fp_forward
  public :: constituted_primal, constituted_adjoint, constituted_tangent

  !===================================================================!
  ! The three operation faces the solvers consume. Each holds the
  ! supports and domains it needs and NO coefficients of its own -
  ! every number it uses comes from the law table below.
  !===================================================================!

  type, extends(graph_operation) :: constituted_primal
     class(relation)  , allocatable :: jq, jp
     class(set), allocatable :: q_dom, y_dom, p_dom
     real(dp)         , allocatable :: p_val(:)
   contains
     procedure :: name   => primal_name
     procedure :: domain => primal_domain
     procedure :: apply  => primal_apply
  end type constituted_primal

  type, extends(graph_operation) :: constituted_adjoint
     class(relation)  , allocatable :: jq, fq
     class(set), allocatable :: q_dom, y_dom, z_dom
   contains
     procedure :: name   => adjoint_name
     procedure :: domain => adjoint_domain
     procedure :: apply  => adjoint_apply
  end type constituted_adjoint

  type, extends(graph_operation) :: constituted_tangent
     class(relation)  , allocatable :: jq, jp
     class(set), allocatable :: q_dom, y_dom, p_dom
     real(dp)         , allocatable :: dp_val(:)
   contains
     procedure :: name   => tangent_name
     procedure :: domain => tangent_domain
     procedure :: apply  => tangent_apply
  end type constituted_tangent

  interface constituted_primal
     module procedure create_primal
  end interface constituted_primal

  interface constituted_adjoint
     module procedure create_adjoint
  end interface constituted_adjoint

  interface constituted_tangent
     module procedure create_tangent
  end interface constituted_tangent

contains

  !===================================================================!
  ! THE LAW TABLE. Every coefficient appears once, keyed by the
  ! members it relates. An incidence no law binds refuses.
  !===================================================================!

  real(dp) function coeff_state(row, col) result(a)

    integer, intent(in) :: row, col      ! row in Y, col in Q

    if      (row .eq. TGT_R1 .and. col .eq. VAR_U) then; a = 2.0_dp
    else if (row .eq. TGT_R1 .and. col .eq. VAR_V) then; a = 1.0_dp
    else if (row .eq. TGT_R2 .and. col .eq. VAR_U) then; a = 3.0_dp
    else if (row .eq. TGT_R2 .and. col .eq. VAR_V) then; a = 4.0_dp
    else
       error stop 'constitution: no law binds this state incidence'
    end if

  end function coeff_state

  !===================================================================!
  ! dR/dp = -b: the residual falls as the parameter rises.
  !===================================================================!

  real(dp) function coeff_param(row, col) result(a)

    integer, intent(in) :: row, col      ! row in Y, col in P

    if      (row .eq. TGT_R1 .and. col .eq. VAR_P) then; a = -4.0_dp
    else if (row .eq. TGT_R2 .and. col .eq. VAR_P) then; a = -11.0_dp
    else
       error stop 'constitution: no law binds this parameter incidence'
    end if

  end function coeff_param

  real(dp) function coeff_response_state(row, col) result(a)

    integer, intent(in) :: row, col      ! row in Z, col in Q

    if      (row .eq. TGT_F .and. col .eq. VAR_U) then; a = 1.0_dp
    else if (row .eq. TGT_F .and. col .eq. VAR_V) then; a = 2.0_dp
    else
       error stop 'constitution: no law binds this response-state incidence'
    end if

  end function coeff_response_state

  real(dp) function coeff_response_param(row, col) result(a)

    integer, intent(in) :: row, col      ! row in Z, col in P

    if (row .eq. TGT_F .and. col .eq. VAR_P) then
       a = 2.0_dp
    else
       error stop 'constitution: no law binds this response-parameter incidence'
    end if

  end function coeff_response_param

  !===================================================================!
  ! The forward state action: Rq v, walking J_Q.
  !===================================================================!

  subroutine rq_forward(jq, y_dom, q_dom, v_q, out_y)

    class(relation)  , intent(in)  :: jq
    class(set), intent(in)  :: y_dom, q_dom
    real(dp)         , intent(in)  :: v_q(:)
    real(dp)         , intent(out) :: out_y(:)

    integer :: i, j, row, col

    if (size(v_q) .ne. q_dom % size() .or. &
         & size(out_y) .ne. y_dom % size()) then
       error stop 'constitution: every vector is sized by its domain'
    end if

    out_y = 0.0_dp
    do i = 1, y_dom % size()
       row = y_dom % member(i)
       do j = 1, q_dom % size()
          col = q_dom % member(j)
          if (jq % has([row, col])) then
             out_y(y_dom % local_index(row)) = &
                  & out_y(y_dom % local_index(row)) + &
                  & coeff_state(row, col) * v_q(q_dom % local_index(col))
          end if
       end do
    end do

  end subroutine rq_forward

  !===================================================================!
  ! The reverse state action: Rq^T lambda, walking the SAME J_Q with
  ! the SAME coefficients, accumulating into the state slots. No
  ! transposed table exists to disagree with the forward one.
  !===================================================================!

  subroutine rq_reverse(jq, y_dom, q_dom, bar_y, out_q)

    class(relation)  , intent(in)  :: jq
    class(set), intent(in)  :: y_dom, q_dom
    real(dp)         , intent(in)  :: bar_y(:)
    real(dp)         , intent(out) :: out_q(:)

    integer :: i, j, row, col

    if (size(bar_y) .ne. y_dom % size() .or. &
         & size(out_q) .ne. q_dom % size()) then
       error stop 'constitution: every vector is sized by its domain'
    end if

    out_q = 0.0_dp
    do i = 1, y_dom % size()
       row = y_dom % member(i)
       do j = 1, q_dom % size()
          col = q_dom % member(j)
          if (jq % has([row, col])) then
             out_q(q_dom % local_index(col)) = &
                  & out_q(q_dom % local_index(col)) + &
                  & coeff_state(row, col) * bar_y(y_dom % local_index(row))
          end if
       end do
    end do

  end subroutine rq_reverse

  !===================================================================!
  ! The parameter action: Rp dp, walking J_P.
  !===================================================================!

  subroutine rp_forward(jp, y_dom, p_dom, dp_p, out_y)

    class(relation)  , intent(in)  :: jp
    class(set), intent(in)  :: y_dom, p_dom
    real(dp)         , intent(in)  :: dp_p(:)
    real(dp)         , intent(out) :: out_y(:)

    integer :: i, j, row, col

    out_y = 0.0_dp
    do i = 1, y_dom % size()
       row = y_dom % member(i)
       do j = 1, p_dom % size()
          col = p_dom % member(j)
          if (jp % has([row, col])) then
             out_y(y_dom % local_index(row)) = &
                  & out_y(y_dom % local_index(row)) + &
                  & coeff_param(row, col) * dp_p(p_dom % local_index(col))
          end if
       end do
    end do

  end subroutine rp_forward

  !===================================================================!
  ! The response's state action and its reverse - again one lookup,
  ! two directions. The reverse with a unit seed is how f_q^T is
  ! obtained without ever writing f_q^T down.
  !===================================================================!

  subroutine fq_forward(fq, z_dom, q_dom, v_q, out_z)

    class(relation)  , intent(in)  :: fq
    class(set), intent(in)  :: z_dom, q_dom
    real(dp)         , intent(in)  :: v_q(:)
    real(dp)         , intent(out) :: out_z(:)

    integer :: i, j, row, col

    out_z = 0.0_dp
    do i = 1, z_dom % size()
       row = z_dom % member(i)
       do j = 1, q_dom % size()
          col = q_dom % member(j)
          if (fq % has([row, col])) then
             out_z(z_dom % local_index(row)) = &
                  & out_z(z_dom % local_index(row)) + &
                  & coeff_response_state(row, col) * &
                  & v_q(q_dom % local_index(col))
          end if
       end do
    end do

  end subroutine fq_forward

  subroutine fq_reverse(fq, z_dom, q_dom, bar_z, out_q)

    class(relation)  , intent(in)  :: fq
    class(set), intent(in)  :: z_dom, q_dom
    real(dp)         , intent(in)  :: bar_z(:)
    real(dp)         , intent(out) :: out_q(:)

    integer :: i, j, row, col

    out_q = 0.0_dp
    do i = 1, z_dom % size()
       row = z_dom % member(i)
       do j = 1, q_dom % size()
          col = q_dom % member(j)
          if (fq % has([row, col])) then
             out_q(q_dom % local_index(col)) = &
                  & out_q(q_dom % local_index(col)) + &
                  & coeff_response_state(row, col) * &
                  & bar_z(z_dom % local_index(row))
          end if
       end do
    end do

  end subroutine fq_reverse

  subroutine fp_forward(fp, z_dom, p_dom, dp_p, out_z)

    class(relation)  , intent(in)  :: fp
    class(set), intent(in)  :: z_dom, p_dom
    real(dp)         , intent(in)  :: dp_p(:)
    real(dp)         , intent(out) :: out_z(:)

    integer :: i, j, row, col

    out_z = 0.0_dp
    do i = 1, z_dom % size()
       row = z_dom % member(i)
       do j = 1, p_dom % size()
          col = p_dom % member(j)
          if (fp % has([row, col])) then
             out_z(z_dom % local_index(row)) = &
                  & out_z(z_dom % local_index(row)) + &
                  & coeff_response_param(row, col) * &
                  & dp_p(p_dom % local_index(col))
          end if
       end do
    end do

  end subroutine fp_forward

  !===================================================================!
  ! The nonlinear-looking faces: R(q,p) = Aq - bp and
  ! f(q,p) = c^T q + dp, both assembled from the same actions.
  !===================================================================!

  subroutine residual_of(jq, jp, y_dom, q_dom, p_dom, q_val, p_val, &
       & out_y)

    class(relation)  , intent(in)  :: jq, jp
    class(set), intent(in)  :: y_dom, q_dom, p_dom
    real(dp)         , intent(in)  :: q_val(:), p_val(:)
    real(dp)         , intent(out) :: out_y(:)

    real(dp) :: from_state(y_dom % size()), from_param(y_dom % size())

    call rq_forward(jq, y_dom, q_dom, q_val, from_state)
    call rp_forward(jp, y_dom, p_dom, p_val, from_param)
    out_y = from_state + from_param

  end subroutine residual_of

  subroutine response_of(fq, fp, z_dom, q_dom, p_dom, q_val, p_val, &
       & out_z)

    class(relation)  , intent(in)  :: fq, fp
    class(set), intent(in)  :: z_dom, q_dom, p_dom
    real(dp)         , intent(in)  :: q_val(:), p_val(:)
    real(dp)         , intent(out) :: out_z(:)

    real(dp) :: from_state(z_dom % size()), from_param(z_dom % size())

    call fq_forward(fq, z_dom, q_dom, q_val, from_state)
    call fp_forward(fp, z_dom, p_dom, p_val, from_param)
    out_z = from_state + from_param

  end subroutine response_of

  !===================================================================!
  ! The three operation faces. Each is pure delegation to the
  ! actions above; none holds a coefficient of its own.
  !===================================================================!

  type(constituted_primal) function create_primal(jq, jp, q_dom, &
       & y_dom, p_dom, p_val) result(this)
    class(relation)  , intent(in) :: jq, jp
    class(set), intent(in) :: q_dom, y_dom, p_dom
    real(dp)         , intent(in) :: p_val(:)
    allocate(this % jq, source=jq)
    allocate(this % jp, source=jp)
    allocate(this % q_dom, source=q_dom)
    allocate(this % y_dom, source=y_dom)
    allocate(this % p_dom, source=p_dom)
    this % p_val = p_val
  end function create_primal

  type(constituted_adjoint) function create_adjoint(jq, fq, y_dom, &
       & q_dom, z_dom) result(this)
    class(relation)  , intent(in) :: jq, fq
    class(set), intent(in) :: y_dom, q_dom, z_dom
    allocate(this % jq, source=jq)
    allocate(this % fq, source=fq)
    allocate(this % y_dom, source=y_dom)
    allocate(this % q_dom, source=q_dom)
    allocate(this % z_dom, source=z_dom)
  end function create_adjoint

  type(constituted_tangent) function create_tangent(jq, jp, q_dom, &
       & y_dom, p_dom, dp_val) result(this)
    class(relation)  , intent(in) :: jq, jp
    class(set), intent(in) :: q_dom, y_dom, p_dom
    real(dp)         , intent(in) :: dp_val(:)
    allocate(this % jq, source=jq)
    allocate(this % jp, source=jp)
    allocate(this % q_dom, source=q_dom)
    allocate(this % y_dom, source=y_dom)
    allocate(this % p_dom, source=p_dom)
    this % dp_val = dp_val
  end function create_tangent

  pure function primal_name(this) result(name)
    class(constituted_primal), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'constituted residual'
  end function primal_name

  pure function adjoint_name(this) result(name)
    class(constituted_adjoint), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'constituted adjoint equation'
  end function adjoint_name

  pure function tangent_name(this) result(name)
    class(constituted_tangent), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'constituted tangent equation'
  end function tangent_name

  subroutine primal_domain(this, input_graph, domain)
    class(constituted_primal), intent(in) :: this
    class(ordinary_graph), intent(in) :: input_graph
    class(set), allocatable, intent(out) :: domain
    associate (u1 => input_graph); end associate
    allocate(domain, source=this % y_dom)
  end subroutine primal_domain

  subroutine adjoint_domain(this, input_graph, domain)
    class(constituted_adjoint), intent(in) :: this
    class(ordinary_graph), intent(in) :: input_graph
    class(set), allocatable, intent(out) :: domain
    associate (u1 => input_graph); end associate
    allocate(domain, source=this % q_dom)
  end subroutine adjoint_domain

  subroutine tangent_domain(this, input_graph, domain)
    class(constituted_tangent), intent(in) :: this
    class(ordinary_graph), intent(in) :: input_graph
    class(set), allocatable, intent(out) :: domain
    associate (u1 => input_graph); end associate
    allocate(domain, source=this % y_dom)
  end subroutine tangent_domain

  !===================================================================!
  ! R(q, p) - a state on Q in, a residual on Y out.
  !===================================================================!

  subroutine primal_apply(this, input_graph, input_data, output)

    class(constituted_primal), intent(in)          :: this
    class(ordinary_graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(field)                    :: out
    class(set), allocatable :: dom
    real(dp), allocatable          :: q(:), r(:)

    associate (u1 => input_graph); end associate

    if (.not. present(input_data)) then
       error stop 'constitution: the residual needs a state to judge'
    end if
    call input_data(1) % domain(dom)
    if (.not. dom % equals(this % q_dom)) then
       error stop 'constitution: the state must live on the state domain'
    end if
    call input_data(1) % get_real_vector(q)

    allocate(r(this % y_dom % size()))
    call residual_of(this % jq, this % jp, this % y_dom, this % q_dom, &
         & this % p_dom, q, this % p_val, r)

    out = field('residual', this % y_dom)
    call out % set_real_vector(r)
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine primal_apply

  !===================================================================!
  ! Rq^T lambda - fq^T - a covector on Y in, a residual on Q out.
  ! Both terms are reverse actions of the ONE law table: the state
  ! block reversed, and the response block reversed with a unit
  ! seed. No transposed coefficients are written anywhere.
  !===================================================================!

  subroutine adjoint_apply(this, input_graph, input_data, output)

    class(constituted_adjoint), intent(in)         :: this
    class(ordinary_graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(field)                    :: out
    class(set), allocatable :: dom
    real(dp), allocatable          :: lam(:), r(:), rhs(:), seed(:)

    associate (u1 => input_graph); end associate

    if (.not. present(input_data)) then
       error stop 'constitution: the adjoint equation needs a covector to judge'
    end if
    call input_data(1) % domain(dom)
    if (.not. dom % equals(this % y_dom)) then
       error stop 'constitution: the covector must live on the residual-row domain'
    end if
    call input_data(1) % get_real_vector(lam)

    allocate(r(this % q_dom % size()), rhs(this % q_dom % size()))
    allocate(seed(this % z_dom % size()))

    call rq_reverse(this % jq, this % y_dom, this % q_dom, lam, r)

    seed = 1.0_dp
    call fq_reverse(this % fq, this % z_dom, this % q_dom, seed, rhs)

    out = field('adjoint residual', this % q_dom)
    call out % set_real_vector(r - rhs)
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine adjoint_apply

  !===================================================================!
  ! Rq q_p + Rp dp - a state direction on Q in, a residual on Y out.
  ! The forward reading of the same law.
  !===================================================================!

  subroutine tangent_apply(this, input_graph, input_data, output)

    class(constituted_tangent), intent(in)         :: this
    class(ordinary_graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(field)                    :: out
    class(set), allocatable :: dom
    real(dp), allocatable          :: qp(:), r(:), from_param(:)

    associate (u1 => input_graph); end associate

    if (.not. present(input_data)) then
       error stop 'constitution: the tangent equation needs a direction to judge'
    end if
    call input_data(1) % domain(dom)
    if (.not. dom % equals(this % q_dom)) then
       error stop 'constitution: the state must live on the state domain'
    end if
    call input_data(1) % get_real_vector(qp)

    allocate(r(this % y_dom % size()), from_param(this % y_dom % size()))

    call rq_forward(this % jq, this % y_dom, this % q_dom, qp, r)
    call rp_forward(this % jp, this % y_dom, this % p_dom, &
         & this % dp_val, from_param)

    out = field('tangent residual', this % y_dom)
    call out % set_real_vector(r + from_param)
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine tangent_apply

end module adjoint_constitution_fixture
