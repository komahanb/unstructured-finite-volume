!=====================================================================!
! LEARNING TOWER . LEVEL 5 . FIELD CALCULUS
!
! The level answers one question: WHAT NUMERICAL VALUES LIVE ON
! WHICH DOMAINS. The first numbers enter - and each has exactly one
! home:
!
!      K     = { y, x }     observed data        [6, 2]
!      Theta = { w }        trainable parameter  w0 = 0
!      U     = { e, yhat }  computed later       NO FIELD AT ALL
!
! Data and parameters differ by DOMAIN and role, not by field class:
! no data_field, no parameter_field, no tensor - one field type,
! three subdomains, pairwise disjoint and covering V exactly once.
! The zero on Theta means one thing only: the initial parameter.
! U carries no field, no zero, no NaN, no sentinel - uncomputed
! slots have no fabricated values, and NOT constructing a field is
! a stronger truth than constructing an empty one.
!
! And the import list is the rung's own proof: a field needs a
! domain - learning_assert, graph_carrier, class_graph_field - and
! no graph topology anywhere. The model law is still unspoken:
! knowing x = 2 is not knowing what predict does with it.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program learning_level_5

  use iso_fortran_env  , only : dp => REAL64
  use learning_assert  , only : report, verdict
  use learning_assert  , only : SLOT_W, SLOT_X, SLOT_YHAT, SLOT_Y, SLOT_E
  use fractal_graph        , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use graph_set_map        , only : set_map
  use graph_inclusion_map  , only : inclusion_map, declared_subobject
  use class_graph_field, only : field

  implicit none

  type(set_graph) :: v
  type(set_graph)  :: k, theta, u
  type(field)       :: qk, theta0
  integer           :: nfail
  type(set_map)     :: sets
  type(inclusion_map)     :: inclusions

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "learning tower . level 5 . field calculus"
  write(*,'(1x,a)') "============================================="

  call v % declare()
  call sets % bind(v, counted_set_representation(5))

  ! The three roles, declared in deliberately nonnumeric order.
  call k % declare()
  call sets       % bind(k, listed_set_representation([SLOT_Y, SLOT_X]))
  call inclusions % include_in(k, v)
  call theta % declare()
  call sets       % bind(theta, listed_set_representation([SLOT_W]))
  call inclusions % include_in(theta, v)
  call u % declare()
  call sets       % bind(u, listed_set_representation([SLOT_E, SLOT_YHAT]))
  call inclusions % include_in(u, v)

  call check_partition(nfail)
  call check_observed_field(nfail)
  call check_parameter_field(nfail)
  call check_domains_distinguish(nfail)

  ! U deliberately receives NO field here: no q_U, no zeros, no
  ! sentinels. The subdomain alone says "computed later", and that
  ! absence is the rung's strongest statement.

  call verdict(nfail, "level 5")

contains

  !===================================================================!
  ! Extensions, embeddings, and the one-home law: every member of V
  ! belongs to exactly one of K, Theta, U - disjointness and
  ! coverage proved together, composed locally from membership.
  !===================================================================!

  subroutine check_partition(nfail)

    integer, intent(inout) :: nfail

    integer :: i, m, homes
    logical :: ok

    call report(sets % size_of(k) .eq. 2 .and. sets % has_in(k, SLOT_Y) .and. &
         &      sets % has_in(k, SLOT_X) .and. .not. sets % has_in(k, SLOT_W), &
         & "K = { y, x }, the observations", nfail)
    call report(sets % size_of(theta) .eq. 1 .and. sets % has_in(theta, SLOT_W) .and. &
         &      .not. sets % has_in(theta, SLOT_X), &
         & "Theta = { w }, the one trainable", nfail)
    call report(sets % size_of(u) .eq. 2 .and. sets % has_in(u, SLOT_E) .and. &
         &      sets % has_in(u, SLOT_YHAT) .and. .not. sets % has_in(u, SLOT_Y), &
         & "U = { e, yhat }, where answers will one day live", nfail)

    call report(declared_subobject(k, v, inclusions) .and. &
         &      declared_subobject(theta, v, inclusions) .and. &
         &      declared_subobject(u, v, inclusions), &
         & "all three stand embedded in the value slots", nfail)

    ok = .true.
    do i = 1, sets % size_of(v)
       m = sets % member_of(v, i)
       homes = count([sets % has_in(k, m), sets % has_in(theta, m), sets % has_in(u, m)])
       ok = ok .and. (homes .eq. 1)
    end do
    call report(ok, &
         & "every slot has exactly one home: disjoint, and covering", &
         & nfail)

  end subroutine check_partition

  !===================================================================!
  ! The observed data: [6, 2] on K = { y, x } - storage follows the
  ! DECLARATION, y's value first, and every read walks local_index.
  !===================================================================!

  subroutine check_observed_field(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: dom
    real(dp), allocatable          :: val(:)

    qk = field('observations', k, sets % size_of(k))
    call qk % set_real_vector([6.0_dp, 2.0_dp])

    dom = qk % domain()
    call report(dom % same_as(k), &
         & "the data field's domain is K, by identity", nfail)
    call report(qk % num_entries() .eq. 2 .and. &
         &      qk % num_components() .eq. 1, &
         & "two entries, one component", nfail)

    call qk % get_real_vector(val)
    call report(abs(val(sets % index_in(k, SLOT_Y)) - 6.0_dp) < 1.0d-14, &
         & "y = 6, read through the domain map", nfail)
    call report(abs(val(sets % index_in(k, SLOT_X)) - 2.0_dp) < 1.0d-14, &
         & "x = 2", nfail)
    call report(abs(val(1) - 6.0_dp) < 1.0d-14 .and. &
         &      abs(val(2) - 2.0_dp) < 1.0d-14, &
         & "storage follows declaration: y's value stands first", nfail)

  end subroutine check_observed_field

  !===================================================================!
  ! The trainable parameter: w0 = 0, and that zero means INITIAL -
  ! nothing else, nowhere else.
  !===================================================================!

  subroutine check_parameter_field(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: dom
    real(dp), allocatable          :: val(:)

    theta0 = field('parameters', theta, sets % size_of(theta))
    call theta0 % set_real_vector([0.0_dp])

    dom = theta0 % domain()
    call report(dom % same_as(theta), &
         & "the parameter field's domain is Theta, by identity", nfail)
    call report(theta0 % num_entries() .eq. 1, &
         & "one trainable entry", nfail)

    call theta0 % get_real_vector(val)
    call report(abs(val(sets % index_in(theta, SLOT_W))) < 1.0d-14, &
         & "w = 0: training has not happened, and nothing pretends it has", &
         & nfail)

  end subroutine check_parameter_field

  !===================================================================!
  ! Data and parameters are told apart by their DOMAINS - no
  ! subtype, no flag, no tensor: the identities do all the work.
  !===================================================================!

  subroutine check_domains_distinguish(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: dk, dt

    dk = qk % domain()
    dt = theta0 % domain()

    call report(.not. dk % same_as(theta) .and. .not. dt % same_as(k), &
         & "observed and trainable are distinguished by domain alone", &
         & nfail)

  end subroutine check_domains_distinguish

end program learning_level_5
