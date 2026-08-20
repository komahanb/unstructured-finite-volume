!=====================================================================!
! The diffusion statement: the constitution speaks, an operator
! answers.
!
! LEVEL 3 OF THE STRATIFICATION. This is where the physics words
! live and stop: a conduction law says what the material carries, a
! set of robin conditions says what the walls hold, and this module
! translates both into the assembly's neutral vocabulary - scales
! per face, values per wall - then delegates. It owns no
! mathematics, no loops over neighbourhoods, no fitting: physics
! words in, one compiled operator out,
!
!      scales ·········· keff * area, the conductivity's answer
!                        through every face
!      wall relation ··· two numbers per tagged face, the eliminated
!                        face value as an affine function of its
!                        cell: a held value, a held gradient, or
!                        anything between
!      the form ········ the caller's choice of shape; polynomials
!                        unless said otherwise
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_diffusion

  use iso_fortran_env      , only : dp => REAL64
  use view_directed  , only : directed_graph
  use graph_fractal      , only : graph
  use map_set      , only : set_map
  use map_label    , only : label_map
  use map_inclusion, only : inclusion_map
  use field_forms          , only : form
  use field_stored    , only : stored_field
  use view_mesh     , only : mesh
  use operation_stencil  , only : stencil
  use operation_fitted_balance , only : fitted_balance_stencil
  use operation_conduction     , only : conduction
  use operation_robin_condition, only : robin_condition
  use field_forms       , only : polynomial_form

  implicit none

  private
  public :: diffusion_stencil

contains

  function diffusion_stencil(m, law, conditions, shape) result(op)

    type(mesh)           , intent(in) :: m
    type(conduction)     , intent(in) :: law
    type(robin_condition), intent(in) :: conditions(:)
    class(form), intent(in), optional :: shape

    type(stencil) :: op

    class(form), allocatable :: chosen
    !----------------------------------------------------------------!
    ! Each condition's faces are carved, read and dropped inside this
    ! call, so their interpretation is local too. A fresh identity per
    ! call is what lets one map hold every condition's faces without
    ! two of them claiming to describe one set.
    !----------------------------------------------------------------!

    type(graph)     :: members
    type(set_map)       :: sets
    type(label_map)     :: labels
    type(inclusion_map) :: inclusions

    type(stored_field) :: fa
    real(dp), allocatable :: keff(:), areas(:), scales(:)
    real(dp), allocatable :: vb(:), wb(:), values(:), weights(:)
    integer :: k, f, e, ne

    ne = m % num_edges()

    ! The material, through every face.
    call law % normal_conductivity(m, keff)
    fa = m % face_area()
    call fa % real_vector(areas)
    scales = keff * areas

    ! The walls, each condition on its own tagged faces. Both
    ! numbers of the wall relation travel: a wall that holds a value
    ! and a wall that holds a gradient are not the same wall, and
    ! one number cannot tell them apart.
    allocate(vb(ne), wb(ne))
    vb = 0.0_dp
    wb = 1.0_dp
    do k = 1, size(conditions)
       call conditions(k) % faces(m, sets, labels, inclusions, members)
       call conditions(k) % wall_relation(m, weights, values)
       do f = 1, sets % num_members_of(members)
          e = sets % member_of(members, f)
          wb(e) = weights(f)
          vb(e) = values(f)
       end do
    end do

    ! The shape, chosen or defaulted; then the assembly does the act.
    if (present(shape)) then
       allocate(chosen, source=shape)
    else
       allocate(chosen, source=polynomial_form())
    end if

    op = fitted_balance_stencil(m, chosen, scales, &
         & boundary_values=vb, boundary_weights=wb)

  end function diffusion_stencil

end module operation_diffusion
