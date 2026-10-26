!=====================================================================!
! Module that states diffusion: coefficients in, operator out.
!
! This is level 3 of the stratification, where the physics vocabulary
! is defined and ends. A conduction law specifies the material's
! conductivity, and a set of robin conditions specifies the boundary.
! This module translates both into the assembly's neutral vocabulary,
! scales per face and values per boundary face, and delegates. It
! defines no mathematics, no loops over neighbourhoods and no fitting:
! physics terms in, one compiled operator out,
!
!      scales ·············· keff * area, the conductivity's value
!                            through every face
!      boundary relation ·· two numbers per tagged face, the
!                            eliminated face value as an affine
!                            function of its cell: a fixed value, a
!                            fixed gradient, or any combination
!      the form ············ the caller's choice of form; polynomials
!                            unless specified otherwise
!
!   conduction ------> scales -------------+
!                                          +--> fitted_balance_stencil --> op
!   robin conditions -> boundary relation -+
!                       and known fluxes
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_diffusion

  use util_precision           , only : dp
  use view_directed            , only : directed_graph
  use graph_fractal            , only : graph
  use map_set_store            , only : set_store
  use field_forms              , only : form
  use view_mesh                , only : mesh
  use operation_stencil        , only : stencil
  use operation_fitted_balance , only : fitted_balance_stencil
  use operation_conduction     , only : conduction
  use operation_robin_condition, only : robin_condition
  use field_forms              , only : polynomial_form

  implicit none

  private

  public :: diffusion_stencil

contains

  !===================================================================!
  ! Returns the diffusion operator of a mesh: the conduction law
  ! through every face, the robin conditions on their tagged faces, and
  ! the form of the fit, the caller's or the polynomials.
  !
  ! At a boundary face that fixes a gradient (a = 0) the flux itself is
  ! fixed, c / b, and enters the balance directly; there is nothing at
  ! such a face to fit:
  !
  !       cell                 boundary face
  !   +-----------+
  !   |           |---->  a = 0:  the flux c / b, known
  !   |   phi_p   |
  !   |           |---o   a /= 0: the face value, an affine function
  !   +-----------+               of phi_p, fitted
  !===================================================================!

  function diffusion_stencil(m, law, conditions, shape, rings) result(op)

    type(mesh)           , intent(in)           :: m
    type(conduction)     , intent(in)           :: law
    type(robin_condition), intent(in)           :: conditions(:)
    class(form)          , intent(in), optional :: shape
    integer              , intent(in), optional :: rings

    type(stencil) :: op

    class(form), allocatable :: approximation_form

    ! Each condition's faces are declared, read and discarded inside
    ! this call, so their interpretation is local too. A new identity
    ! per call lets one map store every condition's faces without two
    ! of them describing one set.
    type(graph)     :: members
    type(set_store) :: sets

    real(dp), allocatable :: scales(:)
    real(dp), allocatable :: vb(:), wb(:), values(:), weights(:), flux(:)
    logical , allocatable :: known(:)
    integer               :: k, f, e, ne
    character(len=250)    :: message

    ne = m % num_edges()

    ! The material, through every face, the headless ones included:
    ! the fit reads the scale on a boundary face
    call law % edge_coefficients(m, .false., scales)

    ! The boundary, each condition on its own tagged faces. Both
    ! numbers of the boundary relation are passed: a boundary that
    ! fixes a value and a boundary that fixes a gradient are not the
    ! same boundary, and one number cannot distinguish them.
    allocate(vb(ne), wb(ne), flux(ne), known(ne))
    vb    = 0.0_dp
    wb    = 1.0_dp
    flux  = 0.0_dp
    known = .false.
    do k = 1, size(conditions)
       call conditions(k) % faces(m, sets, members)
       call conditions(k) % boundary_relation(m, sets, members, weights, values)
       do f = 1, sets % num_members_of(members)
          e = sets % member_of(members, f)
          wb(e) = weights(f)
          vb(e) = values(f)
       end do

       ! A boundary that fixes a gradient, a = 0, fixes the flux itself,
       ! c / b, and there is nothing at such a face to fit: the flux
       ! enters the balance directly. A boundary relation fitted there
       ! is correct for the value and incorrect for the slope, and a
       ! cell on that boundary then has an error that no refinement
       ! removes.
       if (conditions(k) % a == 0.0_dp) then
          if (conditions(k) % b == 0.0_dp) then
             write(message,'(a,i0,a)') 'operation_diffusion: conditions(', k, &
                  & ') fixes neither a value nor a gradient (a = 0 and b = 0); &
                  &a boundary condition must fix a value, a gradient, or both'
             error stop trim(message)
          end if
          do f = 1, sets % num_members_of(members)
             e        = sets % member_of(members, f)
             known(e) = .true.
             flux(e)  = conditions(k) % c / conditions(k) % b
          end do
       end if
    end do

    ! The form, chosen or defaulted; then the assembly
    if (present(shape)) then
       allocate(approximation_form, source=shape)
    else
       allocate(approximation_form, source=polynomial_form(dimension=m % dimension))
    end if

    op = fitted_balance_stencil(m, approximation_form, scales, &
         & boundary_values=vb, boundary_weights=wb, rings=rings, &
         & flux_known=known, boundary_flux=flux)

  end function diffusion_stencil

end module operation_diffusion
