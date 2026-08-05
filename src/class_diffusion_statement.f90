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
!      wall values ····· c / a, each condition's eliminated value
!                        on its own tagged faces
!      the form ········ the caller's choice of shape; polynomials
!                        unless said otherwise
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_diffusion_statement

  use iso_fortran_env      , only : dp => REAL64
  use graph_grammar        , only : graph
  use graph_forms          , only : form
  use class_graph_field    , only : field
  use class_graph_mesh     , only : mesh
  use class_graph_stencil  , only : stencil_operator
  use class_fitted_balance , only : fitted_balance_stencil
  use class_conduction     , only : conduction
  use class_robin_condition, only : robin_condition
  use class_polynomial_form, only : polynomial_form

  implicit none

  private
  public :: diffusion_statement

contains

  function diffusion_statement(m, law, conditions, shape) result(op)

    type(mesh)           , intent(in) :: m
    type(conduction)     , intent(in) :: law
    type(robin_condition), intent(in) :: conditions(:)
    class(form), intent(in), optional :: shape

    type(stencil_operator) :: op

    class(form), allocatable :: chosen
    class(graph), allocatable :: members
    type(field) :: fa
    real(dp), allocatable :: keff(:), areas(:), scales(:), vb(:), bw(:)
    integer :: k, f, e, ne

    ne = m % num_edges()

    ! The material, through every face.
    call law % normal_conductivity(m, keff)
    fa = m % face_area()
    call fa % get_real_vector(areas)
    scales = keff * areas

    ! The walls, each condition on its own tagged faces.
    allocate(vb(ne))
    vb = 0.0_dp
    do k = 1, size(conditions)
       call conditions(k) % faces(m, members)
       call conditions(k) % boundary_values(m, bw)
       do f = 1, members % num_vertices()
          e = members % global_vertex_index(f)
          vb(e) = bw(f)
       end do
    end do

    ! The shape, chosen or defaulted; then the assembly does the act.
    if (present(shape)) then
       allocate(chosen, source=shape)
    else
       allocate(chosen, source=polynomial_form())
    end if

    op = fitted_balance_stencil(m, chosen, scales, boundary_values=vb)

  end function diffusion_statement

end module class_diffusion_statement
