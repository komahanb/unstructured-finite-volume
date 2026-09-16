!=====================================================================!
! Continuous and discrete domains for an operation law.
!
! A continuous_domain is an expression before any graph is chosen. It
! records the mathematical law and its coordinate degrees.
!
! A discrete_domain is that continuous law placed on the vertex set of
! one directed graph. It records the same graph identity and point
! count to every typed field constructed for that law, so a state,
! design, direction, costate, tangent, forcing, residual, and solution
! field cannot silently disagree about their domain.
!
! This module owns no residual assembly and no time-march logic. It
! only states the boundary between a law and the finite set on which a
! caller evaluates it.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_domain

  use graph_fractal       , only : graph
  use view_directed       , only : directed_graph
  use operation_expression, only : expression
  use field_stored        , only : typed_field_domain

  implicit none

  private
  public :: continuous_domain, discrete_domain

  !===================================================================!
  ! A stated law, independent of any point graph.
  !===================================================================!

  type :: continuous_domain

     type(expression), private :: law

   contains

     procedure :: equation_degree      => continuous_domain_equation_degree
     procedure :: num_coordinates      => continuous_domain_num_coordinates
     procedure :: num_components       => continuous_domain_num_components
     procedure :: component_at         => continuous_domain_component_at
     procedure :: highest_degree_along => continuous_domain_highest_degree_along
     procedure :: discrete             => continuous_domain_discrete

  end type continuous_domain

  !===================================================================!
  ! A stated law placed on the vertex set of one directed graph.
  !===================================================================!

  type :: discrete_domain

     type(continuous_domain), private :: continuous
     type(graph)            , private :: point_set
     integer                , private :: points = 0

   contains

     procedure :: law            => discrete_domain_law
     procedure :: num_points     => discrete_domain_num_points
     procedure :: equation_degree => discrete_domain_equation_degree
     procedure :: num_coordinates => discrete_domain_num_coordinates
     procedure :: num_components => discrete_domain_num_components

     procedure :: state_fields            => discrete_domain_state_fields
     procedure :: design_fields           => discrete_domain_design_fields
     procedure :: functional_state_fields => discrete_domain_functional_state_fields

  end type discrete_domain

  interface continuous_domain
     module procedure create_continuous_domain
  end interface continuous_domain

contains

  type(continuous_domain) function create_continuous_domain(law) result(this)

    type(expression), intent(in) :: law

    this % law = law

  end function create_continuous_domain

  pure integer function continuous_domain_equation_degree(this) result(degree)

    class(continuous_domain), intent(in) :: this

    degree = this % law % equation_degree()

  end function continuous_domain_equation_degree

  pure integer function continuous_domain_num_coordinates(this) result(n)

    class(continuous_domain), intent(in) :: this

    n = this % law % num_coordinates()

  end function continuous_domain_num_coordinates

  pure integer function continuous_domain_num_components(this) result(n)

    class(continuous_domain), intent(in) :: this

    n = this % law % num_components()

  end function continuous_domain_num_components

  pure integer function continuous_domain_component_at(this, coordinate, order) result(at)

    class(continuous_domain), intent(in) :: this
    integer                 , intent(in) :: coordinate, order

    at = this % law % component_at(coordinate, order)

  end function continuous_domain_component_at

  pure integer function continuous_domain_highest_degree_along(this, coordinate) result(degree)

    class(continuous_domain), intent(in) :: this
    integer                 , intent(in) :: coordinate

    degree = this % law % highest_degree_along(coordinate)

  end function continuous_domain_highest_degree_along

  type(discrete_domain) function continuous_domain_discrete(this, points) result(domain)

    class(continuous_domain), intent(in) :: this
    class(directed_graph)   , intent(in) :: points

    domain % continuous = this
    domain % point_set  = points % vertex_set()
    domain % points     = points % num_vertices()

  end function continuous_domain_discrete

  type(expression) function discrete_domain_law(this) result(law)

    class(discrete_domain), intent(in) :: this

    law = this % continuous % law

  end function discrete_domain_law

  pure integer function discrete_domain_num_points(this) result(n)

    class(discrete_domain), intent(in) :: this

    n = this % points

  end function discrete_domain_num_points

  pure integer function discrete_domain_equation_degree(this) result(degree)

    class(discrete_domain), intent(in) :: this

    degree = this % continuous % equation_degree()

  end function discrete_domain_equation_degree

  pure integer function discrete_domain_num_coordinates(this) result(n)

    class(discrete_domain), intent(in) :: this

    n = this % continuous % num_coordinates()

  end function discrete_domain_num_coordinates

  pure integer function discrete_domain_num_components(this) result(n)

    class(discrete_domain), intent(in) :: this

    n = this % continuous % num_components()

  end function discrete_domain_num_components

  type(typed_field_domain) function discrete_domain_state_fields(this) result(fields)

    class(discrete_domain), intent(in) :: this

    fields = typed_field_domain(this % point_set, this % points, this % num_components())

  end function discrete_domain_state_fields

  type(typed_field_domain) function discrete_domain_design_fields(this) result(fields)

    class(discrete_domain), intent(in) :: this

    fields = typed_field_domain(this % point_set, this % points)

  end function discrete_domain_design_fields

  type(typed_field_domain) function discrete_domain_functional_state_fields(this) result(fields)

    class(discrete_domain), intent(in) :: this

    fields = this % state_fields()

  end function discrete_domain_functional_state_fields

end module operation_domain
