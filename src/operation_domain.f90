!=====================================================================!
! Module that states the continuous and the discrete domains of an
! operation law.
!
!   continuous_domain               discrete_domain
!   +------------------+            +------------------------------+
!   | law              |  discrete  | the same law                 |
!   |   its degree,    | ---------> | point_set: the vertex set of |
!   |   coordinates,   |  (a graph) |            one directed graph|
!   |   components,    |            | points:    their number      |
!   |   designs        |            +------------------------------+
!   +------------------+                  |
!                                         v
!                    typed_field_domain(point_set, points, width)
!                    for the state, design, direction, costate,
!                    tangent, forcing, residual and solution fields
!
! A continuous_domain is a law before any graph is chosen: it records
! the law and its coordinate degrees. A discrete_domain is that law
! placed on the vertex set of one directed graph. Every typed field
! built for the law reads the same graph identity and point count from
! it, so that no two of its fields can silently disagree about their
! domain.
!
! The module owns no residual assembly and no time march. It only
! states the boundary between a law and the finite set on which a
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

  !-------------------------------------------------------------------!
  ! Type that models a stated law, independent of any graph of points
  !-------------------------------------------------------------------!

  type :: continuous_domain

     type(expression), private :: law

   contains

     ! The law's shape
     procedure :: equation_degree      => continuous_domain_equation_degree
     procedure :: num_coordinates      => continuous_domain_num_coordinates
     procedure :: num_components       => continuous_domain_num_components
     procedure :: num_designs          => continuous_domain_num_designs
     procedure :: component_at         => continuous_domain_component_at
     procedure :: highest_degree_along => continuous_domain_highest_degree_along

     ! The law placed on a graph
     procedure :: discrete             => continuous_domain_discrete

  end type continuous_domain

  !-------------------------------------------------------------------!
  ! Type that models a stated law placed on the vertex set of one
  ! directed graph
  !-------------------------------------------------------------------!

  type :: discrete_domain

     type(continuous_domain), private :: continuous
     type(graph)            , private :: point_set
     integer                , private :: points = 0

   contains

     ! The law and its shape
     procedure :: law             => discrete_domain_law
     procedure :: num_points      => discrete_domain_num_points
     procedure :: equation_degree => discrete_domain_equation_degree
     procedure :: num_coordinates => discrete_domain_num_coordinates
     procedure :: num_components  => discrete_domain_num_components

     ! The domains of the typed fields
     procedure :: state_fields            => discrete_domain_state_fields
     procedure :: design_fields           => discrete_domain_design_fields
     procedure :: functional_state_fields => discrete_domain_functional_state_fields

  end type discrete_domain

  interface continuous_domain
     module procedure create_continuous_domain
  end interface continuous_domain

contains

  !===================================================================!
  ! Creates the continuous domain of a law
  !===================================================================!

  type(continuous_domain) function create_continuous_domain(law) result(this)

    type(expression), intent(in) :: law

    this % law = law

  end function create_continuous_domain

  !===================================================================!
  ! Returns the law's degree along its first coordinate
  !===================================================================!

  pure integer function continuous_domain_equation_degree(this) result(degree)

    class(continuous_domain), intent(in) :: this

    degree = this % law % equation_degree()

  end function continuous_domain_equation_degree

  !===================================================================!
  ! Returns the number of coordinates the law has a degree along
  !===================================================================!

  pure integer function continuous_domain_num_coordinates(this) result(n)

    class(continuous_domain), intent(in) :: this

    n = this % law % num_coordinates()

  end function continuous_domain_num_coordinates

  !===================================================================!
  ! Returns how many components one point of the state stores, the
  ! multipliers excluded
  !===================================================================!

  pure integer function continuous_domain_num_components(this) result(n)

    class(continuous_domain), intent(in) :: this

    n = this % law % num_components()

  end function continuous_domain_num_components

  !===================================================================!
  ! Returns the number of design coordinates the law reads
  !===================================================================!

  pure integer function continuous_domain_num_designs(this) result(n)

    class(continuous_domain), intent(in) :: this

    n = this % law % num_designs()

  end function continuous_domain_num_designs

  !===================================================================!
  ! Returns the position of the component that is the derivative of
  ! the given order along the given coordinate
  !===================================================================!

  pure integer function continuous_domain_component_at(this, coordinate, order) result(at)

    class(continuous_domain), intent(in) :: this
    integer                 , intent(in) :: coordinate, order

    at = this % law % component_at(coordinate, order)

  end function continuous_domain_component_at

  !===================================================================!
  ! Returns the highest order the law names along a coordinate, or -1
  ! when it names none
  !===================================================================!

  pure integer function continuous_domain_highest_degree_along(this, coordinate) result(degree)

    class(continuous_domain), intent(in) :: this
    integer                 , intent(in) :: coordinate

    degree = this % law % highest_degree_along(coordinate)

  end function continuous_domain_highest_degree_along

  !===================================================================!
  ! Places the law on the vertex set of a directed graph
  !===================================================================!

  type(discrete_domain) function continuous_domain_discrete(this, points) result(domain)

    class(continuous_domain), intent(in) :: this
    class(directed_graph)   , intent(in) :: points

    domain % continuous = this
    domain % point_set  = points % vertex_set()
    domain % points     = points % num_vertices()

  end function continuous_domain_discrete

  !===================================================================!
  ! Returns the law placed on the graph
  !===================================================================!

  type(expression) function discrete_domain_law(this) result(law)

    class(discrete_domain), intent(in) :: this

    law = this % continuous % law

  end function discrete_domain_law

  !===================================================================!
  ! Returns the number of points, the vertices of the graph
  !===================================================================!

  pure integer function discrete_domain_num_points(this) result(n)

    class(discrete_domain), intent(in) :: this

    n = this % points

  end function discrete_domain_num_points

  !===================================================================!
  ! Returns the law's degree along its first coordinate
  !===================================================================!

  pure integer function discrete_domain_equation_degree(this) result(degree)

    class(discrete_domain), intent(in) :: this

    degree = this % continuous % equation_degree()

  end function discrete_domain_equation_degree

  !===================================================================!
  ! Returns the number of coordinates the law has a degree along
  !===================================================================!

  pure integer function discrete_domain_num_coordinates(this) result(n)

    class(discrete_domain), intent(in) :: this

    n = this % continuous % num_coordinates()

  end function discrete_domain_num_coordinates

  !===================================================================!
  ! Returns how many components one point of the state stores, the
  ! multipliers excluded
  !===================================================================!

  pure integer function discrete_domain_num_components(this) result(n)

    class(discrete_domain), intent(in) :: this

    n = this % continuous % num_components()

  end function discrete_domain_num_components

  !===================================================================!
  ! Returns the domain of a state field: a value of every component
  ! at every point of the graph
  !
  !              component 1   ...   component c
  !   point 1        .                    .
  !    ...           .                    .
  !   point n        .                    .
  !===================================================================!

  type(typed_field_domain) function discrete_domain_state_fields(this) result(fields)

    class(discrete_domain), intent(in) :: this

    fields = typed_field_domain(this % point_set, this % points, this % num_components())

  end function discrete_domain_state_fields

  !===================================================================!
  ! Returns the domain of a design field: a value of every design
  ! coordinate at every point of the graph
  !===================================================================!

  type(typed_field_domain) function discrete_domain_design_fields(this) result(fields)

    class(discrete_domain), intent(in) :: this

    fields = typed_field_domain(this % point_set, this % points, this % continuous % num_designs())

  end function discrete_domain_design_fields

  !===================================================================!
  ! Returns the domain of a functional's state field, which is that of
  ! a state field
  !===================================================================!

  type(typed_field_domain) function discrete_domain_functional_state_fields(this) result(fields)

    class(discrete_domain), intent(in) :: this

    fields = this % state_fields()

  end function discrete_domain_functional_state_fields

end module operation_domain
