!=====================================================================!
! The refusals that must die at the time graph seat, one per
! invocation:
!
!      novertices    a graph with no vertices at all
!      buildidx0     materializing relation 0
!      buildidxhigh  materializing past the relations
!      novertextuple a relation touching no vertices
!      vertexidx0    a relation naming vertex 0
!      vertexidxhigh a relation naming a vertex past the graph
!      unknown0      a relation whose unknown position is 0
!      unknownhigh   a relation whose unknown position exceeds its
!                    arity
!      norules       a relation whose motif has no rows
!      noweights     a motif row with no weights
!      weightcount   a motif row weighted thrice over two vertices
!      missingq      a referenced vertex whose sample holds no q
!
! Every case must error stop; a case that returns is a failure of
! the suite.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program refusal

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use class_graph_field    , only : field
  use gti_time_local_schemes, only : gti_time_sample
  use gti_time_motif_builders, only : gti_time_motif_builder
  use gti_time_graphs      , only : gti_time_graph

  implicit none

  type(graph) :: states
  type(field) :: q1_field, z2_field

  type(gti_time_motif_builder)       :: builder
  type(gti_time_graph)               :: time_graph
  type(gti_time_sample), allocatable :: samples(:)
  character(len=32) :: which

  call get_command_argument(1, which)

  call states % declare()

  q1_field = field('q1', states, 3)
  call q1_field % set_real_vector([1.0_dp, 2.0_dp, 4.0_dp])
  z2_field = field('placeholder', states, 3)
  call z2_field % set_real_vector([0.0_dp, 0.0_dp, 0.0_dp])

  !-------------------------------------------------------------------!
  ! The lawful two-vertex, one-relation graph every case corrupts.
  !-------------------------------------------------------------------!

  allocate(time_graph % vertex(2))
  allocate(time_graph % vertex(1) % sample % state % component(1))
  time_graph % vertex(1) % sample % state % component(1) % value = q1_field
  allocate(time_graph % vertex(2) % sample % state % component(1))
  time_graph % vertex(2) % sample % state % component(1) % value = z2_field

  allocate(time_graph % relation(1))
  call builder % bdf_uniform(1, 0.5_dp, time_graph % relation(1) % motif)
  time_graph % relation(1) % sample_vertex  = [1, 2]
  time_graph % relation(1) % unknown_sample = 2

  select case (trim(which))

  case ('novertices')

     deallocate(time_graph % vertex)
     call time_graph % validate()

  case ('buildidx0')

     call time_graph % build_samples(0, samples)

  case ('buildidxhigh')

     call time_graph % build_samples(2, samples)

  case ('novertextuple')

     deallocate(time_graph % relation(1) % sample_vertex)
     call time_graph % validate()

  case ('vertexidx0')

     time_graph % relation(1) % sample_vertex = [0, 2]
     call time_graph % validate()

  case ('vertexidxhigh')

     time_graph % relation(1) % sample_vertex = [1, 9]
     call time_graph % validate()

  case ('unknown0')

     time_graph % relation(1) % unknown_sample = 0
     call time_graph % validate()

  case ('unknownhigh')

     time_graph % relation(1) % unknown_sample = 3
     call time_graph % validate()

  case ('norules')

     deallocate(time_graph % relation(1) % motif % rule)
     call time_graph % validate()

  case ('noweights')

     deallocate(time_graph % relation(1) % motif % rule(1) % weight)
     call time_graph % validate()

  case ('weightcount')

     time_graph % relation(1) % motif % rule(1) % weight = [1.0_dp, 0.0_dp, 0.0_dp]
     call time_graph % validate()

  case ('missingq')

     deallocate(time_graph % vertex(2) % sample % state % component)
     call time_graph % validate()

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(which)

end program refusal
