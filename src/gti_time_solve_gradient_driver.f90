!=====================================================================!
! GTI TIME SOLVE-GRADIENT DRIVER
!
! The first end-to-end GTI composition: from a residual form, a
! functional form, a time graph, functional terms, a design, and
! a design direction, one call produces the functional value and
! the total design-gradient action -
!
!      forward solve
!        -> functional value + vertex seeds + explicit F_xi[eta]
!        -> reverse residual contribution
!        -> value and  dF/dxi[eta] = explicit + residual.
!
! This driver COMPOSES the existing seats and does nothing else:
! it never evaluates a form, never builds a Jacobian, never
! touches a relation except through the forward and reverse
! drivers, and never sees the local calculus at all - its imports
! are the three time-level drivers, the graph, and the
! containers, and nothing deeper is nameable from here. Every
! law of the pieces is enforced by the pieces; this seat adds
! only the two scalar checks its own addition needs.
!
! A failed forward march is an answer, not an error: the result
! reports it, the seeding and reverse passes never run, the seed
! array is never allocated, and the graph keeps everything the
! converged steps wrote.
!
! Seeds live in the result: the caller receives the final
! propagated vertex_seed array - the covectors dF/dq_v at every
! time vertex - beside the scalar value and the gradient action.
!
! The driver carries nothing: no form, no graph, no seeds, no
! solver state, no design, no map.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_time_solve_gradient_drivers

  use iso_fortran_env      , only : dp => REAL64
  use gti_value_buffers    , only : gti_value_buffer
  use gti_design_bundles   , only : gti_design_bundle
  use gti_form_interface   , only : gti_differentiable_form
  use gti_time_graphs      , only : gti_time_graph
  use gti_time_forward_drivers, only : gti_time_forward_driver, &
       & gti_time_forward_options, gti_time_forward_result
  use gti_time_reverse_drivers, only : gti_time_reverse_driver, &
       & gti_time_reverse_options, gti_time_reverse_result
  use gti_time_functional_seed_drivers, only : gti_time_functional, &
       & gti_time_functional_seed_result, gti_time_functional_seed_driver

  implicit none

  private
  public :: gti_time_solve_gradient_options
  public :: gti_time_solve_gradient_result
  public :: gti_time_solve_gradient_driver

  !===================================================================!
  ! The knobs of one solve-gradient pass: exactly the forward and
  ! reverse knobs, and nothing of its own.
  !===================================================================!

  type :: gti_time_solve_gradient_options

     type(gti_time_forward_options) :: forward
     type(gti_time_reverse_options) :: reverse

   contains

     procedure :: validate => options_validate

  end type gti_time_solve_gradient_options

  !===================================================================!
  ! What one pass reports: completion, the forward account, the
  ! functional value, the two scalar actions and their sum, the
  ! final propagated seeds, and every sub-result in full.
  !===================================================================!

  type :: gti_time_solve_gradient_result

     logical :: completed = .false.
     logical :: forward_converged = .false.

     real(dp) :: value = 0.0_dp

     type(gti_value_buffer) :: explicit_design_action
     type(gti_value_buffer) :: residual_design_action
     type(gti_value_buffer) :: total_design_gradient_action

     type(gti_value_buffer), allocatable :: vertex_seed(:)

     type(gti_time_forward_result)         :: forward
     type(gti_time_functional_seed_result) :: functional_seed
     type(gti_time_reverse_result)         :: reverse

  end type gti_time_solve_gradient_result

  !===================================================================!
  ! The stateless composition verb. The types keep their public
  ! singular names; Fortran denies a type its host module's name,
  ! so the module speaks in the plural.
  !===================================================================!

  type :: gti_time_solve_gradient_driver

   contains

     procedure :: solve

  end type gti_time_solve_gradient_driver

contains

  !===================================================================!
  ! The pass's knobs are lawful when the forward and reverse
  ! knobs are.
  !===================================================================!

  pure subroutine options_validate(this)

    class(gti_time_solve_gradient_options), intent(in) :: this

    call this % forward % validate()
    call this % reverse % validate()

  end subroutine options_validate

  !===================================================================!
  ! The composition: march, seed, reverse, add. Each stage is the
  ! existing driver, called whole; a failed march ends the pass
  ! with its account on record and nothing downstream run.
  !===================================================================!

  subroutine solve(this, residual_form, functional_form, graph, functional, &
       & design, design_direction, options, result)

    class(gti_time_solve_gradient_driver), intent(in)    :: this
    class(gti_differentiable_form)       , intent(in)    :: residual_form
    class(gti_differentiable_form)       , intent(in)    :: functional_form
    type(gti_time_graph)                 , intent(inout) :: graph
    type(gti_time_functional)            , intent(in)    :: functional
    type(gti_design_bundle)              , intent(in)    :: design
    type(gti_value_buffer)               , intent(in)    :: design_direction
    type(gti_time_solve_gradient_options), intent(in)    :: options
    type(gti_time_solve_gradient_result) , intent(inout) :: result

    type(gti_time_forward_driver)         :: forward
    type(gti_time_functional_seed_driver) :: seeder
    type(gti_time_reverse_driver)         :: reverse

    real(dp) :: explicit_scalar, residual_scalar

    call options % validate()

    !----------------------------------------------------------------!
    ! The default failed/incomplete state: flags down, value zero,
    ! buffers and sub-results fresh, no seed array.
    !----------------------------------------------------------------!

    result % completed         = .false.
    result % forward_converged = .false.
    result % value             = 0.0_dp

    call result % explicit_design_action % clear()
    call result % residual_design_action % clear()
    call result % total_design_gradient_action % clear()

    if (allocated(result % vertex_seed)) deallocate(result % vertex_seed)

    result % functional_seed = gti_time_functional_seed_result()
    result % reverse         = gti_time_reverse_result()

    !----------------------------------------------------------------!
    ! March. A failed march is an answer: report it and stop -
    ! no seeding, no reverse, no seed array.
    !----------------------------------------------------------------!

    call forward % solve_all(residual_form, graph, design, &
         & options % forward, result % forward)

    result % forward_converged = result % forward % converged

    if (.not. result % forward_converged) then
       return
    end if

    !----------------------------------------------------------------!
    ! Seed: the functional terms become the value, the vertex
    ! seeds, and the explicit design action.
    !----------------------------------------------------------------!

    allocate(result % vertex_seed(graph % num_vertices()))

    call seeder % seed_all(functional_form, graph, functional, design, &
         & design_direction, result % vertex_seed, result % functional_seed)

    result % value = result % functional_seed % value

    !----------------------------------------------------------------!
    ! Reverse: the residual contribution, over the same seeds -
    ! which arrive back fully propagated in the result.
    !----------------------------------------------------------------!

    call reverse % reverse_all(residual_form, graph, design, &
         & design_direction, result % vertex_seed, options % reverse, &
         & result % reverse)

    !----------------------------------------------------------------!
    ! Add: total = explicit + residual, both held to scalarness.
    !----------------------------------------------------------------!

    result % explicit_design_action = result % functional_seed % explicit_design_action
    result % residual_design_action = result % reverse % design_gradient_action

    call scalar_from_buffer(result % explicit_design_action, &
         & 'gti_time_solve_gradient_driver: explicit design action is scalar', &
         & explicit_scalar)
    call scalar_from_buffer(result % residual_design_action, &
         & 'gti_time_solve_gradient_driver: residual design action is scalar', &
         & residual_scalar)

    call add_scalar_buffers(result % explicit_design_action, &
         & result % residual_design_action, &
         & result % total_design_gradient_action)

    result % completed = .true.

  end subroutine solve

  !===================================================================!
  ! One scalar out of a buffer, or the given refusal.
  !===================================================================!

  subroutine scalar_from_buffer(buffer, law, value)

    type(gti_value_buffer), intent(in)  :: buffer
    character(len=*)      , intent(in)  :: law
    real(dp)              , intent(out) :: value

    real(dp), allocatable :: values(:)

    call buffer % get_real(values)

    if (buffer % nentries /= 1 .or. buffer % ncomp /= 1 .or. &
         & size(values) /= 1) then
       error stop law
    end if

    value = values(1)

  end subroutine scalar_from_buffer

  !===================================================================!
  ! One scalar into a buffer.
  !===================================================================!

  subroutine set_scalar(buffer, value)

    type(gti_value_buffer), intent(inout) :: buffer
    real(dp)              , intent(in)    :: value

    call buffer % set_real([value])

  end subroutine set_scalar

  !===================================================================!
  ! The one addition this seat owns: total = a + b, both held to
  ! scalarness by their own laws.
  !===================================================================!

  subroutine add_scalar_buffers(a, b, output)

    type(gti_value_buffer), intent(in)    :: a, b
    type(gti_value_buffer), intent(inout) :: output

    real(dp) :: a_value, b_value

    call scalar_from_buffer(a, &
         & 'gti_time_solve_gradient_driver: explicit design action is scalar', &
         & a_value)
    call scalar_from_buffer(b, &
         & 'gti_time_solve_gradient_driver: residual design action is scalar', &
         & b_value)

    call set_scalar(output, a_value + b_value)

  end subroutine add_scalar_buffers

end module gti_time_solve_gradient_drivers
