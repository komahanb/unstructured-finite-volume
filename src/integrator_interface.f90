#include "scalar.fpp"

!=====================================================================!
! Parent class for solving n-th order differential equations. Specific
! integrators extend this class.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module integrator_interface

  use dynamic_physics_interface, only : dynamics
  !use dynamic_analysis_interface, only : dynamic_analysis
  
  implicit none

  private
  public ::  integrator

  !-------------------------------------------------------------------!
  ! Define the type
  !-------------------------------------------------------------------!
  
  type, abstract :: integrator !, extends(dynamic_analysis)

     class(dynamics), allocatable :: system
     type(scalar)   , allocatable :: time(:)  ! time values (steps)
     type(scalar)   , allocatable :: U(:,:,:) ! state varibles (steps, deriv_ord, nvars)
     type(scalar)   , allocatable :: X(:,:)   ! xpoints
     type(scalar)                 :: tinit
     type(scalar)                 :: tfinal
     type(scalar)                 :: h
     type(logical)                :: implicit
     type(integer)                :: num_stages
     type(integer)                :: num_time_steps
     type(integer)                :: total_num_steps

   contains

     procedure                                    :: construct, destruct
     procedure                                    :: solve
     procedure(step_interface), deferred          :: step
     procedure(get_bandwidth_interface), deferred :: get_bandwidth
     procedure                                    :: write_solution
     procedure                                    :: to_string 
     
  end type integrator

  !-------------------------------------------------------------------!
  ! Define interfaces to deferred procedures
  !-------------------------------------------------------------------!
  
  interface

     impure subroutine step_interface(this, t, u, h, p, ierr)

       import integrator

       class(integrator) , intent(inout) :: this
       type(scalar)      , intent(inout) :: t(:)
       type(scalar)      , intent(inout) :: u(:,:,:)
       type(integer)     , intent(in)    :: p
       type(scalar)      , intent(in)    :: h
       type(integer)     , intent(out)   :: ierr

     end subroutine step_interface

     pure type(integer) function get_bandwidth_interface(this, time_index) result(width)

       import integrator

       class(integrator), intent(in) :: this
       type(integer)    , intent(in) :: time_index

     end function get_bandwidth_interface
     
  end interface

contains
   
  !===================================================================!
  ! Base class constructor logic
  !===================================================================!

  subroutine construct(this, system, tinit, tfinal, h, implicit, num_stages)

    class(integrator) , intent(inout) :: this
    class(dynamics)   , intent(in)    :: system
    type(scalar)      , intent(in)    :: tinit, tfinal
    type(scalar)      , intent(in)    :: h
    type(logical)     , intent(in)    :: implicit
    type(integer)     , intent(in)    :: num_stages

    ! Set parameters
    allocate(this % system, source = system)
    this % tinit           = tinit
    this % tfinal          = tfinal
    this % h               = h    
    this % num_time_steps  = floor((this % tfinal - this % tinit)/this % h) + 1
    this % num_stages      = num_stages
    this % total_num_steps = this % num_time_steps*(this % num_stages+1) - this % num_stages ! the initial step does not have stages
    this % implicit        = implicit

  end subroutine construct

  !===================================================================!
  ! Base class destructor
  !===================================================================!
  
  pure subroutine destruct(this)

    class(integrator), intent(inout) :: this

    ! Clear global states and time
    if(allocated(this % U)) deallocate(this % U)
    if(allocated(this % time)) deallocate(this % time)
    if(allocated(this % system)) deallocate(this % system)
    if(allocated(this % X)) deallocate(this % x)
    
  end subroutine destruct

  !===================================================================!
  ! Write solution to file
  !===================================================================!

  subroutine write_solution(this, filename)

    class(integrator)             :: this
    character(len=*), intent(in)  :: filename
    character(len=7), parameter   :: directory = "output/"
    character(len=:), allocatable :: path
    character(len=:), allocatable :: new_name
    integer                       :: k, j, i, ierr

    ! Open resource
    path = trim(filename)

    open(unit=90, file=trim(path), iostat= ierr)
    if (ierr .ne. 0) then
       write(*,'("  >> Opening file ", 39A, " failed")') path
       return
    end if
    
!!$    ! Write data
!!$    loop_time: do k = 1, this % total_num_steps
!!$       write(90, *) this % time(k), this % U (k,1,:)
!!$    end do loop_time
    
    ! Write data
    write(90, *) "time ", "x ", "y ", "z ",  "state "
    loop_vars : do j = 1, this % system % get_num_state_vars()
       loop_time: do k = 1, this % num_time_steps !total_num_steps !
          write(90, *) this % time((this % num_stages+1)*k - this % num_stages), &
               & this % system % x(1,j), &
               & this % system % x(2,j), &
               & this % system % x(3,j), &
               & this % U ((this % num_stages+1)*k - this % num_stages, 1, j)
       end do loop_time
    end do loop_vars
    
    ! Close resource
    close(90)
    
  end subroutine write_solution

  !===================================================================!
  ! Time integration logic
  !===================================================================!

  impure subroutine solve( this )
  
    class(integrator), intent(inout) :: this
    
    integer :: k, p
    integer :: ierr
    
    ! State and time history
    if (allocated(this % time)) deallocate(this%time)
    if (allocated(this % U)) deallocate(this%U)
    
    allocate(this % time(this % total_num_steps))
    this % time = 0.0d0
    this % time(1) = this % tinit
    
    allocate( this % U( &
         & this % total_num_steps, &
         & this % system % get_differential_order() + 1, &
         & this % system % get_num_state_vars() &
         & ))
    this % U = 0.0d0   

    ! Get the initial condition
    call this % system % get_initial_condition(this % U(1,:,:), this % X)
    
    ! March in time
    time: do k = 2, this % total_num_steps

       p = this % get_bandwidth(k)

       call this % step(this % time(k-p:k) , &
            & this % U(k-p:k,:,:), &
            & this % h, &
            & p, &
            & ierr)   

       !! if (k.eq.50) stop

    end do time
    
  end subroutine solve
  
  !===================================================================!
  ! Prints important fields of the class
  !===================================================================!
  
  subroutine to_string(this)
    
    class(integrator), intent(in) :: this
    
    print '("  >> Physical System      : " ,A10)' , this % system % get_description()
    print '("  >> Num state variables  : " ,i8)'  , this % system % get_num_state_vars()
    print '("  >> Equation order       : " ,i4)'  , this % system % get_differential_order()
    print '("  >> Start time           : " ,F8.3)', this % tinit
    print '("  >> End time             : " ,F8.3)', this % tfinal
    print '("  >> Step size            : " ,E9.3)', this % h
    print '("  >> Number of time steps : " ,i10)' , this % num_time_steps
    print '("  >> Number of stages     : " ,i10)' , this % num_stages
    print '("  >> Tot. Number of steps : " ,i10)' , this % total_num_steps

  end subroutine to_string
  
end module integrator_interface
