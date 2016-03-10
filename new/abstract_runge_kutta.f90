!=====================================================================!
! An abstract module for Runge-Kutta integration
!=====================================================================!
! o This module is suitable for the First order ODEs: q'=f(q(t),t)
!
! o User needs to provide the implementation for the function 'f' in
!   the module
! 
! o Yet to implement the multivariate case
!=====================================================================!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module abstract_runge_kutta

  implicit none

  private

  public :: RK
  
  !-------------------------------------------------------------------!
  ! Abstract Runge-Kutta type
  !-------------------------------------------------------------------!

  type, abstract :: RK

     integer :: num_stages = 1  ! default number of stages
     integer :: order           ! order of accuracy

     real(8) :: h = 0.1d0       ! default step size (will reconsider when implementing adaptive step size)
     real(8) :: time            ! scalar to track integration time
     
     ! The Butcher Tableau 
     real(8), dimension(:,:), allocatable :: A ! forms the coeff matrix
     real(8), dimension(:)  , allocatable :: B ! multiplies the state derivatives
     real(8), dimension(:)  , allocatable :: C ! multiplies the time

     ! The stage time and its corresponding derivatives
     real(8), dimension(:)  , allocatable :: T ! the corresponding stage time
     real(8), dimension(:)  , allocatable :: Y ! the corresponding state
     real(8), dimension(:)  , allocatable :: K ! the stage derivatives K = F(T,Y)

     real(8), dimension(:)  , allocatable :: R ! stage residual
     real(8), dimension(:,:), allocatable :: J ! stage jacobian
     
     ! The form of the governing equation
     logical :: descriptor_form = .false.
     
   contains

     ! Implemented common procedures (visible to the user)
     procedure :: initialize, finalize, integrate
     
     ! Implemented procedures (not callable by the user)
     procedure, private :: update_states
     procedure, private :: reset_stage_values
     procedure, private :: check_butcher_tableau

     ! Deferred common procedures
     procedure(compute_stage_values_interface), private, deferred :: compute_stage_values
     procedure(buthcher_interface), private, deferred :: setup_butcher_tableau

  end type RK

  !-------------------------------------------------------------------!
  ! Interfaces for deferred specialized procedures 
  !-------------------------------------------------------------------!

  interface

     !----------------------------------------------------------------!
     ! Interface for finding the stage derivatives at each time step
     !----------------------------------------------------------------!

     subroutine compute_stage_values_interface(this, k, q)
       import RK
       class(RK) :: this
       integer, intent(in) :: k 
       real(8), intent(in), dimension(:) :: q
     end subroutine compute_stage_values_interface

     !----------------------------------------------------------------!
     ! Interface for setting the Butcher tableau for each type of RK
     ! scheme
     !----------------------------------------------------------------!

     subroutine buthcher_interface(this)
       import RK
       class(RK) :: this
     end subroutine buthcher_interface

  end interface

contains

  !-------------------------------------------------------------------!
  ! Initialize the dirk datatype and construct the tableau
  !-------------------------------------------------------------------!

  subroutine initialize(this, tinit, num_stages, h)

    class(RK) :: this
    integer, OPTIONAL, intent(in) :: num_stages
    real(8), OPTIONAL, intent(in) :: tinit
    real(8), OPTIONAL, intent(in) :: h

    !-----------------------------------------------------------------!
    ! set the initial time
    !-----------------------------------------------------------------!
    
    if (present(tinit)) then
       this % time = tinit
    else
       print '("Using default start time : ",F8.3)', this % time
    end if
    
    !-----------------------------------------------------------------!
    ! set the order of integration
    !-----------------------------------------------------------------!
    
    if (present(num_stages)) then
       this % num_stages = num_stages
    else
       print '("Using default number of stages : ",i4)', this % num_stages
    end if

    !-----------------------------------------------------------------!
    ! set the user supplied initial step size
    !-----------------------------------------------------------------!
    
    if (present(h)) then
       this % h = h 
    else
       print '("Using default step size h : ", E9.3)', this % h
    end if
    
    !-----------------------------------------------------------------!
    ! allocate space for the tableau
    !-----------------------------------------------------------------!

    allocate(this % A(this % num_stages, this % num_stages))
    this % A = 0.0d0

    allocate(this % B(this % num_stages))    
    this % B = 0.0d0

    allocate(this % C(this % num_stages))
    this % C = 0.0d0

    !-----------------------------------------------------------------!
    ! allocate space for the stage derivatives
    !-----------------------------------------------------------------!

    allocate(this % K(this % num_stages))
    this % K = 0.0d0

    !-----------------------------------------------------------------!
    ! allocate space for the stage state
    !-----------------------------------------------------------------!

    allocate(this % Y(this % num_stages))
    this % Y = 0.0d0

    !-----------------------------------------------------------------!
    ! allocate space for the stage time
    !-----------------------------------------------------------------!

    allocate(this % T(this % num_stages))
    this % T = 0.0d0

    !-----------------------------------------------------------------!
    ! allocate space for the stage time
    !-----------------------------------------------------------------!

    allocate(this % R(this % num_stages))
    this % R = 0.0d0

    !-----------------------------------------------------------------!
    ! allocate space for the stage time
    !-----------------------------------------------------------------!

    allocate(this % J(this % num_stages, this % num_stages))
    this % J = 0.0d0

    !-----------------------------------------------------------------!
    ! this subroutine puts values into the Butcher tableau
    !-----------------------------------------------------------------!

    call this % setup_butcher_tableau()
    
    !-----------------------------------------------------------------!
    ! sanity check
    !-----------------------------------------------------------------!
    
    call this % check_butcher_tableau()

  end subroutine initialize

  !-------------------------------------------------------------------!
  ! Routine that checks if the Butcher Tableau entries are valid for
  ! the chosen number of stages/order
  !--------------------------------------------------------------------!
  subroutine check_butcher_tableau(this)

    class(RK) :: this
    integer :: i

    do i = 1, this  % num_stages

       if (abs(this % C(i) - sum(this % A(i,:))) .gt. 5.0d-16) then

          print *, "WARNING: sum(A(i,j)) != C(i)", i, this % num_stages

       end if

    end do

    if ((sum(this % B) - 1.0d0) .gt. 5.0d-16) then

       print *, "WARNING: sum(B) != 1", this % num_stages

    end if

  end subroutine check_butcher_tableau

  !-------------------------------------------------------------------!
  ! Deallocate the tableau entries
  !-------------------------------------------------------------------!

  subroutine finalize(this)

    class(RK) :: this

    ! clear butcher's tableau
    if(allocated(this % A)) deallocate(this % A)
    if(allocated(this % B)) deallocate(this % B)
    if(allocated(this % C)) deallocate(this % C)

    ! clear stage values
    if(allocated(this % K)) deallocate(this % K)
    if(allocated(this % T)) deallocate(this % T)
    if(allocated(this % Y)) deallocate(this % Y)

    ! clear the stage residual and jacobian
    if(allocated(this % R)) deallocate(this % R)
    if(allocated(this % J)) deallocate(this % J)

  end subroutine finalize

  !-------------------------------------------------------------------!
  ! Time integration logic
  !-------------------------------------------------------------------!
  ! Input: 
  ! o state arrays q and qdot with initial conditions set at q(1)
  ! o number of steps N
  ! o step size h
  !-------------------------------------------------------------------!
  ! Output:
  ! o q, qdot arrays are modified by the routine
  !-------------------------------------------------------------------!

  subroutine Integrate(this, q, qdot, N)

    class(RK) :: this
    real(8), intent(inout), dimension(:) :: q, qdot
    integer, intent(in) :: N 
    integer :: k

    ! March in time
    march: do k = 2, N + 1

       ! find the stage derivatives at the current step
       call this % compute_stage_values(k, q)

       ! advance the state to the current step
       call this % update_states(k, q, qdot)

       ! set the stage values to zero
       call this % reset_stage_values()

    end do march

  end subroutine Integrate

  !-------------------------------------------------------------------!
  ! Update the states based on RK Formulae
  !-------------------------------------------------------------------!

  subroutine update_states(this, k, q, qdot)

    class(RK) :: this
    integer, intent(in) :: k ! current time step
    real(8), intent(inout), dimension(:) :: q ! actual states
    real(8), intent(inout), dimension(:) :: qdot ! actual state    
    ! real(8), external :: F

    ! increment the time
    this % time = this % time + this % h
    
    ! update q (for first order ODE)
    q(k) = q(k-1) + this % h*sum(this % B * this % K)
    
  end subroutine update_states

  !-------------------------------------------------------------------!
  ! Reset the array to store new stage values at each time step
  !-------------------------------------------------------------------!

  subroutine reset_stage_values(this)

    class(RK) :: this

    ! reset the variables that are computed during each time step    
    this % K = 0.0d0
    this % Y = 0.0d0
    this % T = 0.0d0
    
    this % R = 0.0d0
    this % J = 0.0d0

  end subroutine reset_stage_values

end module abstract_runge_kutta
