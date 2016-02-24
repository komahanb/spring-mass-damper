!=====================================================================!
! Module that extends the abstract system variables and provides
! implementations of some useful initialization routines
!=====================================================================!

module unsteady_problem_class

  use precision
  use interfaces

  implicit none

  private

  public :: system, state
  public :: initialize_simulation, finalize_simulation

  !-------------------------------------------------------------------!
  ! The variables that characterize the unsteady system under study
  ! (any other states and function that are needed other than the
  ! default)
  ! -------------------------------------------------------------------!
  type,  extends(abstract_unsteady_descriptor) :: unsteady_descriptor
   contains
     procedure :: initialize_system_variables  
  end type unsteady_descriptor
  
  !-------------------------------------------------------------------!
  ! The state variables (any other states and function that are needed
  ! other than the default)
  ! -------------------------------------------------------------------!

  type, extends(abstract_state_variables) :: state_variables
   contains
     procedure :: initialize_state_variables
     procedure :: finalize_state_variables
  end type state_variables

  !-------------------------------------------------------------------!
  ! The variables that are shared all over and contains main data
  !-------------------------------------------------------------------!

  type(unsteady_descriptor) :: system
  type(state_variables)     :: state

contains

  !-------------------------------------------------------------------!
  ! Initialization tasks for ALL  variables in the simulation
  !-------------------------------------------------------------------!
  
  subroutine initialize_simulation()
    
    call system % initialize_system_variables()
    call state  % initialize_state_variables()

  end subroutine initialize_simulation

  !-------------------------------------------------------------------!
  ! Finalize all the variables
  !-------------------------------------------------------------------!

  subroutine finalize_simulation()

    call state  % finalize_state_variables()

  end subroutine finalize_simulation

  !-------------------------------------------------------------------!
  ! Initialization tasks for system descriptor variables
  !-------------------------------------------------------------------!
  
  subroutine initialize_system_variables(this)

    class(unsteady_descriptor) :: this
    integer :: num_steps

    ! Use namelist system to read these params

    call system % set_num_vars(1)
    call system % set_start_time(0.0_dp)
    call system % set_end_time(10.0_dp)
    call system % set_step_size(1.0e-3_dp)
    call system % set_unsteady(.true.)

    num_steps = int((system % get_end_time()-system%get_start_time()) &
         &/system % get_step_size()) + 1

    call system % set_num_time_steps(num_steps)

  end subroutine initialize_system_variables

  
  !-------------------------------------------------------------------!
  ! Initialization tasks for state variables in the simulation
  !-------------------------------------------------------------------!

  subroutine initialize_state_variables(this)

    class(state_variables) :: this

    !-----------------------------------------------------------------!
    ! allocate and initialize q
    !-----------------------------------------------------------------!
    if (allocated(this % q)) deallocate(this % q)
    allocate(this % q(system % get_num_time_steps(), system % get_num_vars()))
    this % q = 0.0_dp

    !-----------------------------------------------------------------!
    ! allocate and initialize qdot
    !-----------------------------------------------------------------!
    if (allocated(this % qdot)) deallocate(this % qdot)
    allocate(this % qdot(system % get_num_time_steps(), system % get_num_vars()))
    this % qdot = 0.0_dp

    !-----------------------------------------------------------------!
    ! allocate and initialize qddot
    !-----------------------------------------------------------------!
    if (allocated(this % qddot)) deallocate(this % qddot)
    allocate(this % qddot(system % get_num_time_steps(), system % get_num_vars()))
    this % qddot = 0.0_dp

    !-----------------------------------------------------------------!
    ! allocate and initialize dq (update)
    !-----------------------------------------------------------------!
    if (allocated(this % dq)) deallocate(this % dq)
    allocate(this % dq(system % get_num_vars()))
    this % dq = 0.0_dp

    !-----------------------------------------------------------------!
    ! allocate and initialize R
    !-----------------------------------------------------------------!
    if (allocated(this % R)) deallocate(this % R)
    allocate(this % R(system % get_num_time_steps(), system % get_num_vars()))
    this % R = 0.0_dp

    !-----------------------------------------------------------------!
    ! allocate and initialize dR
    !-----------------------------------------------------------------!
    if (allocated(this % dR)) deallocate(this % dR)
    allocate(this % dR(system % get_num_time_steps(), system % get_num_vars(), system % get_num_vars()))
    this % dR = 0.0_dp

  end subroutine initialize_state_variables

  !-------------------------------------------------------------------!
  ! Finalize the state variables and freeup memory
  !-------------------------------------------------------------------!

  subroutine finalize_state_variables(this)

    class(state_variables) :: this

    write(*, *) "Finalize variables"

    if (allocated(this % q))     deallocate(this % q)
    if (allocated(this % qdot))  deallocate(this % qdot)
    if (allocated(this % qddot)) deallocate(this % qddot)

  end subroutine finalize_state_variables

end module unsteady_problem_class
