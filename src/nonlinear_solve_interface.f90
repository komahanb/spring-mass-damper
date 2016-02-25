!=====================================================================!
! A module that wraps all the data used in nonlinear solve
!
! Author :  Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module nonlinear_solve_interface

  use precision
  ! use abstract_linear_solve

  implicit none

  private

  public :: nonlinear_solve_bean

  type, abstract :: nonlinear_solve_bean

     !     private

     !----------------------------------------------------------------!
     ! Basic setup variables
     !----------------------------------------------------------------!

     integer :: max_nonlinear_iters = 50
     integer :: nvars = 1
     integer :: iter_num = 1

     !----------------------------------------------------------------!
     ! Stopping criteria
     !----------------------------------------------------------------!

     ! absolute tolerances
     real(dp) :: atol_rnrm = 1.0d-12
     real(dp) :: atol_unrm = 1.0d-12

     ! relative tolerances
     real(dp) :: rtol_rnrm = 1.0d-12
     real(dp) :: rtol_unrm = 1.0d-12

     ! actual norm values
     real(dp), dimension(:), allocatable :: rnrm, unrm

     !----------------------------------------------------------------!
     ! Print and solution write control
     !----------------------------------------------------------------!

     integer :: filenum              = 10
     logical :: write_nonlinear_details = .true.
     logical :: exit_on_failure      = .false. 

     logical :: converged = .false. 

   contains

     procedure :: get_max_nonlinear_iters, set_max_nonlinear_iters
     procedure :: get_num_vars, set_num_vars

     procedure :: get_atol_unrm, set_atol_unrm
     procedure :: get_atol_rnrm, set_atol_rnrm

     procedure :: get_rtol_unrm, set_rtol_unrm
     procedure :: get_rtol_rnrm, set_rtol_rnrm

     procedure :: get_file_number, set_file_number

     procedure :: is_write_nonlinear_details, set_write_nonlinear_details
     procedure :: is_exit_on_failure, set_exit_on_failure
     procedure :: is_nonlinear_solve_converged, set_nonlinear_solve_converged

  end type nonlinear_solve_bean

  type, abstract, extends(nonlinear_solve_bean) :: abstract_nonlinear_solve

   contains

     ! will provide implementation for just the solution
     procedure :: solve

     ! define all the deferred procedures
     procedure(init_nonlinear_solve), deferred :: init
     procedure(finish_nonlinear_solve), deferred :: finish
     procedure(work_nonlinear_solve), deferred :: work
     procedure(linear_nonlinear_solve), deferred :: linear_solve
     procedure(check_stop_nonlinear_solve), deferred :: check_stop
     procedure(write_solution_nonlinear_solve), deferred :: write_solution
     procedure(read_inputfile_nonlinear_solve), deferred :: read_input_file

  end type abstract_nonlinear_solve

  !-------------------------------------------------------------------
  ! define all the interfaces here
  !-------------------------------------------------------------------
  
  interface

     subroutine init_nonlinear_solve(this)
       import abstract_nonlinear_solve
       class(abstract_nonlinear_solve) :: this
     end subroutine init_nonlinear_solve

     subroutine finish_nonlinear_solve(this)
       import abstract_nonlinear_solve
       class(abstract_nonlinear_solve) :: this
     end subroutine finish_nonlinear_solve

     subroutine work_nonlinear_solve(this)
       import abstract_nonlinear_solve
       class(abstract_nonlinear_solve) :: this
     end subroutine work_nonlinear_solve

     subroutine linear_nonlinear_solve(this)
       import abstract_nonlinear_solve
       class(abstract_nonlinear_solve) :: this
     end subroutine linear_nonlinear_solve

     subroutine check_stop_nonlinear_solve(this)
       import abstract_nonlinear_solve
       class(abstract_nonlinear_solve) :: this
     end subroutine check_stop_nonlinear_solve

     subroutine write_solution_nonlinear_solve(this)
       import abstract_nonlinear_solve
       class(abstract_nonlinear_solve) :: this
     end subroutine write_solution_nonlinear_solve
     
     subroutine read_inputfile_nonlinear_solve(this)
       import abstract_nonlinear_solve
       class(abstract_nonlinear_solve) :: this
     end subroutine read_inputfile_nonlinear_solve
     
  end interface

contains

  !-------------------------------------------------------------------!
  ! Get max_nonlinear_iters
  !-------------------------------------------------------------------!

  integer function get_max_nonlinear_iters(this)

    class(nonlinear_solve_bean) :: this

    get_max_nonlinear_iters = this % max_nonlinear_iters

  end function get_max_nonlinear_iters

  !-------------------------------------------------------------------!
  ! Set max_nonlinear_iters
  !-------------------------------------------------------------------!

  subroutine set_max_nonlinear_iters(this, max_nonlinear_iters)

    class(nonlinear_solve_bean) :: this
    integer :: max_nonlinear_iters

    this % max_nonlinear_iters =  max_nonlinear_iters

  end subroutine set_max_nonlinear_iters

  !-------------------------------------------------------------------!
  ! Get number of variables
  !-------------------------------------------------------------------!

  integer function get_num_vars(this)

    class(nonlinear_solve_bean) :: this

    get_num_vars = this % nvars

  end function get_num_vars

  !-------------------------------------------------------------------!
  ! Set number of variables
  !-------------------------------------------------------------------!

  subroutine set_num_vars(this, nvars)

    class(nonlinear_solve_bean) :: this
    integer :: nvars

    this % nvars = nvars

  end subroutine set_num_vars

  !-------------------------------------------------------------------!
  ! Get absolute tolerance of update norm
  !-------------------------------------------------------------------!

  real(dp) function get_atol_unrm(this)

    class(nonlinear_solve_bean) :: this

    get_atol_unrm  = this % atol_unrm

  end function get_atol_unrm

  !-------------------------------------------------------------------!
  ! Set absolute tolerance of update norm
  !-------------------------------------------------------------------!

  subroutine set_atol_unrm(this, atol_unrm)

    class(nonlinear_solve_bean) :: this
    real(dp) :: atol_unrm

    this % atol_unrm = atol_unrm

  end subroutine set_atol_unrm

  !-------------------------------------------------------------------!
  ! Get absolute tolerance of residual norm
  !-------------------------------------------------------------------!

  real(dp) function get_atol_rnrm(this)

    class(nonlinear_solve_bean) :: this

    get_atol_rnrm  = this % atol_rnrm

  end function get_atol_rnrm

  !-------------------------------------------------------------------!
  ! Set absolute tolerance of residual norm
  !-------------------------------------------------------------------!

  subroutine set_atol_rnrm(this, atol_rnrm)

    class(nonlinear_solve_bean) :: this
    real(dp) :: atol_rnrm

    this % atol_rnrm = atol_rnrm

  end subroutine set_atol_rnrm

  !-------------------------------------------------------------------!
  ! Get relative tolerance of update norm
  !-------------------------------------------------------------------!

  real(dp) function get_rtol_unrm(this)

    class(nonlinear_solve_bean) :: this

    get_rtol_unrm  = this % rtol_unrm

  end function get_rtol_unrm

  !-------------------------------------------------------------------!
  ! Set relative tolerance of update norm
  !-------------------------------------------------------------------!

  subroutine set_rtol_unrm(this, rtol_unrm)

    class(nonlinear_solve_bean) :: this
    real(dp) :: rtol_unrm

    this % rtol_unrm = rtol_unrm

  end subroutine set_rtol_unrm

  !-------------------------------------------------------------------!
  ! Get relative tolerance of residual norm
  !-------------------------------------------------------------------!

  real(dp) function get_rtol_rnrm(this)

    class(nonlinear_solve_bean) :: this

    get_rtol_rnrm  = this % rtol_rnrm

  end function get_rtol_rnrm

  !-------------------------------------------------------------------!
  ! Get the output stream filenum
  !-------------------------------------------------------------------!
  
  integer function get_file_number(this)

    class(nonlinear_solve_bean) :: this

    get_file_number = this % filenum

  end function get_file_number

  !-------------------------------------------------------------------!
  ! Should write the nonlinear iteration details at each time step?
  !-------------------------------------------------------------------!
  
  logical function is_write_nonlinear_details(this)

    class(nonlinear_solve_bean) :: this

    is_write_nonlinear_details = this % write_nonlinear_details
    
  end function is_write_nonlinear_details

  !-------------------------------------------------------------------!
  ! Should write the nonlinear iteration details at each time step?
  !-------------------------------------------------------------------!
  
  logical function is_exit_on_failure(this)
    
    class(nonlinear_solve_bean) :: this

    is_exit_on_failure = this % exit_on_failure

  end function is_exit_on_failure

  !-------------------------------------------------------------------!
  ! Is the nonlinear solution converged at this time step
  !-------------------------------------------------------------------!
  
  logical function is_nonlinear_solve_converged(this)
    
    class(nonlinear_solve_bean) :: this

    is_nonlinear_solve_converged = this % converged

  end function is_nonlinear_solve_converged

  !-------------------------------------------------------------------!
  ! Set relative tolerance of residual norm
  !-------------------------------------------------------------------!

  subroutine set_rtol_rnrm(this, rtol_rnrm)

    class(nonlinear_solve_bean) :: this
    real(dp) :: rtol_rnrm

    this % rtol_rnrm = rtol_rnrm

  end subroutine set_rtol_rnrm

  !-------------------------------------------------------------------!
  ! Set the output file number
  !-------------------------------------------------------------------!
  
  subroutine set_file_number(this, filenum)

    class(nonlinear_solve_bean) :: this
    integer :: filenum

    this % filenum = filenum

  end subroutine set_file_number

  !-------------------------------------------------------------------!
  ! Set the print control for writing the details of nonlinear solve
  !-------------------------------------------------------------------!

  subroutine set_write_nonlinear_details(this, write_nonlinear_details)

    class(nonlinear_solve_bean) :: this
    logical :: write_nonlinear_details

    this % write_nonlinear_details = write_nonlinear_details

  end subroutine set_write_nonlinear_details

  !-------------------------------------------------------------------!
  ! Set whether or not to exit when not converged
  !-------------------------------------------------------------------!

  subroutine set_exit_on_failure(this, exit_on_failure)

    class(nonlinear_solve_bean) :: this
    logical :: exit_on_failure

    this % exit_on_failure = exit_on_failure

  end subroutine set_exit_on_failure

  !-------------------------------------------------------------------!
  ! Did the nonlinear solve converge
  !-------------------------------------------------------------------!

  subroutine set_nonlinear_solve_converged(this, converged)

    class(nonlinear_solve_bean) :: this
    logical :: converged

    this % converged = converged

  end subroutine set_nonlinear_solve_converged

  !-------------------------------------------------------------------!
  ! Routine that performs the nonlinear solution process
  !-------------------------------------------------------------------!

  subroutine solve(this)

    class(abstract_nonlinear_solve) :: this

    ! perform initialization tasks
    call this % init()

    ! perform newton solve
    call this % work()

    ! perform finalization tasks
    call this % finish()

  end subroutine solve

!!$  ! shouln't be here
!!$  subroutine read_input_file(this)
!!$    class(nonlinear_solve) :: this
!!$  end subroutine read_input_file
  
end module nonlinear_solve_interface
