!=====================================================================!
! Module that manages the precision of real variables
!=====================================================================!

module precision

  !-------------------------------------------------------------------!
  !  Define constants to manage precision [TUNABLE]
  !-------------------------------------------------------------------!

  integer, parameter :: sp = kind(0.0)    ! single precision
  integer, parameter :: dp = kind(0.0d0)  ! double precision

end module precision

!=====================================================================!
! A module that wraps all the data used in Newton solve
!
! Author :  Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module newton_solve_bean_class

  use precision

  implicit none

  private
  public :: newton_solve_bean

  type newton_solve_bean

     !----------------------------------------------------------------!
     ! Basic setup variables
     !----------------------------------------------------------------!

     integer :: max_newton_iters = 50
     integer :: nvars = 1

     !----------------------------------------------------------------!
     ! Solution variables
     !----------------------------------------------------------------!

     real(dp) :: init_q
     real(dp) :: init_qdot

     real(dp) :: q
     real(dp) :: qdot
     real(dp) :: qddot

     !----------------------------------------------------------------!
     ! Stopping criteria
     !----------------------------------------------------------------!

     ! absolute tolerances
     real(dp) :: atol_rnrm = 1.0d-12
     real(dp) :: atol_unrm = 1.0d-12

     ! relative tolerances
     real(dp) :: rtol_rnrm = 1.0d-6
     real(dp) :: rtol_unrm = 1.0d-6

     ! actual norm values
     real(dp) :: rnrm      = 0.0_dp
     real(dp) :: unrm      = 0.0_dp 

     !----------------------------------------------------------------!
     ! Print and solution write control
     !----------------------------------------------------------------!

     integer :: filenum              = 6 
     logical :: write_newton_details = .true.
     logical :: exit_on_failure      = .false. 

   contains

     ! getters and setters for the above data
     procedure :: get_max_newton_iters, set_max_newton_iters
     procedure :: get_num_vars, set_num_vars

     procedure :: get_atol_unrm, set_atol_unrm
     procedure :: get_atol_rnrm, set_atol_rnrm

     procedure :: get_rtol_unrm, set_rtol_unrm
     procedure :: get_rtol_rnrm, set_rtol_rnrm

     procedure :: set_file_number
     procedure :: set_write_newton_details
     procedure :: set_exit_on_failure

     procedure :: set_init_x, set_init_xdot

  end type newton_solve_bean

contains

  ! max_newton_iters
  integer function get_max_newton_iters(this)

    class(newton_solve_bean) :: this

    get_max_newton_iters = this % max_newton_iters

  end function get_max_newton_iters

  subroutine set_max_newton_iters(this, max_newton_iters)

    class(newton_solve_bean) :: this
    integer :: max_newton_iters

    this % max_newton_iters =  max_newton_iters

  end subroutine set_max_newton_iters

  ! number of variables
  integer function get_num_vars(this)

    class(newton_solve_bean) :: this

    get_num_vars = this % nvars

  end function get_num_vars

  subroutine set_num_vars(this, nvars)

    class(newton_solve_bean) :: this
    integer :: nvars

    this % nvars = nvars

  end subroutine set_num_vars

  ! absolute tolerance of update norm

  real(dp) function get_atol_unrm(this)

    class(newton_solve_bean) :: this

    get_atol_unrm  = this % atol_unrm

  end function get_atol_unrm

  subroutine set_atol_unrm(this, atol_unrm)

    class(newton_solve_bean) :: this
    real(dp) :: atol_unrm

    this % atol_unrm = atol_unrm

  end subroutine set_atol_unrm
  
  ! absolute tolerance of residual norm
  
  real(dp) function get_atol_rnrm(this)

    class(newton_solve_bean) :: this

    get_atol_rnrm  = this % atol_rnrm

  end function get_atol_rnrm

  subroutine set_atol_rnrm(this, atol_rnrm)

    class(newton_solve_bean) :: this
    real(dp) :: atol_rnrm

    this % atol_rnrm = atol_rnrm

  end subroutine set_atol_rnrm

  ! relative tolerance of update norm

  real(dp) function get_rtol_unrm(this)

    class(newton_solve_bean) :: this

    get_rtol_unrm  = this % rtol_unrm

  end function get_rtol_unrm

  subroutine set_rtol_unrm(this, rtol_unrm)

    class(newton_solve_bean) :: this
    real(dp) :: rtol_unrm

    this % rtol_unrm = rtol_unrm

  end subroutine set_rtol_unrm

  ! relative tolerance of residual norm

  real(dp) function get_rtol_rnrm(this)

    class(newton_solve_bean) :: this

    get_rtol_rnrm  = this % rtol_rnrm

  end function get_rtol_rnrm

  subroutine set_rtol_rnrm(this, rtol_rnrm)

    class(newton_solve_bean) :: this
    real(dp) :: rtol_rnrm

    this % rtol_rnrm = rtol_rnrm

  end subroutine set_rtol_rnrm

  ! set the output file number

  subroutine set_file_number(this, filenum)

    class(newton_solve_bean) :: this
    integer :: filenum

    this % filenum = filenum

  end subroutine set_file_number

  ! set the print control for writing the details of newton solve

  subroutine set_write_newton_details(this, write_newton_details)

    class(newton_solve_bean) :: this
    logical :: write_newton_details

    this % write_newton_details = write_newton_details

  end subroutine set_write_newton_details

  ! set whether or not to exit when not converged

  subroutine set_exit_on_failure(this, exit_on_failure)

    class(newton_solve_bean) :: this
    logical :: exit_on_failure

    this % exit_on_failure = exit_on_failure

  end subroutine set_exit_on_failure

  ! set the initial value of x (starting point)

  subroutine set_init_x(this, init_q)

    class(newton_solve_bean) :: this
    real(dp) :: init_q

    this % init_q = init_q

  end subroutine set_init_x

  ! set the initial value of x (starting point)

  subroutine set_init_xdot(this, init_qdot)

    class(newton_solve_bean) :: this
    real(dp) :: init_qdot

    this % init_qdot = init_qdot

  end subroutine set_init_xdot

end module newton_solve_bean_class

!=====================================================================!
! A module that uses Newton's root finding method to approximate the
! solution to equations of the form R(q,qdot,qddot) = 0
!
! The module is generic, which means it can be used to find solutions
! (a) R = 0 (using R and dR values)
! (b) R'= 0 (using R' and d2R values)
!
! Usage : In a main program,
!
! type(newton_solve) :: newton ! create an instance of the solve
!
! call newton%init()
! call newton%solve()
!
! Author :  Komahan Boopathy (komahan@gatech.edu)
! =====================================================================!

module newton_solve_class

  ! import dependencies
  use newton_solve_bean_class

  ! No implicit varaible definitions
  implicit none

  ! All routines and varaibles are private by default
  private

  ! Expose datatypes
  public :: newton_solve

  ! A type that encapsulates the data that is used by the Newton's
  ! method
  type newton_solve
     type(newton_solve_bean) :: bean
   contains
     procedure :: init  => newton_init_
     procedure :: solve => newton_solve_
  end type newton_solve

!!$  ! Interface for the initialization of Newton root finding
!!$  interface init
!!$     subroutine newton_init_(this)
!!$       import newton_solve
!!$       type(newton_solve):: this
!!$     end subroutine newton_init_
!!$  end interface init
!!$
!!$  ! Interface for the initialization of Newton root finding
!!$  interface solve
!!$     subroutine newton_solve_(this)
!!$       import newton_solve
!!$       type(newton_solve) :: this
!!$     end subroutine newton_solve_
!!$  end interface solve

contains

  !-------------------------------------------------------------------!
  !                                                                   
  !-------------------------------------------------------------------!

  subroutine newton_init_(this)
    class(newton_solve) :: this
  end subroutine newton_init_

  !-------------------------------------------------------------------!
  !                                                                   
  !-------------------------------------------------------------------!

  subroutine newton_solve_(this)
    class(newton_solve) :: this
  end subroutine newton_solve_

end module newton_solve_class

! Program to test newton_solve_class module
program test_newton_solve_class

  use newton_solve_class
  implicit none  

  type(newton_solve) :: newton

  !  call newton%init()
  !  call newton%solve()

end program test_newton_solve_class

