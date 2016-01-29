!=====================================================================!
! A module that wraps all the data used in Newton solve
!
! Author :  Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module newton_solve_bean_class

  implicit none

  private
  public :: newton_solve_bean

  !-------------------------------------------------------------------!
  !  Define constants to manage precision [TUNABLE]
  !-------------------------------------------------------------------!

  integer, parameter :: sp = kind(0.0)    ! single precision
  integer, parameter :: dp = kind(0.0d0)  ! double precision

  type newton_solve_bean

     !----------------------------------------------------------------!
     ! Basic setup variables
     !----------------------------------------------------------------!

     integer :: max_newton_iters = 50
     integer :: nvars = 1

     !----------------------------------------------------------------!
     ! Solution variables
     !----------------------------------------------------------------!

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

  end type newton_solve_bean

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

  ! Expose routines
  public :: init, solve

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
    type(newton_solve) :: this
  end subroutine newton_init_

  !-------------------------------------------------------------------!
  !                                                                   
  !-------------------------------------------------------------------!

  subroutine newton_solve_(this)
    type(newton_solve) :: this
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
