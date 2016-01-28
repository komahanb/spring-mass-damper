!=====================================================================!
! A module that wraps all the data used in Newton solve
!
! Author :  Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module newton_solve_bean_class

  implicit none

  private
  public :: newton_solve_bean

  type newton_solve_bean     
     logical :: store_data
     logical :: write_data
     real(8) :: atol, rtol
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
! type(newton_solve) :: solve ! create an instance of the solve
!
! call solve%init()
! call solve%solve()
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
     ! all type bound procedures go here

  end type newton_solve

  ! Interface for the initialization of Newton root finding
  interface init
     module procedure newton_init_
  end interface init

  ! Interface for the initialization of Newton root finding
  interface solve
     module procedure newton_solve_
  end interface solve

contains

  !-------------------------------------------------------------------!
  !                                                                   
  !-------------------------------------------------------------------!

  subroutine newton_init_

  end subroutine newton_init_

  !-------------------------------------------------------------------!
  !                                                                   
  !-------------------------------------------------------------------!

  subroutine newton_solve_

  end subroutine newton_solve_

end module newton_solve_class

! Program to test newton_solve_class module
program test_newton_solve_class
  use newton_solve_class
  implicit none  
  type(newton_solve) :: newton
end program test_newton_solve_class
