!=====================================================================!
! Module that implements explicit, implicit and semi-implicit
! Runge-Kutta integration schemes
!=====================================================================!

module runge_kutta_class
  
  implicit none

  private

  public :: DIRK, IRK, ERK

  ! Abstract Runge-Kutta type
  type, abstract :: RK
     ! put all the common variables
   contains
     ! put all the common procedures
     ! procedure :: integrate
  end type RK

  ! Explicit Runge-Kutta
  type, extends(RK) :: ERK
  end type ERK

  ! Implicit Runge-Kutta  
  type, extends(RK) :: IRK
  end type IRK
  
  ! Diagonally implicit Runge-Kutta
  type, extends(RK) :: DIRK

     integer :: dirk_stages
     real(8), dimension(:,:), allocatable :: A ! forms the coeff matrix
     real(8), dimension(:)  , allocatable :: B ! multiplies the state derivatives
     real(8), dimension(:)  , allocatable :: C ! multiplies the time

   contains

     procedure :: init, finalize
     procedure :: get_first_stage_deriv
     procedure :: approx_q
     !procedure :: get_second_stage_deriv
     !procedure :: approx_qdot

  end type DIRK

contains

  !-------------------------------------------------------------------!
  ! Get the stage derivatives by solving the nonlinear system using
  ! Newton's method
  !-------------------------------------------------------------------!

  function get_first_stage_deriv(this) result(ydot)
    
    class(dirk) :: this
    real(8) :: q, qdot(this % dirk_stages)
    real(8) :: h
    integer :: r
    real(8) :: ydot(this % dirk_stages)

    ! solve the nonlinear system to get the stage derivatives
    ! call newton1(1, time(I), q(i), K(1,i), b, c)


  end function get_first_stage_deriv

  ! Approximate q based on the runge-kutta scheme
  real(8) function approx_q(this)
    
    class(dirk) :: this
    real(8) :: q, qdot(this % dirk_stages)
    real(8) :: h
    integer :: r, i, j
    real(8) :: ydot(this % dirk_stages)
    

    ydot = this % get_first_stage_deriv()
    
    do i = 1, this % dirk_stages
       do  j = 1, this % dirk_stages
          q = q + h * this % A(j,i)*ydot(j)
       end do
    end do
    
  end function approx_q

  ! Initialize the dirk datatype and construct the tableau
  subroutine init(this, dirk_stages)

    class(dirk) :: this
    integer :: dirk_stages

    ! set the order of integration
    this % dirk_stages = dirk_stages

    ! allocate space for the tableau
    allocate(this % A(dirk_stages, dirk_stages))
    allocate(this % B(dirk_stages))    
    allocate(this % C(dirk_stages))

    ! put the entries into the tableau
    if (this % dirk_stages .eq. 1) then 
       this % A(1,1) = 0.5d0
       this % B(1)   = 1.0d0
       this % C(1)   = 0.5d0
    else
       stop"yet to add entries into the tableau"
    end if
    
  end subroutine init

  ! Deallocate the tableau entries
  subroutine finalize(this)

    class(dirk) :: this

    if(allocated(this % A)) deallocate(this % A)
    if(allocated(this % B)) deallocate(this % B)
    if(allocated(this % C)) deallocate(this % C)

  end subroutine finalize

end module runge_kutta_class

program main

  use runge_kutta_class

  implicit none

  type(DIRK) :: dirk1, dirk2

  print *, "Beginning execution"
  
  ! Initialize DIRK instance and create the Butcher tableau
  call dirk1 % init(1)
  
  ! call dirk % integrate()

  ! Deallocate the Butcher tableau
  call dirk % finalize()

  print *, "Executing complete"

end program main
