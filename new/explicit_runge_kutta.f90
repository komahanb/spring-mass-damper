!=====================================================================!
! Explicit Runge-Kutta Implementation
!=====================================================================!

module explicit_runge_kutta
  
  use abstract_runge_kutta

  implicit none
  
  private

  public :: ERK
  
  type, extends(RK) :: ERK

   contains

     procedure :: compute_stage_values  => compute_stage_valuesERK
     procedure :: setup_butcher_tableau => ButcherERK

  end type ERK

contains
  
  !-------------------------------------------------------------------!
  ! Butcher's tableau for ERK 
  !-------------------------------------------------------------------!
  
  subroutine ButcherERK(this)

    class(ERK) :: this
    real(8), parameter :: half = 1.0d0/2.0d0
    real(8), parameter :: onethird = 1.0d0/3.0d0
    real(8), parameter :: onesixth = 1.0d0/6.0d0
    real(8), parameter :: oneeight = 1.0d0/8.0d0
    real(8), parameter :: one = 1.0d0
    real(8), parameter :: alpha = 2.0d0/3.0d0

    ! put the entries into the tableau
    if (this % num_stages .eq. 1) then 

       ! Explicit Euler

       this % B(1) = one
       this % C(1) = 0.0d0

       this % order = 1

    else if (this % num_stages .eq. 2) then 

       ! Midpoint method alpha = 1.0d0/2.0d0
       ! Raltson method alpha = 2.0d0/3.0d0
       ! Heun's Method alpha = 1.0d0

       this % A(2,1) = alpha

       this % B(1) = (one-half/alpha)
       this % B(2) = half/alpha

       this % C(1) = 0.0d0
       this % C(2) = alpha

       this % order = 2 

    else if (this % num_stages .eq. 3) then 

       ! Kutta's third order method

       this % A(2,1) = half
       this % A(3,1) = -one
       this % A(3,2)  = 2.0d0

       this % B(1) = onesixth
       this % B(2) = 2.0d0*onethird
       this % B(3) = onesixth

       this % C(1) = 0.0d0
       this % C(2) = half
       this % C(3) = one

       this % order = 3

    else if (this % num_stages .eq. 4) then 

       ! Kutta's formula more accurate (not the classical Runge formula)

       this % A(2,1) = onethird
       this % A(3,1) = -onethird
       this % A(4,1) = one
       this % A(3,2) = one
       this % A(4,2) = -one       
       this % A(4,3) = one

       this % B(1) = oneeight
       this % B(2) = 3.0d0*oneeight
       this % B(3) = this % B(2)
       this % B(4) = oneeight

       this % C(1) = 0.0d0
       this % C(2) = onethird
       this % C(3) = 2.0d0 * this % C(2)
       this % C(4) = one

       this % order = 4

!!$       ! The classical RK tableau
!!$       this % A(2,1) =  half
!!$       this % A(3,2) =  half
!!$       this % A(4,3) =  one
!!$
!!$       this % B(1) = onesixth
!!$       this % B(2) = onethird
!!$       this % B(3) = onethird
!!$       this % B(4) = onesixth
!!$
!!$       this % C(1) = 0.0d0
!!$       this % C(2) = half
!!$       this % C(3) = half
!!$       this % C(4) = one

    else

       print *, this % num_stages
       stop "ERK Butcher tableau is not implemented for the requested order"

    end if

  end subroutine ButcherERK

  !-------------------------------------------------------------------!
  ! Get the stage derivative array for the current step and states ERK
  !-------------------------------------------------------------------!

  subroutine compute_stage_valuesERK(this, k, q)

    class(ERK) :: this
    integer, intent(in) :: k 
    real(8), intent(in), dimension(:) :: q
    real(8) :: tmp
    integer :: j, i 
    real(8), external :: F
    
    ! Stage derivatives are explicitly found at each iteration

    do j = 1, this % num_stages

       ! stage time
       this % T(j) = this % time + this % C(j)*this % h

       ! stage Y
       this % Y(j) = q(k-1) + this % h*sum(this % A(j,:)*this % K(:))

       ! stage derivative
       this % K(j) =  F(this % T(j), this % Y(j))

       ! print*, "j, t, y, F(t,y)",j, this % T(j), this % Y(j), this %K(j)

    end do

  end subroutine compute_stage_valuesERK

end module explicit_runge_kutta
