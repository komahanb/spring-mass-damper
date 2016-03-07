!=====================================================================!
! Implicit Runge-Kutta Implementation
!=====================================================================!

module implicit_runge_kutta
  
  use abstract_runge_kutta

  implicit none  

  private

  public :: IRK, DIRK

  !-------------------------------------------------------------------!  
  ! Implicit Runge-Kutta
  !-------------------------------------------------------------------!  
  
  type, extends(RK) :: IRK

   contains

     ! common implicit routines
     procedure :: compute_stage_values => compute_stage_values_implicit
     procedure :: newton_solve
     procedure :: compute_stage_residual
     procedure :: compute_stage_jacobian

     ! routines specialized for this type
     procedure :: setup_butcher_tableau => ButcherIRK

  end type IRK

  !-------------------------------------------------------------------!  
  ! Diagonally implicit Runge-Kutta
  !-------------------------------------------------------------------!  

  type, extends(IRK) :: DIRK

   contains

     procedure :: setup_butcher_tableau => ButcherDIRK

  end type DIRK

contains


  !-------------------------------------------------------------------!
  ! Butcher's tableau for DIRK 
  !-------------------------------------------------------------------!

  subroutine ButcherDIRK(this)

    class(DIRK) :: this
    real(8), parameter :: PI = 22.0d0/7.0d0
    real(8), parameter :: tmp  = 1.0d0/(2.0d0*dsqrt(3.0d0))
    real(8), parameter :: half = 1.0d0/2.0d0
    real(8), parameter :: one  = 1.0d0
    real(8), parameter :: alpha = 2.0d0*cos(PI/18.0d0)/dsqrt(3.0d0)

    ! put the entries into the tableau (ROGER ALEXANDER 1977)
    if (this % num_stages .eq. 1) then 

       ! Implicit mid-point rule (A-stable)

       this % A(1,1) = half
       this % B(1)   = one
       this % C(1)   = half
       
       this % order = 2

!!$       ! Implicit Euler (Backward Euler) but first order accurate
!!$       this % A(1,1) = one
!!$       this % B(1)   = one
!!$       this % C(1)   = one
!!$       this % order = 1
       

    else if (this % num_stages .eq. 2) then

       ! Crouzeix formula (A-stable)

       this % A(1,1) = half + tmp
       this % A(2,1) = -one/dsqrt(3.0d0)
       this % A(2,2) = this % A(1,1)

       this % B(1)   = half
       this % B(2)   = half

       this % C(1)   = half + tmp
       this % C(2)   = half - tmp

       this % order = 3

    else if (this % num_stages .eq. 3) then

       ! Crouzeix formula (A-stable)

       this % A(1,1) = (one+alpha)*half
       this % A(2,1) = -half*alpha
       this % A(3,1) =  one + alpha

       this % A(2,2) = this % A(1,1)
       this % A(3,2) = -(one + 2.0d0*alpha)
       this % A(3,3) = this % A(1,1)
       
       this % B(1)   = one/(6.0d0*alpha*alpha)
       this % B(2)   = one - one/(3.0d0*alpha*alpha)
       this % B(3)   = this % B(1)

       this % C(1) = (one + alpha)*half
       this % C(2) = half
       this % C(3) = (one - alpha)*half
       
       this % order = 4

    else if (this % num_stages .eq. 4) then

       stop "Four stage DIRK formula does not exist"
       
    else
       
       print *, this % num_stages
       stop "DIRK Butcher tableau is not implemented for the requested order"

    end if

  end subroutine ButcherDIRK


  !-------------------------------------------------------------------!
  ! Butcher's tableau for IRK 
  !-------------------------------------------------------------------!

  subroutine ButcherIRK(this)

    class(IRK) :: this
    
    ! put the entries into the tableau
    if (this % num_stages .eq. 1) then 

       ! Implicit mid-point rule (A-stable)

       this % A(1,1) = 0.5d0
       this % B(1)   = 1.0d0
       this % C(1)   = 0.5d0

       this % order = 2
       
!!$       ! Implicit Euler (Backward Euler) but first order accurate
!!$       this % A(1,1) = one
!!$       this % B(1)   = one
!!$       this % C(1)   = one
!!$       this % order = 1
       

    else if (this % num_stages .eq. 2) then 

       ! Radau II A scheme (2 step)

       this % A(1,1) = 5.0d0/12.0d0
       this % A(2,1) = 3.0d0/4.0d0
       this % A(1,2) = -1.0d0/12.0d0
       this % A(2,2) = 1.0d0/4.0d0

       this % B(1) = 3.0d0/4.0d0
       this % B(2) = 1.0d0/4.0d0

       this % C(1) = 1.0d0/3.0d0
       this % C(2) = 1.0d0

       this % order = 3

    else if (this % num_stages .eq. 3) then 

       ! Radau II A scheme (3 step)

       this % A(1,1) = (88.0d0 -7.0d0*dsqrt(6.0d0))/360.0d0
       this % A(2,1) = (296.0d0 + 169.0d0*dsqrt(6.0d0))/1800.0d0
       this % A(3,1) = (16.0d0 - dsqrt(6.0d0))/36.0d0

       this % A(1,2) = (296.0d0 - 169.0d0*dsqrt(6.0d0))/1800.0d0
       this % A(2,2) = (88.0d0 + 7.0d0*dsqrt(6.0d0))/360.0d0
       this % A(3,2) = (16.0d0 + dsqrt(6.0d0))/36.0d0

       this % A(1,3) = (-2.0d0 + 3.0d0*dsqrt(6.0d0))/225.0d0
       this % A(2,3) = (-2.0d0 - 3.0d0*dsqrt(6.0d0))/225.0d0
       this % A(3,3) =  1.0d0/9.0d0

       this % B(1) = (16.0d0 - dsqrt(6.0d0))/36.0d0
       this % B(2) = (16.0d0 + dsqrt(6.0d0))/36.0d0
       this % B(3) = 1.0d0/9.0d0

       this % C(1) = (4.0d0 - dsqrt(6.0d0))/10.0d0
       this % C(2) = (4.0d0 + dsqrt(6.0d0))/10.0d0
       this % C(3) = 1.0d0

       this % order = 5

!!$
!!$       this % A(1,1) = 11.0d0/45.0d0 - 7.0d0*dsqrt(6.0d0)/360.0d0
!!$       this % A(2,1) = 37.0d0/225.0d0 + 169.0d0*dsqrt(6.0d0)/1800.0d0
!!$       this % A(3,1) = 4.0d0/9.0d0 - dsqrt(6.0d0)/36.0d0
!!$
!!$       this % A(1,2) = 37.0d0/225.0d0 - 169.0d0*dsqrt(6.0d0)/1800.0d0
!!$       this % A(2,2) = 11.0d0/45.0d0 + 7.0d0*dsqrt(6.0d0)/360.0d0
!!$       this % A(3,2) = 4.0d0/9.0d0 + dsqrt(6.0d0)/36.0d0
!!$
!!$       this % A(1,3) = -2.0d0/225.0d0 + dsqrt(6.0d0)/75.0d0
!!$       this % A(2,3) = -2.0d0/225.0d0 - dsqrt(6.0d0)/75.0d0 
!!$       this % A(3,3) = 1.0d0/9.0d0
!!$
!!$       this % B(1) = 4.0d0/9.0d0 - dsqrt(6.0d0)/36.0d0
!!$       this % B(2) = 4.0d0/9.0d0 + dsqrt(6.0d0)/36.0d0
!!$       this % B(3) = 1.0d0/9.0d0
!!$
!!$       this % C(1) = 2.0d0/5.0d0 - dsqrt(6.0d0)/10.0d0
!!$       this % C(2) = 2.0d0/5.0d0 + dsqrt(6.0d0)/10.0d0
!!$       this % C(3) = 1.0d0
!!$
!!$       this % order = 4

    else

       print *, this % num_stages
       stop "IRK Butcher tableau is not implemented for the requested order/stage"

    end if

  end subroutine ButcherIRK


  !-------------------------------------------------------------------!
  ! Get the stage derivative array for the current step and states 
  !-------------------------------------------------------------------!
  
  subroutine compute_stage_values_implicit(this, k, q)

    class(IRK) :: this
    integer, intent(in) :: k 
    real(8), intent(in), dimension(:) :: q
    integer :: j

    ! Find the stage times
    do j = 1, this % num_stages
       this % T(j) = dble(k-2)*this % h + this % C(j)*this % h
    end do
    
    ! Guess the solution for stage states
    this % Y(:) = 1.0d0
    
    ! solve the non linear stage equations using Newton's method for
    ! the actual stage states 
    ! S function calls
    ! S jacobian calls
    ! Times the number of newton iterations (N)
    call this % newton_solve(q(k-1))

    ! at this point Y, K, T are finalized--just like ERK
    
  end subroutine compute_stage_values_implicit

  !-------------------------------------------------------------------!
  ! Solve nonlinear stage equations using Newton's method at each time
  ! step
  !-------------------------------------------------------------------!

  subroutine newton_solve(this, qk)
    
    class(IRK) :: this
    real(8), intent(in) :: qk
    real(8) :: dQ(this % num_stages)
    integer :: max_newton = 20
    integer :: n, jj
    integer :: ipiv(this %num_stages), info
    logical :: conv = .false.
    
    newton: do n = 1, max_newton
       
       ! Get the residual of the function
       call this % compute_stage_residual(qk)
       
       ! Get the jacobian matrix
       call this % compute_stage_jacobian()

       ! call lapack to solve the stage values system
       dq = - this % R
       call DGESV(this % num_stages, 1, this % J, this % num_stages, &
            & IPIV, dq, this % num_stages, INFO)

       ! check stop (change this to norm2 when implementing
       ! multivariate case)
       if (norm2(this % R) .le. 1.0d-12 .or. norm2(dq) .le. 1.0d-12) then
          conv = .true.
          exit newton
       end if
       
       ! update q
       this % Y = this % Y + dQ
       
    end do newton

    if (.not.conv) print *, "Newton solve failed after", n

  end subroutine newton_solve
  
  !-------------------------------------------------------------------!
  ! Computes the stage residual and sets into the same instance
  !-------------------------------------------------------------------!
  
  subroutine compute_stage_residual(this, qk)

    class(IRK) :: this
    real(8) :: qk, tmp
    integer :: i, j
    real(8), external :: F

    do i = 1, this % num_stages

       ! this could be the last newton iteration, so we store the
       ! function call value; may be there is better way to handle
       ! this logic

       tmp = 0.0d0

       do j = 1, this % num_stages

          this % K(j) = F(this % T(j), this % Y(j))

          tmp =  tmp + this % A(i,j) * this % K(j)

       end do

       ! yields num_stages equations

       this % R(i) = this % Y(i) - qk - this % h * tmp

    end do

  end subroutine compute_stage_residual

  !-------------------------------------------------------------------!
  ! Computes the stage jacobian and sets into the same instance  
  !          J[s x s] = [I(s)-h A(i,i) DFDY(T(i),Y(i))]
  !-------------------------------------------------------------------!

  subroutine compute_stage_jacobian(this)

    class(IRK) :: this
    integer :: i, j
    real(8), external :: DFDY

    select type (this)

    type is (DIRK)

       ! Jacobian is a lower triangle matrix
       do i = 1, this % num_stages

          do j = 1, i 

             ! Evaluate only when the coeff is nonzero
             if (this % A(i,j) .ne. 0.0d0) then

                if (i .eq. j) then

                   this % J(i,j) = 1.0d0 - this % h * this % A(i,j) &
                        &* DFDY(this % T(j), this % Y(j))

                else

                   this % J(i,j) = - this % h * this % A(i,j) &
                        &* DFDY(this % T(j), this % Y(j))

                end if ! diagonal or not

             end if ! non-zero

          end do

       end do

    type is (IRK)

       ! Jacobian is a FULL  matrix

       do i = 1, this % num_stages

          do j = 1, this % num_stages

             ! Evaluate only when the coeff is nonzero
             if (this % A(i,j) .ne. 0.0d0) then

                if (i .eq. j) then

                   this % J(i,j) = 1.0d0 - this % h * this % A(i,j) &
                        &* DFDY(this % T(j), this % Y(j))

                else

                   this % J(i,j) = - this % h * this % A(i,j) &
                        &* DFDY(this % T(j), this % Y(j))                   

                end if ! diagonal or not

             end if ! non-zero

          end do

       end do

    end select

  end subroutine compute_stage_jacobian

end module implicit_runge_kutta
