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

     integer :: max_newton = 25
     real(8) :: tol = 1.0d-12
     
   contains

     private

     ! implement/override the abstract class routines
     procedure :: setup_butcher_tableau => ButcherIRK
     procedure :: compute_stage_values => compute_stage_values_implicit

     ! Routines that are common to IRK and its subtypes
     procedure :: newton_solve
     procedure :: compute_stage_residual
     procedure :: compute_stage_jacobian

  end type IRK

  !-------------------------------------------------------------------!  
  ! Diagonally implicit Runge-Kutta
  !-------------------------------------------------------------------!  

  type, extends(IRK) :: DIRK

   contains

     private

     ! implement/override the abstract class routines
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
       stop "DIRK Butcher tableau is not implemented for the requested&
            & order/stages"

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
       stop "IRK Butcher tableau is not implemented for the requested &
            & order/stage"

    end if

  end subroutine ButcherIRK


  !-------------------------------------------------------------------!
  ! Get the stage derivative array for the current step and states 
  !-------------------------------------------------------------------!
  
  subroutine compute_stage_values_implicit(this, k, q)

    class(IRK) :: this
    integer, intent(in) :: k 
    real(8), intent(in), dimension(:,:) :: q
    integer :: j

    ! Find the stage times
    do j = 1, this % num_stages
       this % T(j) = this % time + this % C(j)*this % h
    end do
    
    ! Guess the solution for stage states
    if (.not. this % descriptor_form) then
       this % Q = 1.0d0
    else 
       this % QDOT = 1.0d0 
    end if
    
    ! solve the non linear stage equations using Newton's method for
    ! the actual stage states 

    call this % newton_solve(q(k-1,:))

    ! at this point Q, QDOT, T are finalized--just like ERK
    
  end subroutine compute_stage_values_implicit

  !-------------------------------------------------------------------!
  ! Solve nonlinear stage equations using Newton's method at each time
  ! step.
  !
  ! q_{k,i} = q_{k} + h \sum_{j=1}^s {a_{i}j f(t_{k,j}, q_{k,j})
  ! i = 1,\ldots,s 
  !
  ! This yields $s$ equations and $s$ unknown stage values, $q_{k,i}$,
  ! that are solved using Newton's method at each time step
  ! -------------------------------------------------------------------!

  subroutine newton_solve(this, qk)
    
    class(IRK) :: this
    real(8), intent(in) :: qk
    real(8) :: dQ(this % num_stages)
    integer :: n, jj
    integer :: ipiv(this %num_stages), info
    logical :: conv = .false.
    
    newton: do n = 1, this % max_newton

       ! Get the residual of the function
       call this % compute_stage_residual(qk)
       this % fcnt = this % fcnt + 1

       ! Get the jacobian matrix
       call this % compute_stage_jacobian()
       this % fgcnt = this % fgcnt + 1
       
       ! call lapack to solve the stage values system
       dq = this % R
       call DGESV(this % num_stages, 1, this % J, this % num_stages, &
            & IPIV, dq, this % num_stages, INFO)

       ! check stop (change this to norm2 when implementing
       ! multivariate case)
       if (norm2(this % R) .le. this % tol .or. &
            & norm2(dq) .le. this % tol) then
          conv = .true.
          exit newton
       end if
       
       if (.not. this % descriptor_form) then
          ! update q(k,i)
          this % Q = this % Q - dQ
       else
          ! update qdot(k,i)
          this % QDOT = this % QDOT - dQ
       end if

    end do newton

    if (.not. conv) then
       print '("Newton solve failed : iters = ", i3," |R| = ",E10.3," &
            &|dq| = ",E10.3)',&
            & n, norm2(this % R), norm2(dq)
    end if

    print*, this % time, n, int(dble(this % fcnt)/dble(n)), int(dble(this % fgcnt)/dble(n))

  end subroutine newton_solve
  
  !-------------------------------------------------------------------!
  ! Computes the stage residual for the set stage state Y (comes from
  ! Newton's iteration) and sets into the same instance
  !
  ! R_{i}= q_{k,i} - q_{k} - h \sum_{j=1}^s {a_{ij} f(t_{k,j}, q_{k,j})
  ! i = 1,\ldots,s 
  !
  !-------------------------------------------------------------------!
  
  subroutine compute_stage_residual(this, qk)

    class(IRK) :: this
    real(8) :: qk, tmp
    integer :: i, j
    real(8), external :: F, R

    if (.not. this % descriptor_form) then

       do i = 1, this % num_stages
          
          this % QDOT(i) = F(this % T(i), this % Q(i))
          
          !          this % fcnt = this % fcnt + 1
          
          ! compute the stage residuals
          this % R(i) = this % Q(i) - qk - this % h * sum(this % A(i,:)*this % QDOT(:))
          
       end do

    else
       
       do i = 1, this % num_stages
          
          ! compute the stage states for the guessed QDOT
          this % Q(i) = qk + this % h*sum(this % A(i,:)*this % QDOT)
          
          ! compute the stage residuals
          this % R(i) = R(this % T(i), this % Q(i), this % QDOT(i))

!          this % fcnt = this % fcnt + 1
          
       end do
       
    end if

  end subroutine compute_stage_residual

  !-------------------------------------------------------------------!
  ! Computes the stage jacobian and sets into the same instance
  !          J(i,j) = [ 1 - h A(i,j) DFDQ(T(j),Y(j))]
  !-------------------------------------------------------------------!

  subroutine compute_stage_jacobian(this)
    
    class(IRK) :: this
    integer :: i, j, loop
    real(8), external :: DFDQ, DRDQDOT
    
    if (.not. this % descriptor_form) then

       do i = 1, this % num_stages

          select type (this)
          type is (DIRK)
             loop = i
          type is (IRK)
             loop = this % num_stages
          end select

          do j = 1, loop

             if (i .eq. j) then

                ! compute the diagonal entry
                this % J(i,j) = 1.0d0 - this % h * this % A(i,j) &
                     &* DFDQ(this % T(j), this % Q(j))

                this % fgcnt = this % fgcnt + 1

             else

                ! off diagonal entries

                ! compute only when the coeff is nonzero
                if (this % A(i,j) .ne. 0.0d0) then

                   this % J(i,j) = - this % h * this % A(i,j) &
                        &* DFDQ(this % T(j), this % Q(j))

                   this % fgcnt = this % fgcnt + 1

                end if ! non-zero

             end if  ! diagonal or not

          end do

       end do
       
    else

       do i = 1, this % num_stages

          select type(this)
          type is (DIRK)
             loop = i
          type is (IRK)
             loop = this % num_stages
          end select

          do j = 1, loop

             if (i .eq. j) then

                ! compute the diagonal entry
                this % J(i,j) = DRDQDOT(this % T(j), this % Q(j), this % QDOT(j), &
                     & this % h, this %A(i,j))
                
                this % fgcnt = this % fgcnt + 1
                
             else

                ! off diagonal entries

                ! compute only when the coeff is nonzero
                if (this % A(i,j) .ne. 0.0d0) then

                   this % J(i, j) = DRDQDOT(this % T(j), this % Q(j), this % QDOT(j), &
                        & this % h, this % A(i,i))

                   this % fgcnt = this % fgcnt + 1

                end if ! non-zero

             end if  ! diagonal or not

          end do

       end do

    end if

  end subroutine compute_stage_jacobian

end module implicit_runge_kutta







