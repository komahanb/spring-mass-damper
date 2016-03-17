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
     procedure :: update_newton_state
     procedure :: setup_linear_system
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
    real(8), intent(in) :: qk(:)
    real(8), allocatable, dimension(:) :: res, dq
    real(8), allocatable, dimension(:,:) :: jac
    integer, allocatable, dimension(:) :: ipiv
    integer :: n, info, size
    logical :: conv = .false.
    
    size = this % num_stages * this % nvars
    
    if (.not.allocated(ipiv)) allocate(ipiv(size))
    if (.not.allocated(res)) allocate(res(size))
    if (.not.allocated(dq)) allocate(dq(size))
    if (.not.allocated(jac)) allocate(jac(size,size))
    
    newton: do n = 1, this % max_newton
       
       ! Get the residual of the function
       call this % compute_stage_residual(qk)
       this % fcnt = this % fcnt + 1

       ! Get the jacobian matrix
       call this % compute_stage_jacobian()
       this % fgcnt = this % fgcnt + 1

       ! setup linear system in lapack format
       call this % setup_linear_system(res, jac)

       ! check stopping
       if (norm2(res) .le. this % tol) then
          conv = .true.
          exit newton
       end if

       ! call lapack to solve the stage values system
       dq = -res
       call DGESV(size, 1, jac, size, IPIV, dq, size, INFO)

       ! check stopping
       if (norm2(dq) .le. this % tol) then
          conv = .true.
          exit newton
       end if

       ! update the solution
       call this % update_newton_state(dq)
       
    end do newton

    ! print warning message if not converged
    if (.not. conv) then
       print '("Newton solve failed : iters = ", i3," |R| = ",E10.3," &
            &|dq| = ",E10.3)',&
            & n, norm2(res), norm2(dq)
       stop
    end if
    
    if (allocated(ipiv)) deallocate(ipiv)
    if (allocated(res)) deallocate(res)
    if (allocated(dq)) deallocate(dq)
    if (allocated(jac)) deallocate(jac)

  end subroutine newton_solve

  !-------------------------------------------------------------------!
  ! Routine that packs the matrix in a form that is used in lapack
  !-------------------------------------------------------------------!

  subroutine setup_linear_system(this, res, jac)

    implicit none

    class(irk) :: this
    
    real(8), intent(inout), dimension(:) :: res
    real(8), intent(inout), dimension(:,:) :: jac
    
    integer :: size
    integer :: i, j
    integer :: istart, iend, cnt, istart2, iend2, cnt2

    size = this % nvars * this % num_stages
    
    cnt = 0 

    !-----------------------------------------------------------------!
    ! convert the residual into vector form
    !-----------------------------------------------------------------!

    do i = 1, this % num_stages
             
       cnt = cnt + 1 

       istart = (cnt-1)*this % nvars + 1 
       iend = (cnt)*this % nvars
       
       res (istart:iend) = this % R (i,:) 
       
    end do

    !-----------------------------------------------------------------!
    ! convert the jacobian into matrix form
    ! 
    ! We have a (nvar x nvar) block stored in each i,j location of this % J
    !-----------------------------------------------------------------!
    
    cnt = 0    
    do i = 1, this % num_stages

       cnt = cnt + 1 
       istart = (cnt-1)*this % nvars + 1 
       iend = (cnt)*this % nvars

       cnt2 = 0
       do j = 1, this % num_stages

          cnt2 = cnt2 + 1
          istart2 = (cnt2-1)*this % nvars + 1 
          iend2 = (cnt2)*this % nvars

          jac(istart:iend,istart2:iend2) = this % J(i,j,:,:)

       end do

    end do
    
  end subroutine setup_linear_system
  
  !-------------------------------------------------------------------!
  ! After the solution of stage equations, we update the states using
  ! this call
  ! -------------------------------------------------------------------!
  
  subroutine update_newton_state(this, sol)
    
    implicit none
    
    class(irk) :: this
    real(8) :: sol(:)
    integer :: i
    integer :: istart, iend, cnt

    if (.not. this % descriptor_form) then

       ! update q(k,i)       

       cnt = 0
       do i = 1, this % num_stages

          cnt = cnt + 1 
          
          istart = (cnt-1)*this % nvars + 1 
          iend = (cnt)*this % nvars
          
          this % Q(i,:) = this % Q(i,:) + sol(istart:iend)

       end do

    else
       
       ! update qdot(k,i)

       cnt = 0
       do i = 1, this % num_stages
          
          cnt = cnt + 1 
          
          istart = (cnt-1)*this % nvars + 1 
          iend = (cnt)*this % nvars
          
          this % QDOT(i,:) = this % QDOT(i,:) + sol(istart:iend)
          
       end do

    end if

  end subroutine update_newton_state
  
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
    real(8) :: qk(:)
    integer :: i, m
    external :: F, R
    
    if (.not. this % descriptor_form) then

       do i = 1, this % num_stages
          
          ! compute qdot
          call F(this % nvars, this % T(i), this % Q(i,:),this % QDOT(i,:))
          
          ! compute the stage residuals
          forall(m = 1 : this % nvars)
             this % R(i,m) = this % Q(i,m) - qk(m) &
                  & - this % h * sum(this % A(i,:)*this % QDOT(:,m))
          end forall

          ! print *, this % num_stages, this % R(i,:)

       end do

    else
       
       do i = 1, this % num_stages
          
          ! compute the stage states for the guessed QDOT
          forall(m = 1 : this % nvars)
             this % Q(i,m) = qk(m) &
                  & + this % h*sum(this % A(i,:)*this % QDOT(:, m))
          end forall

          ! compute the stage residuals
          call R(this % nvars, this % T(i), this % Q(i,:), &
               & this % QDOT(i,:), this % R(i,:))

          ! print *, this % num_stages, this % R(i,:)

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
    external :: DFDQ, DRDQDOT
    
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
!                this % J(i,j) = 1.0d0 - this % h * this % A(i,j) &
!                     &* DFDQ(this % T(j), this % Q(j))
                
                call DFDQ(this % nvars, this % T(j), this % Q(j,:),&
                     & this % h, this %A(i,j), this % J(i,j,:,:))
                
                this % J(i,j,:,:) = 1.0d0 - this % h * this % A(i,j) &
                     &* this % J(i,j,:,:)
                
             else

                ! off diagonal entries

                ! compute only when the coeff is nonzero
                if (this % A(i,j) .ne. 0.0d0) then

!                   this % J(i,j) = - this % h * this % A(i,j) &
!                        &* DFDQ(this % T(j), this % Q(j))

                   call DFDQ(this % nvars, this % T(j), this % Q(j,:),&
                        & this % h, this %A(i,j), this % J(i,j,:,:))

                   this % J(i,j,:,:) = - this % h * this % A(i,j) &
                        &* this % J(i,j,:,:)

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
                ! this % J(i,j) = DRDQDOT(this % T(j), this % Q(j), this % QDOT(j), &
                !     & this % h, this %A(i,j))

                call DRDQDOT(this % nvars, this % T(j), &
                     & this % Q(j,:), this % QDOT(j,:), &
                     & this % J(i,j,:,:))
                
                ! multiply with coeffs
                this % J(i,j,:,:) = 1.0d0 + this % h * this % A(i,i) &
                     &* this % J(i,j,:,:)

             else

                ! off diagonal entries

                ! compute only when the coeff is nonzero
                if (this % A(i,j) .ne. 0.0d0) then

                   call DRDQDOT(this % nvars, this % T(j), &
                        & this % Q(j,:), this % QDOT(j,:), &
                        & this % J(i,j,:,:))

                   ! multiply with coeffs
                   this % J(i,j,:,:) = this % h * this % A(i,j) &
                        &* this % J(i,j,:,:)

                end if ! non-zero

             end if  ! diagonal or not
             
          end do

       end do

    end if

  end subroutine compute_stage_jacobian

end module implicit_runge_kutta












