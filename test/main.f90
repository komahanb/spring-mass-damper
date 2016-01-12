!=====================================================================!
! Program that solves a spring-mass-damper system
!
! M uddot + C udot +k u = 0
!
! using Backward Differencing.
!
! Author: Komahan Boopathy 
!=====================================================================!
program spring_mass_damper

  implicit none

  !-------------------------------------------------------------------!
  !  Define constants to manage precision [TUNABLE]
  !-------------------------------------------------------------------!

  integer, parameter :: sp = kind(0.0)    ! single precision
  integer, parameter :: dp = kind(0.0d0)  ! double precision

  !-------------------------------------------------------------------!
  !  Integration constants
  !-------------------------------------------------------------------!  

  integer, parameter :: max_bdf_order = 3       ! maximum BDF integ order
  integer, parameter :: num_time_steps = 1000    ! number of time steps
  integer, parameter :: max_newton_iters = 50   ! number of newton iters

  integer, parameter :: max_fo_terms = max_bdf_order*(max_bdf_order+1)/2
  integer, parameter :: max_so_terms = (2*max_bdf_order-1)*(2*max_bdf_order-1) + 1

  real(dp)  :: alpha(max_bdf_order, max_fo_terms) = 0.0_dp
  real(dp)  :: beta(max_bdf_order, max_so_terms) = 0.0_dp

  real(dp), parameter :: tinit  = 0.0_dp
  real(dp), parameter :: tfinal = 10.0_dp

  real(dp), parameter :: dt = (tfinal-tinit)/dble(num_time_steps)
  real(dp), parameter :: dt2 = dt*dt

  !-------------------------------------------------------------------!
  ! The system parameters 
  !-------------------------------------------------------------------!

  real(dp) :: M = 1.0_dp
  real(dp) :: C = 0.02_dp
  real(dp) :: K = 5.0_dp

  !-------------------------------------------------------------------!
  ! The state variable and derivatives
  !-------------------------------------------------------------------!

  real(dp) :: u(num_time_steps+1)     = 0.0_dp
  real(dp) :: udot(num_time_steps+1)  = 0.0_dp
  real(dp) :: uddot(num_time_steps+1) = 0.0_dp

  real(dp) :: dR = 0.0_dp  ! jacobian
  real(dp) :: R = 0.0_dp   ! residual
  real(dp) :: du = 0.0_dp  ! state update

  real(dp) :: rnrm = 0.0_dp, unrm  = 0.0_dp 
  real(dp) :: rnrm_tol = 1.0d-8
  real(dp) :: unrm_tol = 1.0d-8

  integer  :: order = 0
  integer  :: i, n ! loop variable

  logical  :: newton_details = .false.

  !-------------------------------------------------------------------!
  ! Setup BDF coefficients
  !-------------------------------------------------------------------!

  ! set the BDF coefficeints for first derivative
  alpha(1, 1:2) = (/ 1.0, -1.0 /)
  alpha(2, 1:3) = (/ 1.5_dp, -2.0_dp, 0.5_dp /)
  alpha(3, 1:4) = (/ 11.0_dp/6.0_dp, -3.0_dp, 1.5_dp, -1.0_dp/3.0_dp /)

  ! set the BDF coefficient for second derivative
  beta(1, 1:3) = (/ 1.0_dp, -2.0_dp, 1.0_dp /)
  beta(2, 1:5) = (/ 2.25_dp, -6.0_dp, 5.5_dp, -2.0_dp, 0.25_dp /)

  !-------------------------------------------------------------------!
  ! Set the initial condition (known)
  !-------------------------------------------------------------------!

  u(1) = 1.0_dp
  udot(1) = 0.0_dp
  uddot(1) = -(C*udot(1) + K*u(1))/M   ! rearrange the governing eqn

  print*, dble(0)*dt, 0, u(1), udot(1), uddot(1), du, R, &
       & exact_solution(dble(0)*dt,u(1),udot(1))

!!$  ! extrapolate to the first time step
!!$  u(2) = u(1) + udot(1)*dt + uddot(1)*dt2/2.0d0
!!$  ! approximate udot using first order BDF
!!$  udot(2) = (u(2) - u(1))/dt
!!$  ! solve for uddot(2)
!!$  uddot(2) = -(C*udot(2) + K*u(2))/M

!!$  print*, dble(1)*dt, u(2), udot(2), uddot(2), du, R, &
!!$       &u(2) - exact_solution(dble(0)*dt,u(1),udot(1))
  
  ! get the jacobian (system matrix)
  dR = K + C/dt + M/dt2

  time: do i = 2, num_time_steps + 1 

     ! header
     if (newton_details) then

        write(*, '(2x, a)')  "-------------------------------------------&
             &------------&
             &-----------------------------------------------------------&
             &--------------"
        write(*, '(5x,a,25x,a,25x,a,25x,a,25x,a,25x)')  "u",  "udot",  &
             &"uddot",  "|R|",  "|du|"
        write(*, '(2x, a)')  "-------------------------------------------&
             &------------&
             &-----------------------------------------------------------&
             &--------------"
     end if

     ! extrapolate to the next time step (good starting point for Newton iters)
     u(i) = u(i-1) + udot(i-1)*dt + uddot(i-1)*dt2/2.0_dp

     ! approximate udot using first order BDF
     ! udot(i) = (u(i) - u(i-1))/dt
     ! find uddot(3) algebraically
     ! uddot(i) = -(C*udot(i) + K*u(i))/M

     newton: do n = 2, max_newton_iters

        !-------------------------------------------------------------!
        ! FD approximation to I derivative
        !-------------------------------------------------------------!

        if (i .eq. 2) then
           ! first order approximation to I derivative
           udot(i) &
                &= alpha(1, 1)*u(i)/dt &
                &+ alpha(1, 2)*u(i-1)/dt
        else if (i .eq. 3) then
           ! second order approximation to I derivative
           udot(i) &
                &= alpha(2, 1)*u(i)/dt &
                &+ alpha(2, 2)*u(i-1)/dt &
                &+ alpha(2, 3)*u(i-2)/dt
        else
           ! third order approximation to I derivative
           udot(i) &
                &= alpha(3, 1)*u(i)/dt &
                &+ alpha(3, 2)*u(i-1)/dt &
                &+ alpha(3, 3)*u(i-2)/dt &
                &+ alpha(3, 4)*u(i-3)/dt
        end if

        !-------------------------------------------------------------!
        ! FD approximation to II derivative
        !-------------------------------------------------------------!

!!$        if (i .ge. 50000) then
!!$           ! second order approx to II derivative
!!$           uddot(i) &
!!$                &= beta(2, 1)*u(i)/dt2 &
!!$                &+ beta(2, 2)*u(i-1)/dt2 &
!!$                &+ beta(2, 3)*u(i-2)/dt2 &
!!$                &+ beta(2, 4)*u(i-3)/dt2 &
!!$                &+ beta(2, 5)*u(i-4)/dt2
!!$        else

        if (i .eq. 2) then
           ! first order approx to II derivative (startup)
           uddot(i) = (u(i) - u(i-1))/dt2 - udot(i-1)/dt

        else
           ! first order approx to II derivative
           uddot(i) &
                &= beta(1, 1)*u(i)/dt2 &
                &+ beta(1, 2)*u(i-1)/dt2 &
                &+ beta(1, 3)*u(i-2)/dt2
        end if

        ! get the residual
        R = M*uddot(i) + C*udot(i) + K*u(i)

        ! find the update
        du = -R/dR

        ! apply the updates
        u(i) = u(i) + du

        if (newton_details) print*, u(i), udot(i), uddot(i), du, R

        rnrm = sqrt(R*R)
        unrm = sqrt(du*du)

        if ((rnrm .le. rnrm_tol) .or. (unrm .le. unrm_tol)) exit newton
        
        if (n .eq. max_newton_iters) then
           print*, "res=",r, "du=", du
           stop"Newton failed to converge"
        end if

     end do newton

     if (newton_details)  write(*, '(2x, a)')  "----------------------&
          &---------------------&
          &------------&
          &-----------------------------------------------------------&
          &--------------"

     print*, dble(i-1)*dt, n, u(i), udot(i), udot(i), uddot(i), du, R, &
          & exact_solution(dble(i-1)*dt,u(1),udot(1))

  end do time

contains

  !===================================================================!
  ! Exact solution to the spring mass damper system
  !===================================================================!

  function exact_solution(t, x0, v0) result (x)

    real(dp) :: t, x, x0, v0
    complex(dp) :: mul, a, b, term1, term2, term3

    a = 0.02_dp
    b = 5.0_dp

    mul = exp(-a*t*0.5_dp)/sqrt(a*a - 4.0_dp*b)
    term1 = a*sinh(0.5_dp*t*sqrt(a*a - 4.0_dp*b))
    term2 = sqrt(a*a - 4.0_dp*b)*cosh(0.5_dp*t*sqrt(a*a - 4.0_dp*b))
    term3 = 2.0_dp*v0*sinh(0.5_dp*t*sqrt(a*a - 4.0_dp*b))

    x = real(mul*((term1 + term2)*x0 + term3))

  end function exact_solution

end program spring_mass_damper

