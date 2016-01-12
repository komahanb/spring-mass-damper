program main

  implicit none

  !-------------------------------------------------------------------!
  !  Define constants to manage precision [TUNABLE]
  !-------------------------------------------------------------------!
  integer, parameter :: sp = kind(0.0)    ! single precision
  integer, parameter :: dp = kind(0.0d0)  ! double precision

  !-------------------------------------------------------------------!
  !  Integration constants
  !-------------------------------------------------------------------!  
  integer, parameter :: max_bdf_order = 2 ! maximum BDF integ order
  integer, parameter :: num_time_steps = 100   ! number of time steps
  integer, parameter :: max_newton_iters = 50   ! number of newton iters

  integer, parameter :: max_fo_terms = max_bdf_order*(max_bdf_order+1)/2
  integer, parameter :: max_so_terms = (2*max_bdf_order-1)*(2*max_bdf_order-1) + 1

  !real(dp), parameter :: beta(2*(max_bdf_order+1)-1) = (/ 2.25_dp, -6.0_dp, 5.5_dp, -2.0_dp, 0.25_dp /)
  !alpha(3, 1:4) = (/ 11.0_dp/6.0_dp, -3.0_dp, 1.5_dp, -1.0_dp/3.0_dp /)

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
  real(dp)  :: alpha(max_bdf_order, max_fo_terms) = 0.0_dp
  real(dp)  :: beta(max_bdf_order, max_so_terms) = 0.0_dp

  real(dp) :: u(num_time_steps+1)     = 0.0_dp
  real(dp) :: udot(num_time_steps+1)  = 0.0_dp
  real(dp) :: uddot(num_time_steps+1) = 0.0_dp

  real(dp) :: dR = 0.0_dp  ! jacobian
  real(dp) :: R = 0.0_dp   ! residual
  real(dp) :: du = 0.0_dp  ! state update
  real(dp) :: rnrm = 0.0_dp, unrm  = 0.0_dp 
  real(dp) :: rnrm_tol = 1.0d-10
  real(dp) :: unrm_tol = 1.0d-10

  integer  :: order = 0
  integer  :: i, n ! loop variable

  logical  :: newton_details = .false.

  ! set the BDF coefficeints for first derivative
  alpha(1, 1:2) = (/ 1.0, -1.0 /)
  alpha(2, 1:3) = (/ 1.5_dp, -2.0_dp, 0.5_dp /)

  ! set the BDF coefficient for second derivative
  beta(1, 1:3) = (/ 1.0_dp, -2.0_dp, 1.0_dp /)
  beta(2, 1:5) = (/ 2.25_dp, -6.0_dp, 5.5_dp, -2.0_dp, 0.25_dp /)

  ! set the initial condition (known)
  u(1) = 1.0_dp
  udot(1) = 0.0_dp
  ! rearrange the governing equation and solve for uddot(1)
  uddot(1) = -(C*udot(1) + K*u(1))/M

  ! extrapolate to the first time step
  u(2) = u(1) + udot(1)*dt + uddot(1)*dt2/2.0d0
  ! approximate udot using first order BDF
  udot(2) = (u(2) - u(1))/dt
  ! solve for uddot(2)
  uddot(2) = -(C*udot(2) + K*u(2))/M

  ! get the jacobian
  dR = K + C/dt + M/dt2

  print*, dble(0)*dt, u(1), udot(1), uddot(1), du, R
  print*, dble(1)*dt, u(2), udot(2), uddot(2), du, R

  time: do i = 3, num_time_steps + 1 
     
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

     newton: do n = 1, max_newton_iters
        
        ! find udot and uddot        
        
        ! second order approximation to I derivative
        udot(i) &
             &= alpha(2, 1)*u(i)/dt &
             &+ alpha(2, 2)*u(i-1)/dt &
             &+ alpha(2, 3)*u(i-2)/dt

!!$        if (i .ge. 50000) then
!!$           ! second order approx to II derivative
!!$           uddot(i) &
!!$                &= beta(2, 1)*u(i)/dt2 &
!!$                &+ beta(2, 2)*u(i-1)/dt2 &
!!$                &+ beta(2, 3)*u(i-2)/dt2 &
!!$                &+ beta(2, 4)*u(i-3)/dt2 &
!!$                &+ beta(2, 5)*u(i-4)/dt2
!!$        else
        ! first order approx to II derivative
        uddot(i) &
             &= beta(1, 1)*u(i)/dt2 &
             &+ beta(1, 2)*u(i-1)/dt2 &
             &+ beta(1, 3)*u(i-2)/dt2

!!$        end if

        ! get the residual
        R = M*uddot(i) + C*udot(i) + K*u(i)

        ! find the update
        du = -R/dR

        ! apply the updates
        u(i) = u(i) + du
        
        if (newton_details) print*, u(i), udot(i), uddot(i), du, R
        
        rnrm = sqrt(R*R)
        unrm = sqrt(du*du)

        if ((rnrm .le. rnrm_tol) .and. (unrm .le. unrm_tol)) exit 
        
     end do newton
     
     if (newton_details)  write(*, '(2x, a)')  "-------------------------------------------&
          &------------&
          &-----------------------------------------------------------&
          &--------------"

     print*, dble(i-1)*dt, u(i), udot(i), uddot(i), du, R

  end do time

end program main
