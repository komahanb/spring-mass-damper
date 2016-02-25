!=====================================================================!
! Module that wraps the backward difference integration logic
!=====================================================================!

module backward_difference_class

  use precision
  use unsteady_problem_class
  use steady_solve_class

  implicit none

  private                                                               

  public aa, bb, nonlinear_solver
  public update_states, extrapolate, approximate_derivatives, integrate
  
  type, extends(abstract_integrator) :: backward_difference

     ! Coefficients for linearising the residual
     real(dp) :: alpha(6+1), beta(2*6+1)
     real(dp) :: aa = 0.0_dp, bb = 0.0_dp

   contains

  end type backward_difference

contains

  !-------------------------------------------------------------------!
  ! Returns the bdf coeffs (unscaled with respect to the step size h)
  !-------------------------------------------------------------------!
  ! Input:
  !-------------------------------------------------------------------!
  ! d : d-th derivative
  ! m : order of accuracy
  !-------------------------------------------------------------------!
  ! Output:
  !-------------------------------------------------------------------!
  ! Array of size (d+m) containing the coefficients
  !-------------------------------------------------------------------!

  function get_bdf_coeffs(d, m) result(c)

    integer, intent(in)    :: d
    integer, intent(in)    :: m
    integer                :: n, k

    real(dp)               :: c(d+m), ctmp(d+m)
    real(dp)               :: x(d+m)
    real(dp), parameter    :: h = 1.0_dp

    ! number of data points needed for desired derivative degree and
    ! order of accuracy
    call differ_backward ( h, d, m, ctmp, x )

    n = m + d

    ! flip the array
    forall(k=1:n) c(k) = ctmp(n-k+1)

    if(d .eq. 1) aa = c(1)
    if(d .eq. 2) bb = c(1)

  end function get_bdf_coeffs

  !-------------------------------------------------------------------!
  ! Updates the state vector and its derivatives at the k-th time step
  !-------------------------------------------------------------------!

  subroutine update_states()

    integer :: k
    real(dp) ::dt, dt2

    k = system % get_current_time_step()

    dt = system % get_step_size()
    dt2 = dt * dt

    state % q(k,:) = state % q(k,:) + state % dq(:) 
    state % qdot(k, :) = state % qdot(k,:) + aa*state % dq(:)/dt
    state % qddot(k, :) = state % qddot(k,:) + bb*state % dq(:)/dt2  

  end subroutine update_states

  !-------------------------------------------------------------------!
  ! Approximate first and second derivative for the next iteration
  !-------------------------------------------------------------------!

  subroutine approximate_derivatives()

    implicit none

    real(dp) :: dt, dt2
    integer  :: forder, sorder, k, i

    k = system % get_current_time_step()

    dt = system % get_step_size()
    dt2 = dt * dt

    !-----------------------------------------------------------------!
    ! Approximate state at the current step (copy the previous state)
    !-----------------------------------------------------------------!

    ! state % q(k, :) = state % q(k-1, :)

    !-----------------------------------------------------------------!
    ! FD approximation to I derivative
    !-----------------------------------------------------------------!

    ! Determine the order of accuracy for first derivative
    forder = k - 1
    if (forder .gt. 6) forder = 6

    ! Get the BDF coefficients for first derivative
    alpha(1:forder+1) = get_bdf_coeffs(1, forder)

    ! Apply the BDF formula
    do i = 0, forder ! m+1 points
       state % qdot(k,:) = state % qdot(k,:) + &
            & alpha(i+1) * state % q(k-i,:)/dt 
    end do

    !-----------------------------------------------------------------!
    ! FD approximation to II derivative
    !-----------------------------------------------------------------!

    ! Order of accuracy for second derivative
    sorder = (k-1)/2
    if (sorder .gt. 6) sorder = 6

    if (sorder .gt. 0) then

       ! Get the BDF coefficients for second derivative   
       beta(1:2*sorder+1) = get_bdf_coeffs(2, sorder)

       ! Apply the BDF formula for second derivative
       do i = 0, 2*sorder ! 2m+1 points
          state % qddot(k,:) = state % qddot(k,:) + &
               & beta(i+1) * state % q(k-i,:)/dt2
       end do

    else

       !  we dont have enought points yet
       bb = 1.0_dp
       state % qddot(k,:) = (state%qdot(k-1,:) - state%qdot(k,:))/dt

    end if

  end subroutine approximate_derivatives

  !-------------------------------------------------------------------!
  ! Extrapolate to next time step
  !-------------------------------------------------------------------!

  subroutine extrapolate()

    integer  :: k
    real(dp) :: dt, dt2

    k = system % get_current_time_step()

    dt = system % get_step_size()
    dt2 = dt*dt

    ! Extrapolate using gradient and Hessian information
    if (k .gt. 1) then
       state % q(k,:) = state % q(k-1,:) &
            &+ state % qdot(k-1,:)*dt &
            &+ state % qddot(k-1,:)*dt2/2.0_dp
    end if

  end subroutine extrapolate

  subroutine integrate(qinit, qdotinit)

    real(dp), dimension(:) :: qinit, qdotinit
    integer  :: k
    
    ! set the initial conditions for the problem
    call state % set_initial_state(qinit, qdotinit)

    time: do k = 2, system % get_num_time_steps()

       ! set the current time step
       call system % set_current_time_step(k)

       ! extrapolate from previous time step
       call extrapolate()

       ! approximate the derivatives
       call approximate_derivatives()

       ! solve the linearized system at the current step
       call steady % solve()

       if (.not. steady % is_steady_solve_converged()) exit time

    end do time

  end subroutine integrate

end module backward_difference_class
