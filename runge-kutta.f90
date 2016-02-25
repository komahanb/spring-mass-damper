
program runge_kutta

  implicit none

  integer, parameter :: N = 10
  real(8), parameter :: h = 1.0d0
  real(8) :: q(N+1)=0.0d0
  integer :: i

  ! set initial condition
  q(1) = 15.0

  ! Solve using Explicit RK 4
  !call explicit_runge_kutta(q, n, h)
  !write (*, '(2F15.6)', advance='yes') (dble(i-1)*h, q(i), i = 1, N+1)
  !write(*,*)

  call explicit_euler(q, n, h)
  write (*, '(2F15.6)', advance='yes') (dble(i-1)*h, q(i), i = 1, N+1)

  ! Solve using Implicit RK 4

contains
  
  !-------------------------------------------------------------------!
  ! Explicit Runge Kutta for first order system
  !-------------------------------------------------------------------!

  subroutine explicit_runge_kutta(q, N, h)

    implicit none

    integer, parameter :: order = 4
    integer :: N
    real(8) :: h
    real(8) :: q(N+1), time(N+1), err(n+1)
    real(8), dimension(order, N) :: K
    integer :: i

    ! march in time
    do i = 1, N

       K(1,i) = ydot(time(i), q(i))
       K(2,i) = ydot(time(i) + h/2.0d0, q(i) + K(1,i)/2.0d0)
       K(3,i) = ydot(time(i) + h/2.0d0, q(i) + K(2,i)/2.0d0)
       K(4,i) = ydot(time(i) + h, q(i) + K(3,i))

       q(i+1) = q(i) + h*(K(1,i)/6.0d0 + K(2,i)/3.0d0 + K(3,i)/3.0d0 &
            &+ K(4,i)/6.0d0)

       time(i+1) = time(i) + h

       ! compute the error in solution
       err(i+1) = q(i+1) - (-1.0d0 + q(1) + exp(time(i+1)))

    end do

  end subroutine explicit_runge_kutta

  !-------------------------------------------------------------------!
  ! Explicit Runge Kutta for first order system
  !-------------------------------------------------------------------!

  subroutine explicit_euler(q, N, h)

    implicit none

    integer, parameter :: order = 1
    integer :: N
    real(8) :: h
    real(8) :: q(N+1), time(N+1), err(n+1)
    real(8), dimension(order, N) :: K
    integer :: i

    ! march in time
    do i = 1, N

       K(1,i) = ydot(time(i), q(i))

       q(i+1) = q(i) + h*(K(1,i))

       time(i+1) = time(i) + h

    end do

  end subroutine explicit_euler

  !--------------------------------------------------------------------
  ! First derivative of the state with time
  !--------------------------------------------------------------------

  real function ydot(t,q)

    real(8) :: t, q

    ydot = cos(t)

  end function ydot

end program runge_kutta
