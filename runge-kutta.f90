
program runge_kutta

  implicit none

  integer, parameter :: N = 100
  real(8), parameter :: h = 1.0d-1
  real(8) :: q(N+1)=0.0d0
  integer :: i

  ! set initial condition
  q(1) = 15.0

  ! Solve using Explicit RK 4
  call explicit_runge_kutta1(q, n, h)
  write (*, '(2F15.6)', advance='yes') (dble(i-1)*h, &
       &( (-1.0d0 + q(1) + cos(dble(i-1)*h)) - q(i)), i = 1, N+1)

  write(*,*)

!  call explicit_euler(q, n, h)
!  write (*, '(2F15.6)', advance='yes') (dble(i-1)*h, q(i), i = 1, N+1)

contains

  !-------------------------------------------------------------------!
  ! Explicit Runge Kutta for first order system
  !-------------------------------------------------------------------!

  subroutine explicit_runge_kutta1(q, N, h)

    implicit none

    integer, parameter :: order = 4
    integer :: N
    real(8) :: h
    real(8) :: q(N+1), time(N+1), err(n+1)
    real(8), dimension(order, N) :: K
    integer :: i, j

    real(8) :: A(order-1,order-1), B(order), C(order)

    A = 0.0d0
    B = 0.0d0
    C = 0.0d0

    A(1,1) = 1.0d0/2.0d0
    A(2,2) = 1.0d0/2.0d0
    A(3,3) = 1.0d0

    B(2) = 1.0d0/2.0d0
    B(3) = 1.0d0/2.0d0    
    B(4) = 1.0d0

    C(1) = 1.0d0/6.0d0
    C(2) = 1.0d0/3.0d0
    C(3) = 1.0d0/3.0d0    
    C(4) = 1.0d0/6.0d0

    ! March in time
    do i = 1, N

       do j = 1, order
          if (j .eq. 1) then
             K(j,i) = ydot(time(i), q(i))
          else
             K(j,i) = ydot(time(i) + h*B(j) , q(i) + A(j-1,j-1)*K(j-1,i))
          end if
       end do

       q(i+1) = q(i)   
       do j = 1, order
          q(i+1) = q(i+1) + h*C(j)*K(j,i)
       end do

       time(i+1) = time(i) + h

       ! compute the error in solution
       err(i+1) = q(i+1) - (-1.0d0 + q(1) + cos(time(i+1)))

    end do

  end subroutine explicit_runge_kutta1

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



  !--------------------------------------------------------------------
  ! First derivative of the state with time
  !--------------------------------------------------------------------

  real function fprime(t,q)

    real(8) :: t, q
    
    ydot = cos(t)

  end function fprime

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

       q(i+1) = q(i) + h*(&
            &  K(1,i)/6.0d0 &
            &+ K(2,i)/3.0d0 &
            &+ K(3,i)/3.0d0 &
            &+ K(4,i)/6.0d0 &
            &)

       time(i+1) = time(i) + h

       ! compute the error in solution
       err(i+1) = q(i+1) - (-1.0d0 + q(1) + cos(time(i+1)))

    end do

  end subroutine explicit_runge_kutta

end program runge_kutta





