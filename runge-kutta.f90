
program runge_kutta

  implicit none

  integer, parameter :: N = 100
  real(8), parameter :: h = 1.0d-1
  real(8) :: q(N+1)=0.0d0
  integer :: i

  ! set initial condition
  q(1) = 15.0

  call irk(q,n,h)

  stop

  ! Solve using Explicit RK 4
  call explicit_runge_kutta1(q, n, h)
  write (*, '(2F15.6)', advance='yes') (dble(i-1)*h, &
       &( (-1.0d0 + q(1) + exp(dble(i-1)*h)) - q(i)), i = 1, N+1)

  write(*,*)

!  call explicit_euler(q, n, h)
!  write (*, '(2F15.6)', advance='yes') (dble(i-1)*h, q(i), i = 1, N+1)

contains




  subroutine newton_raphson(order, time, q, b, c, K)

    implicit none

    integer :: order
    real(8) :: K(order, order), b(order), c(order)
    real(8) :: time, q
    integer :: i, max_newton = 20
    real(8) :: res = 0.0d0, jac = 0.0d0
    real(8) :: tmp = 0.0d0, fd = 0.0d0, deltak = 0.0d0, small = 1.0d-8

    newton: do  i = 1, max_newton

       ! find residual
       Res = (k(1,1) - ydot(time(1) + h*b(1), q(1) + h*c(1)*k(1,1)))

       ! test with a FD approx of residual
       tmp = (k(1,1) + small - ydot(time(1) + h*b(1), q(1) + h*c(1)*k(1,1)+small))
       fd = (tmp-res)/small

       ! find jacobian (wrong)
       Jac = 1.0d0 - ydot_k(time(1) + h*b(1), q(1) + h*c(1)*k(1,1))*h*b(1)

       !print*, fd, jac

       !solve linear system
       deltak = -Res/Jac

       print*, res, deltak

       ! check stop
       if(abs(res) .le. 1.0d-12 &
            & .or. abs(deltak).le.1.0d-12) exit newton

       ! apply update
       k(1,1) =  k(1,1) + deltak

    end do newton

    print * , k(1,1), ydot(time(1) + h*b(1) , q(1) + h*c(1)*k(1,1))


  end subroutine newton_raphson

  !-------------------------------------------------------------------!
  ! Explicit Runge Kutta for first order system
  !-------------------------------------------------------------------!
  
  subroutine irk(q, N, h)
    
    implicit none

    integer, parameter :: order = 4
    integer :: N
    integer :: i, j
    real(8) :: h
    real(8) :: q(N+1), time(N+1), err(n+1)
    real(8), dimension(order, N) :: K

    real(8) :: A(order-1,order-1), B(order), C(order)
    real(8) :: deltak(50), res(50), jac(50)
    real(8) :: fd=0.0d0, tmp=0.0d0, small = 1.0d-8

    A = 0.0d0
    B = 0.0d0
    C = 0.0d0

!!$    A(1,1) = 1.0d0/2.0d0
!!$    A(2,2) = 1.0d0/2.0d0
!!$    A(3,3) = 1.0d0
!!$
!!$    B(2) = 1.0d0/2.0d0
!!$    B(3) = 1.0d0/2.0d0    
!!$    B(4) = 1.0d0
!!$
!!$    C(1) = 1.0d0/6.0d0
!!$    C(2) = 1.0d0/3.0d0
!!$    C(3) = 1.0d0/3.0d0    
!!$    C(4) = 1.0d0/6.0d0
    
    A(1,1) = 0.5d0
    B(1) =  0.5d0
    C(1) = 1.0d0

    deltak = 0.0d0

    K(1,1) = 1.0d0

    newton: do  i = i, 50
       
       ! find residual
       Res(i) = (k(1,1) - ydot(time(1) + h*b(1), q(1) + h*c(1)*k(1,1)))
       
       ! test with a FD approx of residual
       tmp = (k(1,1) + small - ydot(time(1) + h*b(1), q(1) + h*c(1)*k(1,1)+small))
       fd = (tmp-res(i))/small

       ! find jacobian (wrong)
       Jac(i) = 1.0d0 - ydot_k(time(1) + h*b(1), q(1) + h*c(1)*k(1,1))*h*b(1)
       
       !print*, fd, jac(i)
       
       !solve linear system
       deltak(i) = -Res(i)/Jac(i)

       print*, res(i), deltak(i)

       ! check stop
       if(abs(res(i)) .le. 1.0d-12 &
            & .or. abs(deltak(i)).le.1.0d-12) exit newton
       
       ! apply update
       k(1,1) =  k(1,1) + deltak(i)
       
    end do newton

    print * , k(1,1), ydot(time(1) + h*b(1) , q(1) + h*c(1)*k(1,1))

    
    ! March in time
    do i = 1, N

       
       ! call non linear solver for obtaining the K's by solving the
       ! linearized non-linear system
       ! call newton(K,time,q)



!!$       do j = 1, order

!!$          if (j .eq. 1) then
!!$             K(j,i) = ydot(time(i), q(i))
!!$          else
!!$             K(j,i) = ydot(time(i) + h*B(j) , q(i) + A(j-1,j-1)*K(j-1,i))
!!$          end if
!!$       end do
!!$
!!$       q(i+1) = q(i)   
!!$       do j = 1, order
!!$          q(i+1) = q(i+1) + h*C(j)*K(j,i)
!!$       end do
!!$
!!$       time(i+1) = time(i) + h
!!$
!!$       ! compute the error in solution
!!$       err(i+1) = q(i+1) - (-1.0d0 + q(1) + exp(time(i+1)))
!!$
    end do

  end subroutine irk

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
       err(i+1) = q(i+1) - (-1.0d0 + q(1) + exp(time(i+1)))

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

    ydot = exp(t)

  end function ydot

  !--------------------------------------------------------------------
  ! First derivative of the state with time
  !--------------------------------------------------------------------

  real function ydot_k(t,q)

    real(8) :: t, q
    
    ydot_k = 0.0d0

  end function ydot_k

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
       !err(i+1) = q(i+1) - (-1.0d0 + q(1) + exp(time(i+1)))

    end do

  end subroutine explicit_runge_kutta

end program runge_kutta





