program runge_kutta

  implicit none

  integer, parameter :: N = 100
  real(8), parameter :: h = 0.1d0
  real(8) :: q(N+1)=0.0d0,qdot(N+1)=0.0d0
  integer :: i

  ! set initial condition
  q(1) = 15.0d0

  ! Solve using IRK
  q = 0.0d0; q(1) = 15.0d0;
  call irk(q, qdot, n,h)
  write (*, '(4F15.6)', advance='yes') (dble(i-1)*h, &
       &(q(i)), (-1.0d0 + q(1) + exp(dble(i-1)*h)),qdot(i), i = 1, N+1)

  ! Solve using ERK4
  !  call explicit_runge_kutta1(q, n, h)
  !  write (*, '(2F15.6)', advance='yes') (dble(i-1)*h, &
  !       &( (-1.0d0 + q(1) + exp(dble(i-1)*h)) - q(i)), i = 1, N+1)
  !  write(*,*)
  
  !  call explicit_euler(q, n, h)
  !  write (*, '(2F15.6)', advance='yes') (dble(i-1)*h, q(i), i = 1, N+1)
  
contains

  !-------------------------------------------------------------------!
  ! Explicit Runge Kutta for first order system
  !-------------------------------------------------------------------!

  subroutine irk(q, qdot, N, h)

    implicit none

    integer, parameter :: order = 1
    integer :: N
    integer :: i, j
    real(8) :: h
    real(8) :: q(N+1), time(N+1),qdot(N+1)
    real(8) :: K(order, N+1)
    real(8) :: A(order,order) = 0.0d0, B(order)=0.0d0, C(order)=0.0d0
    real(8) :: deltak(50), res(50), jac(50)
    real(8) :: fd=0.0d0, tmp=0.0d0, small = 1.0d-8,ktmp

!!$    A(1,1) = 1.0d0/2.0d0
!!$    A(2,2) = 1.0d0/2.0d0
!!$    A(3,3) = 1.0d0
!!$
!!$    B(2) = 1.0d0/2.0d0
!!$    B(3) = 1.0d0/2.0d0    
!!$    B(4) = 1.0d0

!!$    C(1) = 1.0d0/6.0d0
!!$    C(2) = 1.0d0/3.0d0
!!$    C(3) = 1.0d0/3.0d0    
!!$    C(4) = 1.0d0/6.0d0

    K = 0.0d0
    time = 0.0d0

    A(1,1) = 0.5d0
    B(1) = 0.5d0
    C(1) = 1.0d0

    deltak = 0.0d0

    ! March in time
    do i = 2, N + 1

       q(i) = q(i-1)   

       time(i) = time(i-1) + h

       ! Call non linear solver for obtaining the K's by solving the
       ! linearized non-linear system

       ! call newton_raphson(1, time(i), q(i), B, c, qdot(i))
       ! find qdot values at different stages
       call newton1(1, time(I), q(i), qdot(i), b, c)

       K(1,i) = qdot(i)

       ! update the state using the qdot values
       q(i) = q(i-1)   
       do j = 1, order
          q(i) = q(i) + h*C(j)*K(j,i-1)
       end do

       time(i) = time(i-1) + h

    end do

  end subroutine irk

  !-------------------------------------------------------------------!
  ! nonlinear solution at each time step for getting K's
  !-------------------------------------------------------------------!

  subroutine newton_raphson(order, time, q, b, c, K)

    implicit none

    integer :: order
    real(8) :: K, b(order), c(order)
    real(8) :: time, q
    integer :: i, max_newton = 20
    real(8) :: res = 0.0d0, jac = 0.0d0
    real(8) :: tmp = 0.0d0, fd = 0.0d0, deltak = 0.0d0, small = 1.0d-8

    newton: do  i = 1, max_newton

       ! find residual
       Res = (k - ydot(time + h*b(1), q + h*c(1)*k))

       ! test with a FD approx of residual
       tmp = (k + small - ydot( time + h*b(1), q + h*c(1)*(k+small) ))
       fd = (tmp-res)/small

       ! find jacobian
       Jac = 1.0d0 - ydot_k(time + h*b(1), q + h*c(1)*k)*h*b(1)

       print*, "RES:",res, jac, fd

       !solve linear system
       deltak = -Res/Jac

       ! check stop
       if(abs(res) .le. 1.0d-12 &
            & .or. abs(deltak).le.1.0d-12) exit newton

       ! apply update
       k =  k + deltak

    end do newton

    !    print * , k, ydot(time + h*b(1) , q + h*c(1)*k)

  end subroutine newton_raphson

  !-------------------------------------------------------------------!
  ! nonlinear solution at each time step for getting K's
  !-------------------------------------------------------------------!

  subroutine newton1(order, t, q, qdot, b, c)
    
    implicit none

    integer :: order
    real(8) :: K, b(order), c(order)
    real(8) :: t, q, qdot
    integer :: i, max_newton = 20
    real(8) :: res = 0.0d0, jac = 0.0d0
    real(8) :: tmp = 0.0d0, fd = 0.0d0, deltak = 0.0d0, small = 1.0d-8

    newton: do  i = 1, max_newton

       ! find residual
       ! Res = (k - ydot(time + h*b(1), q + h*c(1)*k))

       call residual(qdot, (q+h*c(1)*qdot), (t + h*b(1)), res)

       ! perturb qdot
       call residual(qdot+small, q + h*c(1)*(qdot+small), t+h*b(1), tmp)

       ! call jacobian()
       ! Jac = 1.0d0 - ydot_k(time + h*b(1), q + h*c(1)*k)*h*b(1)

       !       tmp = (k + small - ydot(time + h*b(1), q + h*c(1)*k+small))
       !       fd = (tmp-res)/small
       
       jac = (tmp - res)/small ! dR/dqdot
       
       !print*, res, jac

       ! find jacobian
       !Jac = 1.0d0 - ydot_k(time + h*b(1), q + h*c(1)*k)*h*b(1)
       
       !solve linear system
       deltak = -Res/Jac
       
       ! check stop
       if(abs(res) .le. 1.0d-12 &
            & .or. abs(deltak).le.1.0d-12) exit newton

       ! apply update
       qdot =  qdot + deltak

    end do newton

    print*, "number of  newton: ", i

    !    print * , k, ydot(time + h*b(1) , q + h*c(1)*k)

  end subroutine newton1

  ! residual of the governing equations
  subroutine residual(qdot, q, t, res)

    implicit none

    real(8) :: res
    real(8) :: q, qdot, t

    ! Res = qdot + q
    res = qdot - exp(t)

  end subroutine residual

  ! residual of the governing equations
!!$  subroutine jacobian(qdot,q,t,jac)
!!$    
!!$    implicit none
!!$    
!!$    real(8) :: jac
!!$    real(8) :: q, qdot, t
!!$    
!!$    jac = qdot + q
!!$
!!$  end subroutine jacobian

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

       ! get approximated qdot
       do j = 1, order
          if (j .eq. 1) then
             K(j,i) = ydot(time(i), q(i))
          else
             K(j,i) = ydot(time(i) + h*B(j) , q(i) + A(j-1,j-1)*K(j-1,i))
          end if
       end do

       ! update the states
       q(i+1) = q(i)   
       do j = 1, order
          q(i+1) = q(i+1) + h*C(j)*K(j,i)
       end do

       ! advance in time
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

!!$
!!$    newton: do  i = i, 50
!!$       
!!$       ! find residual
!!$       Res(i) = (k(1,1) - ydot(time(1) + h*b(1), q(1) + h*c(1)*k(1,1)))
!!$       
!!$       ! test with a FD approx of residual
!!$       tmp = (k(1,1) + small - ydot(time(1) + h*b(1), q(1) + h*c(1)*k(1,1)+small))
!!$       fd = (tmp-res(i))/small
!!$
!!$       ! find jacobian (wrong)
!!$       Jac(i) = 1.0d0 - ydot_k(time(1) + h*b(1), q(1) + h*c(1)*k(1,1))*h*b(1)
!!$       
!!$       !print*, fd, jac(i)
!!$       
!!$       !solve linear system
!!$       deltak(i) = -Res(i)/Jac(i)
!!$
!!$       print*, res(i), deltak(i)
!!$
!!$       ! check stop
!!$       if(abs(res(i)) .le. 1.0d-12 &
!!$            & .or. abs(deltak(i)).le.1.0d-12) exit newton
!!$       
!!$       ! apply update
!!$       k(1,1) =  k(1,1) + deltak(i)
!!$       
!!$    end do newton
!!$
!!$    print * , k(1,1), ydot(time(1) + h*b(1) , q(1) + h*c(1)*k(1,1))
