!=====================================================================!
! Main program to test the Runge Kutta Module
!=====================================================================!

program main

  use abstract_runge_kutta
  use explicit_runge_kutta
  use implicit_runge_kutta

  implicit none

  integer :: i, kk,jj

  integer :: N 
  real(8) :: h = 10.0
  real(8), parameter :: tinit = 0.0d0
  real(8), parameter :: tfinal = 100.0d0

  real(8) :: dt(3)

  type(DIRK):: DIRKOBJ
  type(IRK) :: IRKOBJ
  type(ERK) :: ERKOBJ

  real(8), allocatable :: q(:, :), &
       & qdot(:, :), &
       & error(:, :, :)

  ! three different time steps
  dt(1) = 1.0d0
  dt(2) = 1.0d-1
  dt(3) = 1.0d-2

  !-------------------------------------------------------------------!
  ! Explicit Runge Kutta
  !-------------------------------------------------------------------!

  do jj = 1, -2

     h = h/10.0d0 !dt(jj)
     N = int(tfinal-tinit)/h

     allocate (q(4,N+1))
     allocate (qdot(4,N+1))
     allocate (error(4,2,N+1))

     do kk = 1, 4

        q = 0.0d0; qdot =0.0d0;

        q(:,1) = 0.0d0; qdot(:,1)=1.0d0

        ! Test for each RK stage
        call ERKOBJ  % initialize(h=h,tinit=tinit,num_stages=kk)
        call ERKOBJ  % integrate(q(kk,:), qdot(kk,:), N)
        call ERKOBJ  % finalize()

        ! Find the error
        do i = 1, N + 1
           error(kk, 1, i) = abs(q(kk,i) - exact(i, tinit, h))
        end do

        do i = 1, N + 1
           error(kk, 2, i) = abs(qdot(kk,i) - cos(dble(i-1)*ERKOBJ % h))
        end do

     end do ! all stages

     write (*, '(5E15.8)') h, (norm2(error(i, 1, :)), i = 1, 4)

     deallocate(q,qdot,error)

  end do

  !  write (*, '(a,4E15.8)') "ERK :",(norm2(error(i, 2, :)), i = 1, 4)

  !-------------------------------------------------------------------!
  ! Implicit Runge Kutta
  !-------------------------------------------------------------------!

  do jj = 1, 2

     h = h/10.0d0 !dt(jj)
     N = int(tfinal-tinit)/h

     allocate (q(4,N+1))
     allocate (qdot(4,N+1))
     allocate (error(4,2,N+1))

     do kk = 1, 3

        q = 0.0d0; qdot =0.0d0;

        q(:,1) = 0.0d0; qdot(:,1)=1.0d0

        ! Test for each RK stage
        call IRKOBJ  % initialize(h=h,tinit=tinit,num_stages=kk)
        call IRKOBJ  % integrate(q(kk,:), qdot(kk,:), N)
        call IRKOBJ  % finalize()

        ! Find the error
        do i = 1, N + 1
           error(kk, 1, i) = abs(q(kk,i) - exact(i, tinit, h))
        end do

        do i = 1, N + 1
           error(kk, 2, i) = abs(qdot(kk,i) - cos(dble(i-1)*IRKOBJ % h))
        end do

     end do ! all stages

     write (*, '(5E15.8)') h, (norm2(error(i, 1, :)), i = 1, 3)

     deallocate(q,qdot,error)

  end do

  !-------------------------------------------------------------------!
  ! Diagonally Implicit Runge Kutta
  !-------------------------------------------------------------------!

  do jj = 1, -2

     h = h/10.0d0 !dt(jj)
     N = int(tfinal-tinit)/h

     allocate (q(4,N+1))
     allocate (qdot(4,N+1))
     allocate (error(4,2,N+1))

     do kk = 1, 3

        q = 0.0d0; qdot =0.0d0;

        q(:,1) = 0.0d0; qdot(:,1)=1.0d0

        ! Test for each RK stage
        call DIRKOBJ  % initialize(h=h,tinit=tinit,num_stages=kk)
        call DIRKOBJ  % integrate(q(kk,:), qdot(kk,:), N)
        call DIRKOBJ  % finalize()

        ! Find the error
        do i = 1, N + 1
           error(kk, 1, i) = abs(q(kk,i) - exact(i, tinit, h))
        end do

        do i = 1, N + 1
           error(kk, 2, i) = abs(qdot(kk,i) - cos(dble(i-1)*DIRKOBJ % h))
        end do

     end do ! all stages

     write (*, '(5E15.8)') h, (norm2(error(i, 1, :)), i = 1, 3)

     deallocate(q,qdot,error)

  end do

  !stop

contains

  real(8) function exact(k, tinit, h)

    integer, intent(in) ::  k
    real(8), intent(in) ::  h, tinit
    real(8) :: t, c1

    t = tinit + dble(k-1)*h
    exact = sin(t)
    ! exact = 2.0d0*exp(1.0d0-t*t) 

  end function exact

end program main

!!$        print *, tinit + dble(i-1)*h, q(kk,i), &
!!$             & exact(i, tinit, h) , &
!!$             & q(kk,i) - exact(i, tinit, h), &
!!$             & qdot(kk,i) - cos(dble(i-1)*h)
