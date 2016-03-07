!=====================================================================!
! Main program to test the Runge Kutta Module
!=====================================================================!

program main
  
  use abstract_runge_kutta
  use explicit_runge_kutta
  use implicit_runge_kutta
  
  implicit none

  integer :: i, kk

  integer, parameter :: N = 10
  real(8), parameter :: h =  1.0d0
  real(8), parameter :: tinit = 0.0d0

  type(DIRK):: DIRKOBJ
  type(IRK) :: IRKOBJ
  type(ERK) :: ERKOBJ
  
  real(8) :: q(4, N+1) = 0.0d0, &
       & qdot(4, N+1) = 0.0d0, &
       & error(4, 2, N+1) = 0.0d0
  
  !-------------------------------------------------------------------!
  ! Explicit Runge Kutta
  !-------------------------------------------------------------------!

  q = 0.0d0; q(:,1) = 0.0d0; qdot(:,1)=1.0d0
  
  do kk = 1, 4

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
     
  end do

  write (*, '(a,4E15.8)') "ERK :",(norm2(error(i, 1, :)), i = 1, 4)
  write (*, '(a,4E15.8)') "ERK :",(norm2(error(i, 2, :)), i = 1, 4)

  !-------------------------------------------------------------------!
  ! Implicit Runge Kutta
  !-------------------------------------------------------------------!
  
  q = 0.0d0; q(:,1) = 0.0d0; qdot(:,1)=1.0d0

  do kk = 1, 3

     ! Test for each RK stage
     call IRKOBJ % initialize(h=h,tinit=tinit,num_stages=kk)
     call IRKOBJ % integrate(q(kk,:), qdot(kk,:), N)
     call IRKOBJ % finalize()
     
     ! Find the error
     do i = 1, N +1
        error(kk, 1, i) = abs(q(kk,i) - exact(i, tinit, h))
     end do

     do i = 1, N +1
        error(kk, 2, i) = abs(qdot(kk,i) - cos(dble(i-1)*ERKOBJ % h))
     end do

  end do

  write (*, '(a,4E15.8)') "IRK :", (norm2(error(i, 1, :)), i = 1, 3)
  write (*, '(a,4E15.8)') "IRK :", (norm2(error(i, 2, :)), i = 1, 3)

  !-------------------------------------------------------------------!
  ! Diagonally Implicit Runge Kutta
  !-------------------------------------------------------------------!
  
  q = 0.0d0; q(:,1) = 0.0d0; qdot(:,1)=1.0d0
  
  do kk = 1, 3

     ! Test for each RK stage
     call DIRKOBJ  % initialize(h=h,tinit=tinit,num_stages=kk)
     call DIRKOBJ  % integrate(q(kk,:), qdot(kk,:), N)
     call DIRKOBJ  % finalize()

     ! Find the error
     do i = 1, N +1
        error(kk, 1, i) = abs(q(kk,i) - exact(i, tinit, h))
     end do

     do i = 1, N +1
        error(kk, 2, i) = abs(qdot(kk,i) - cos(dble(i-1)*ERKOBJ % h))
     end do

  end do

  write (*, '(a,4E15.8)') "DIRK:", (norm2(error(i, 1, :)), i = 1, 3)
  write (*, '(a,4E15.8)') "DIRK:", (norm2(error(i, 2, :)), i = 1, 3)

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
