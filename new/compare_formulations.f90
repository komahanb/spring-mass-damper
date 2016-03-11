!=====================================================================!
! Main program to test the Runge Kutta Module
!=====================================================================!

program main
  
  use abstract_runge_kutta
  use explicit_runge_kutta
  use implicit_runge_kutta
  
  implicit none

  integer :: i, kk
  integer :: N 
  real(8), parameter :: h =  0.1d0
  real(8), parameter :: tinit = 0.0d0
  real(8), parameter :: tfinal = 5.0d0
  type(DIRK):: DIRKOBJ
  type(IRK) :: IRKOBJ
  type(ERK) :: ERKOBJ
  integer :: nargs, j
  character(len=5) :: str
  logical :: descriptor = .true.
  real(8), allocatable, dimension(:,:) :: q, qdot

  ! get the command line arguments
  nargs = command_argument_count()
  do j = 1, nargs
     call get_command_argument(j, str)
     if (str .eq. "false") then
        descriptor = .false.
     end if
  end do

  ! find the number of steps needed for time marchine
  N = (tfinal - tinit)/h
  
  allocate(q(4, N+1))
  allocate(qdot(4, N+1))
  
!!$  !-------------------------------------------------------------------!
!!$  ! Explicit Runge Kutta
!!$  !-------------------------------------------------------------------!
!!$
!!$  open(unit=90, file='erk.dat')
!!$
!!$  q = 0.0d0; q(:,1) = 0.0d0; qdot(:,1)=1.0d0
!!$
!!$  do kk = 1, 4
!!$
!!$     ! Test for each RK stage
!!$     
!!$     ERKOBJ % descriptor_form = descriptor
!!$     
!!$     call ERKOBJ  % initialize(h=h,tinit=tinit,num_stages=kk)
!!$     call ERKOBJ  % integrate(q(kk,:), qdot(kk,:), N)
!!$     call ERKOBJ  % finalize()
!!$
!!$     ! Find the error
!!$     do i = 1, N + 1
!!$        write(90, *)  tinit + dble(i-1)*h, q(kk,i)
!!$     end do
!!$
!!$  end do
!!$
!!$  close(90)
!!$  
!!$  !-------------------------------------------------------------------!
!!$  ! Implicit Runge Kutta
!!$  !-------------------------------------------------------------------!
!!$  
!!$  open(unit=90, file='irk.dat')
!!$
!!$  q = 0.0d0; q(:,1) = 0.0d0; qdot(:,1)=1.0d0
!!$
!!$  do kk = 1, 3
!!$
!!$     ! Test for each RK stage
!!$     
!!$     IRKOBJ % descriptor_form = descriptor
!!$
!!$     call IRKOBJ % initialize(h=h,tinit=tinit,num_stages=kk)
!!$     call IRKOBJ % integrate(q(kk,:), qdot(kk,:), N)
!!$     call IRKOBJ % finalize()
!!$
!!$     ! Find the error
!!$     do i = 1, N +1
!!$        write(90, *)  tinit + dble(i-1)*h, q(kk,i)
!!$     end do
!!$
!!$  end do
!!$
!!$  close(90)

  !-------------------------------------------------------------------!
  ! Diagonally Implicit Runge Kutta
  !-------------------------------------------------------------------!
 
  q = 0.0d0; q(:,1) = 0.0d0; qdot(:,1)=1.0d0
  
  do kk = 1, 3
     
     if (kk.eq.1) open(unit=90, file='dirk1.dat')
     if (kk.eq.2) open(unit=90, file='dirk2.dat')
     if (kk.eq.3) open(unit=90, file='dirk3.dat')

     ! Test for each RK stage

     DIRKOBJ % descriptor_form = descriptor

     call DIRKOBJ  % initialize(h=h,tinit=tinit,num_stages=kk)
     call DIRKOBJ  % integrate(q(kk,:), qdot(kk,:), N)
     call DIRKOBJ  % finalize()

     ! Find the error
     do i = 1, N +1
        write(90, *)  tinit + dble(i-1)*h, q(kk,i)
     end do

     close(90)

  end do

  deallocate(q,qdot)

end program main
