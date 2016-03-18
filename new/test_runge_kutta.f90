!=====================================================================!
! Main program to test the Runge Kutta Module
!=====================================================================!

program main
  
  use explicit_runge_kutta
  use implicit_runge_kutta
  
  implicit none

  integer :: i, kk

  integer, parameter :: M = 2
  real(8), parameter :: h =  0.1d-1
  real(8), parameter :: tinit = 3.0d0
  real(8), parameter :: tfinal = 5.0d0

  real(8), allocatable, dimension(:,:,:) :: q, qdot

  integer :: nargs, j, N

  logical :: descriptor = .true.
  character(len=5) :: str
  
  type(ERK) :: ERKOBJ
  type(IRK) :: IRKOBJ
  type(DIRK) :: DIRKOBJ

  !-------------------------------------------------------------------!
  ! get the command line arguments
  !-------------------------------------------------------------------!
  
  nargs = command_argument_count()
  do j = 1, nargs
     call get_command_argument(j, str)
     if (str .eq. "false") then
        descriptor = .false.
     end if
  end do
  
  !-------------------------------------------------------------------!
  ! find the number of steps needed for time marching
  !-------------------------------------------------------------------!

  N = int((tfinal - tinit)/h)

  allocate(q(4, N+1,M))
  allocate(qdot(4, N+1,M))
  
  !-------------------------------------------------------------------!
  ! Explicit Runge Kutta
  !-------------------------------------------------------------------!

  q = 0.0d0
  qdot = 0.0d0

  q(:,1,1) = 0.0d0
  q(:,1,2) = 1.0d0
  
  do kk = 1, 3
     
     if (kk.eq.1) open(unit=90, file='erk1.dat')
     if (kk.eq.2) open(unit=90, file='erk2.dat')
     if (kk.eq.3) open(unit=90, file='erk3.dat')
     
     ! Test for each RK 
     ERKOBJ % descriptor_form = descriptor

     call ERKOBJ  % initialize(nvars=M,h=h,tinit=tinit,num_stages=kk)
     call ERKOBJ  % integrate(q(kk,:,:), qdot(kk,:,:), N)
     call ERKOBJ  % finalize()
     
     ! Find the error
     do i = 1, N + 1
        write(90, *)  tinit + dble(i-1)*h, (q(kk,i,j),j=1,M)
     end do

  close(90)

  end do

  !-------------------------------------------------------------------!
  ! Implicit Runge Kutta
  !-------------------------------------------------------------------!

  q = 0.0d0
  qdot = 0.0d0
  
  q(:,1,1) = 0.0d0
  q(:,1,2) = 1.0d0
  
  do kk = 1, 1
     
     if (kk.eq.1) open(unit=90, file='irk1.dat')
     if (kk.eq.2) open(unit=90, file='irk2.dat')
     if (kk.eq.3) open(unit=90, file='irk3.dat')
     
     ! Test for each RK 
     IRKOBJ % descriptor_form = descriptor
     
     call IRKOBJ  % initialize(nvars=M,h=h,tinit=tinit,num_stages=kk)
     call IRKOBJ  % integrate(q(kk,:,:), qdot(kk,:,:), N)
     call IRKOBJ  % finalize()
     
     ! Find the error
     do i = 1, N + 1
        write(90, *)  tinit + dble(i-1)*h, (q(kk,i,j),j=1,M)
     end do

  close(90)

  end do

  !-------------------------------------------------------------------!
  ! Diagonally Implicit Runge Kutta
  !-------------------------------------------------------------------!

  q = 0.0d0
  qdot = 0.0d0
  
  q(:,1,1) = 0.0d0
  q(:,1,2) = 1.0d0
  
  do kk = 1, -3
     
     if (kk.eq.1) open(unit=90, file='dirk1.dat')
     if (kk.eq.2) open(unit=90, file='dirk2.dat')
     if (kk.eq.3) open(unit=90, file='dirk3.dat')
     
     ! Test for each RK 
     DIRKOBJ % descriptor_form = descriptor
     
     call DIRKOBJ  % initialize(nvars=M,h=h,tinit=tinit,num_stages=kk)
     call DIRKOBJ  % integrate(q(kk,:,:), qdot(kk,:,:), N)
     call DIRKOBJ  % finalize()
     
     ! Find the error
     do i = 1, N + 1
        write(90, *)  tinit + dble(i-1)*h, (q(kk,i,j),j=1,M)
     end do

  close(90)

  end do
  
contains

  real(8) function exact(k, tinit, h)
    
    integer, intent(in) ::  k
    real(8), intent(in) ::  h, tinit
    real(8) :: t
    
    t = tinit + dble(k-1)*h
    exact = sin(t)
    ! exact = 2.0d0*exp(1.0d0-t*t) 
    
  end function exact

end program main










