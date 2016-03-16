!-------------------------------------------------------------------!
! Function in Explicit form qdot = f(t, q)
!-------------------------------------------------------------------!
subroutine F(nvars,time, q, qdot)

  implicit none

  integer :: nvars
  real(8)  :: time
  real(8) :: q(nvars)
  real(8) :: qdot(nvars)

!  qdot(1) = sin(q(2)) + cos(time)
!  qdot(2) = sin(time) + cos(q(1))
  
  qdot(1) = q(2)
  qdot(2) = -0.5d0*q(1) + 2.5d0*q(2)
  
end subroutine F

!-------------------------------------------------------------------!
! DFDQ of the function
!-------------------------------------------------------------------!
!!$
subroutine DFDQ(nvars, time, q,  h, a, J)

  implicit none

  integer, intent(in) :: nvars  
  real(8), intent(in) :: time, h, a
  real(8), intent(in) :: q(nvars)
  real(8), intent(inout) :: J(nvars, nvars)

  stop "return a nvar nvar block"
  
end subroutine DFDQ

!!$
!!$real(8) pure function DFDQ(time, q)
!!$
!!$  implicit none
!!$
!!$  real(8), intent(in)  :: time
!!$  real(8), intent(in)  :: q
!!$
!!$  DFDQ = cos(q)
!!$
!!$end function DFDQ

!-------------------------------------------------------------------!
! Function in Implicit form R(t, q, qot) = 0 
!-------------------------------------------------------------------!

subroutine R(nvars, time, q, qdot, res)

  implicit none

  integer, intent(in) :: nvars
  real(8), intent(in) :: time
  real(8), intent(in) :: q(nvars), qdot(nvars)
  real(8), intent(inout) :: res(nvars)

  res(1) = qdot(1) - q(2)
  res(2) = qdot(2) + 0.5d0*q(1) - 2.5d0*q(2)

end subroutine R

!---------------------------------------------------------------------!
! DRDQDOT of the function
!---------------------------------------------------------------------!

subroutine DRDQDOT(nvars, time, q, qdot, h, a, J)

  implicit none

  integer, intent(in) :: nvars  
  real(8), intent(in) :: time
  real(8), intent(in) :: q(nvars), qdot(nvars), h,a
  real(8), intent(inout) :: J(nvars,nvars)  
  !real(8), parameter  :: small = 1.0d-6
  !real(8),  external  :: R

  stop"DRDQDOT return nvar nvar block"

end subroutine DRDQDOT

!!$real(8)  function DRDQDOT(time, q, qdot, h, a)
!!$
!!$  implicit none
!!$
!!$  real(8), intent(in)  :: time
!!$  real(8), intent(in)  :: q, qdot, h,a
!!$  real(8), parameter   :: small = 1.0d-6
!!$  real(8),  external   :: R
!!$
!!$  DRDQDOT = 1.0d0 - cos(q)*h*a! (R(time, q, qdot+small) -R(time, q, qdot))/small
!!$
!!$  !  print*, DRDQDOT
!!$
!!$end function DRDQDOT

!!$real(8) pure function R(time, q, qdot)
!!$
!!$  implicit none
!!$
!!$  real(8), intent(in)  :: time
!!$  real(8), intent(in)  :: q, qdot
!!$
!!$  R = qdot - sin(q) - cos(time)
!!$
!!$end function R

