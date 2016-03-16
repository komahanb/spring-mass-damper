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
! Function in Implicit form R(t, q, qot) = 0 
!-------------------------------------------------------------------!

real(8) pure function R(time, q, qdot)

  implicit none

  real(8), intent(in)  :: time
  real(8), intent(in)  :: q, qdot

  R = qdot - sin(q) - cos(time)

end function R

!-------------------------------------------------------------------!
! DFDQ of the function
!-------------------------------------------------------------------!

real(8) pure function DFDQ(time, q)

  implicit none

  real(8), intent(in)  :: time
  real(8), intent(in)  :: q

  DFDQ = cos(q)

end function DFDQ

!-------------------------------------------------------------------!
! DRDQDOT of the function
!-------------------------------------------------------------------!

real(8)  function DRDQDOT(time, q, qdot, h, a)

  implicit none

  real(8), intent(in)  :: time
  real(8), intent(in)  :: q, qdot, h,a
  real(8), parameter   :: small = 1.0d-6
  real(8),  external   :: R

  DRDQDOT = 1.0d0 - cos(q)*h*a! (R(time, q, qdot+small) -R(time, q, qdot))/small

  !  print*, DRDQDOT

end function DRDQDOT




