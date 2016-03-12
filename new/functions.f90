!-------------------------------------------------------------------!
! Function in Explicit form qdot = f(t, q)
!-------------------------------------------------------------------!

real(8) pure function F(time, q)
  
  implicit none

  real(8), intent(in)  :: time
  real(8), intent(in)  :: q

  !  F = exp(time)
  !  F = cos(time)

  F = sin(q) + cos(time)

  ! F = tan(q) + 1.0d0 
  ! F = -2.0d0*time*q
  ! F = cos(q) - sin(time)
  ! F = qdot - cos(time)
  ! F = cos(q) - sin(time)
  ! F = exp(time)

end function F

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




