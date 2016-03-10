!-------------------------------------------------------------------!
! Function in Explicit form qdot = f(t, q)
!-------------------------------------------------------------------!

real(8) pure function F(time, q)

  real(8), intent(in)  :: time
  real(8), intent(in)  :: q

  !  F = exp(time)

  F = cos(time)

  !  F = tan(q) + 1.0d0 
  !  F = -2.0d0*time*q

  ! F = cos(q) - sin(time)
  ! F = qdot - cos(time)
  ! F = cos(q) - sin(time)
  ! F = exp(time)

end function F

!-------------------------------------------------------------------!
! Function in Implicit form R(t, q, qot) = 0 
!-------------------------------------------------------------------!

real(8) pure function R(time, q, qdot)

  real(8), intent(in)  :: time
  real(8), intent(in)  :: q, qdot

  R = qdot - cos(time)

end function R

!-------------------------------------------------------------------!
! DFDQ of the function
!-------------------------------------------------------------------!

real(8) pure function DFDQ(time, q)

  real(8), intent(in)  :: time
  real(8), intent(in)  :: q

  DFDY = 0.0d0

end function DFDQ

!-------------------------------------------------------------------!
! DRDQDOT of the function
!-------------------------------------------------------------------!

real(8) pure function DRDQDOT(time, q, qdot)

  real(8), intent(in)  :: time
  real(8), intent(in)  :: q, qdot

  DFDQDOT = 1.0d0 
  
end function DRDQDOT




