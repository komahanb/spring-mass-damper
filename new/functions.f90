!-------------------------------------------------------------------!
! F of the govenrning equations
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
! DFDY of the function
!-------------------------------------------------------------------!

real(8) pure function DFDY(time, Y)

  real(8), intent(in)  :: time
  real(8), intent(in)  :: Y

  DFDY = 0.0d0

  ! F = qdot + cos(q) - sin(time)
  ! F = qdot - cos(time)
  ! F = cos(q) - sin(time)
  ! F = exp(time)

end function DFDY
