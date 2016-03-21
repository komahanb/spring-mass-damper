!-------------------------------------------------------------------!
! Function in Explicit form qdot = f(t, q)
!-------------------------------------------------------------------!

subroutine F(nvars, time, q, qdot)

  implicit none

  integer :: nvars
  real(8) :: time
  real(8) :: q(nvars)
  real(8) :: qdot(nvars)

  if (nvars .eq. 1) then

     qdot = sin(q)

  else

     qdot(1) = q(2)
     qdot(2) = -0.5d0*q(1) + 2.5d0*q(2)

  end if

end subroutine F

!-------------------------------------------------------------------!
! DFDQ of the function
!-------------------------------------------------------------------!

subroutine DFDQ(nvars, time, q, J)

  implicit none

  integer, intent(in) :: nvars  
  real(8), intent(in) :: time
  real(8), intent(in) :: q(nvars)
  real(8), intent(inout) :: J(nvars, nvars)
  
  if (nvars .eq. 1) then

     J(1,1) = J(1,1) + cos(q(1))

  else

     ! derivative of first equation

     J(1,1) = J(1,1) + 0.0d0
     J(1,2) = J(1,2) + 1.0d0

     ! derivative of second equation

     J(2,1) = J(2,1) - 0.5d0
     J(2,2) = J(2,2) + 2.5d0

  end if

end subroutine DFDQ

!-------------------------------------------------------------------!
! DFDQ of the function
!-------------------------------------------------------------------!

subroutine DFDQDOT(nvars, time, q, J)

  implicit none

  integer, intent(in) :: nvars  
  real(8), intent(in) :: time
  real(8), intent(in) :: q(nvars)
  real(8), intent(inout) :: J(nvars, nvars)
  
  if (nvars .eq. 1) then

     J(1,1) = J(1,1) + 1.0d0

  else

     ! derivative of first equation

     J(1,1) = J(1,1) + 1.0d0
     J(1,2) = J(1,2) + 0.0d0

     ! derivative of second equation

     J(2,1) = J(2,1) + 0.0d0
     J(2,2) = J(2,2) + 1.0d0

  end if

end subroutine DFDQDOT

!-------------------------------------------------------------------!
! Function in Implicit form R(t, q, qot) = 0 
!-------------------------------------------------------------------!

subroutine R(nvars, time, q, qdot, res)

  implicit none

  integer, intent(in) :: nvars
  real(8), intent(in) :: time
  real(8), intent(in) :: q(nvars), qdot(nvars)
  real(8), intent(inout) :: res(nvars)

  if (nvars .eq. 1) then

     res(1) = qdot(1) - sin(q(1))

  else

     res(1) = qdot(1) - q(2)
     res(2) = qdot(2) + 0.5d0*q(1) - 2.5d0*q(2)

!!$     res(1) = qdot(1) - q(2) + 2.0d0 * q(1)
!!$     res(2) = qdot(2) + 0.5d0*q(1) + 2.5d0*q(2)

  end if

end subroutine R

!---------------------------------------------------------------------!
! DRDQ of the function
!---------------------------------------------------------------------!

subroutine DRDQ(nvars, time, q, qdot, J)

  implicit none

  integer, intent(in) :: nvars  
  real(8), intent(in) :: time
  real(8), intent(in) :: q(nvars), qdot(nvars)
  real(8), intent(inout) :: J(nvars,nvars)  
  
  if (nvars .eq. 1) then

     ! res(1) = qdot(1) - sin(q(1))
     
     J(1,1) = J(1,1) - cos(q(1))

  else


     !res(1) = qdot(1) - q(2)
     !res(2) = qdot(2) + 0.5d0*q(1) - 2.5d0*q(2)
     
     ! derivative of first equation

     J(1,1) = J(1,1) + 0.0d0
     J(1,2) = J(1,2) - 1.0d0

     ! derivative of second equation

     J(2,1) = J(2,1) + 0.5d0
     J(2,2) = J(2,2) - 2.5d0

  end if

end subroutine DRDQ

!---------------------------------------------------------------------!
! DRDQDOT of the function
!---------------------------------------------------------------------!

subroutine DRDQDOT(nvars, time, q, qdot, J)

  implicit none

  integer, intent(in) :: nvars  
  real(8), intent(in) :: time
  real(8), intent(in) :: q(nvars), qdot(nvars)
  real(8), intent(inout) :: J(nvars,nvars)  

  if (nvars .eq. 1) then

     ! res(1) = qdot(1) - sin(q(1))

     J(1,1) = J(1,1) + 1.0d0
     
  else 

     !res(1) = qdot(1) - q(2)
     !res(2) = qdot(2) + 0.5d0*q(1) - 2.5d0*q(2)

     ! derivative of first equation

     J(1,1) = J(1,1) + 1.0d0
     J(1,2) = J(1,2) + 0.0d0

     ! derivative of second equation

     J(2,1) = J(2,1) + 0.0d0
     J(2,2) = J(2,2) + 1.0d0

  end if

end subroutine DRDQDOT
