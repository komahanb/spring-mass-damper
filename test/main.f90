program main

  implicit none
    
  !-------------------------------------------------------------------!
  !  Define constants to manage precision [TUNABLE]
  !-------------------------------------------------------------------!
  integer, parameter :: sp = kind(0.0)    ! single precision
  integer, parameter :: dp = kind(0.0d0)  ! double precision

  integer, parameter :: max_bdf_order = 3

  integer :: i, j ! loop variables

  !-------------------------------------------------------------------!
  !  System parameters
  !-------------------------------------------------------------------!  
  real(dp) :: m = 1.0_dp
  real(dp) :: c = 0.02_dp
  real(dp) :: k = 5.0_dp

  real(dp), parameter :: alpha(max_bdf_order) = (/ 1.5_dp, -2.0_dp, -0.5_dp /)
  real(dp), parameter :: beta(2*max_bdf_order-1) = (/ 2.25_dp, -6.0_dp, 5.5_dp, -2.0_dp, 0.25_dp /)

 
  
end program main
