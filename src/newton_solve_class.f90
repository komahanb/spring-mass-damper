!=====================================================================!
! A module that uses Newton's root finding method to approximate the
! solution to equations of the form R(q,qdot,qddot) = 0
!
! The module is generic, which means it can be used to find solutions
! (a) R = 0 (using R and dR values)
! (b) R'= 0 (using R' and d2R values)
!
! Usage : In a main program,
!
! type(newton_solve) :: newton ! create an instance of the solve
!
! call newton % init()
! call newton % solve()
!
! Author :  Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module newton_solve_class

  ! import dependencies
  use precision
  use system_class

  use residual_class
  use jacobian_class

  use steady_solve_class
  use backward_difference, only : update_states

  ! no implicit varaible definitions
  implicit none

  ! all routines and varaibles are private by default
  private

  ! expose datatypes
  public :: newton_solve

  ! a type that contains the logic for Newton's method
  type, extends(steady_solve) :: newton_solve

   contains

     ! public procedures

     !   procedure, public :: solve => solve

     ! private procedures
     procedure :: init => init
     procedure :: finish => finish
     procedure :: work => work

     procedure :: linear_solve => linear_solve
     procedure :: check_stop => check_stop
     procedure :: write_solution => write_solution

  end type newton_solve

contains

  !-------------------------------------------------------------------!
  ! Solve the linear system to find the newton update
  !-------------------------------------------------------------------!
  
  subroutine linear_solve(this)

    class(newton_solve) :: this

    type(residual) :: res
    type(jacobian_matrix) :: jac

    real(dp):: rtmp, drtmp
    integer :: k

    k = system % get_current_time_step()

    if (system % get_num_vars() .eq. 1) then

       rtmp  = state % R(k, 1)
       drtmp = state % dR(k, 1, 1)

       state % dq = -rtmp/drtmp

    else

       !   state % dq    = - res % R(k,:)
       !  / jac % dR (k, :, :)
       stop"linear solve not implemented for ndim>1"
       !  lapack ()

    end if

  end subroutine linear_solve
  
  !-------------------------------------------------------------------!
  ! Check if the desired tolerance has been reached
  !-------------------------------------------------------------------!
  
  subroutine check_stop(this)

    class(newton_solve) :: this

    type(residual) :: res
    integer :: k, n

    k =  system % get_current_time_step()
    n = this % iter_num

    this % rnrm(n) = norm2(state % R(k,:))
    this % unrm(n) = norm2(state % dq)

    if ((this % rnrm(n).le.this % get_atol_rnrm()) .or.&
         & (this % unrm(n).le.this % get_atol_unrm()))  &
         & this % converged = .true.

  end subroutine check_stop

  !-------------------------------------------------------------------!
  ! Write the solution to output file or screen
  !-------------------------------------------------------------------!
  
  subroutine write_solution(this)
    
    class(newton_solve) :: this
    integer :: k, n
    
    if (this % is_write_steady_details()) then
       k = system % get_current_time_step()
       write(this% get_file_number(),*) dble(k)*system%get_step_size(),&
            &this % iter_num, state % q(k,:),  state % qdot(k,:),&
            &state % qddot(k,:),  state % R(k,:), state % dR(k,:,:),&
            &state % dq
    end if

  end subroutine write_solution

  !-------------------------------------------------------------------!
  ! Routine for initialization tasks
  !-------------------------------------------------------------------!

  subroutine init(this)

    class(newton_solve) :: this

    this % converged = .false.

    this % iter_num = 1

    if (allocated(this % rnrm)) deallocate(this % rnrm)
    allocate(this % rnrm(this % get_max_steady_iters()))
    this%rnrm = 0.0_dp

    if (allocated(this % unrm)) deallocate(this % unrm)
    allocate(this % unrm(this % get_max_steady_iters()))
    this%unrm = 0.0_dp

    ! Open a file for dumping newton iteration output
    if (this % is_write_steady_details()) then
       if(this % get_file_number() .ne. 6) then
          open(unit = this% get_file_number(), file="NewtonHistory.dat")
       end if
    end if

  end subroutine init
  
  !-------------------------------------------------------------------!
  ! Routine that wraps the logic of newton solve                                                          
  !-------------------------------------------------------------------!
  
  subroutine work(this)

    class(newton_solve) :: this

    type(jacobian_matrix) :: jac
    type(residual) :: res
    integer :: n

    newton: do n = 1, this % get_max_steady_iters()

       call res % get_residual()

       call jac % get_jacobian()

       call this % linear_solve()

       call this % write_solution()
       
       call this % check_stop()

       if (this % converged) then
          exit newton
       else
          call  update_states()
       end if

       this % iter_num = this % iter_num + 1

    end do newton

  end subroutine work
  
  !-------------------------------------------------------------------!
  ! Routine that performs finishing tasks
  !-------------------------------------------------------------------!
  
  subroutine finish(this)

    class(newton_solve) :: this

    if (allocated(this%rnrm)) deallocate(this%rnrm)
    if (allocated(this%unrm)) deallocate(this%unrm)
    
    ! close the newton iteration output file if opened
    if (this % is_write_steady_details()) then
       if(this% get_file_number() .ne. 6) then
          close(this% get_file_number())
       end if
    end if

  end subroutine finish

end module newton_solve_class
