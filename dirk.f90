module runge_kutta
  
  implicit none

  private

  public :: dirk

  ! place the coefficients here
  type dirk
     integer :: dirk_order
     real(8), dimension(:,:),  allocatable :: A
     real(8), dimension(:)  ,  allocatable :: B
     real(8), dimension(:)  ,  allocatable :: C
   contains
     procedure :: init, finalize
  end type dirk

contains
  
  ! Initialize the dirk datatype and construct the tableau
  subroutine init(this, dirk_order)

    class(dirk) :: this
    integer :: dirk_order

    ! set the order of integration
    this % dirk_order = dirk_order

    ! allocate space for the tableau
    allocate(this % A(dirk_order, dirk_order))
    allocate(this % B(dirk_order))    
    allocate(this % C(dirk_order))

    ! put the entries into the tableau

    stop"yet to add entries into the tableau"

  end subroutine init

  ! Deallocate the tableau entries
  subroutine finalize(this)
    
    class(dirk) :: this

    if(allocated(this % A)) deallocate(this % A)
    if(allocated(this % B)) deallocate(this % B)
    if(allocated(this % C)) deallocate(this % C)

  end subroutine finalize


  real(8) function approxq(this, q, qdot, h)
    
    class(dirk) :: this
    real(8) :: q, qdot(:)
    real(8) :: h

    !q = q + h * sum(A(1,:)
    
  end function approxq

end module runge_kutta
