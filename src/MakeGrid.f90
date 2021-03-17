module MakeGridModule
    use NumberKinds
    implicit none
    save
    private 
    public :: New, Grid, Delete

    type Grid
        real(KREAL), allocatable :: gridPoints(:)
        integer(KINT) :: numberOfPoints
        real(KREAL) :: interval(2) ! with interval(1) = lowerbound, interval(2) = upperbound
        real(KREAL) :: h
    end type

    interface New 
        module procedure NewPrivate
        module procedure NewPrivateUserInput
    end interface

    
contains

subroutine NewPrivate(self, numberOfPoints, interval)
    type (Grid) :: self
    integer(KINT), intent(in) :: numberOfPoints
    real(KREAL), intent(in) :: interval(2)
    real (KREAL) :: h
    integer (KINT) :: i

    ! save the number of gridpoints and the interval
    self%numberOfPoints = numberOfPoints
    self%interval = interval

    ! allocate the grid points array and initialize
    allocate(self%gridPoints(numberOfPoints))
    self%gridPoints = 0

    ! calculate h (could be a function!!)
    h = (self%interval(2)-self%interval(1))/self%numberOfPoints

    self%gridPoints(1) = self%interval(1)

    do i = 2, self%numberOfPoints
        self%gridPoints(i) = self%gridpoints(i-1) + (i-1)*h
    enddo
end subroutine

subroutine NewPrivateUserInput(self)
    type(Grid) :: self
    integer(KINT) :: i

    ! ask for user input
    print *, "Put in number of gridpoints"
    read *, self%numberOfPoints
    print *, "Put in the interval"
    read *, self%interval(:)

        ! allocate the grid points array and initialize
    allocate(self%gridPoints(self%numberOfPoints))
    self%gridPoints = 0

    ! calculate h (could be a function!!)
    self%h = (self%interval(2)-self%interval(1))/self%numberOfPoints

    self%gridPoints(1) = self%interval(1)

    do i = 2, self%numberOfPoints
        self%gridPoints(i) = self%gridpoints(1) + (i-1)*self%h
    enddo
end subroutine

subroutine Delete(self)
    type(Grid) :: self

    deallocate(self%gridPoints)
end subroutine
    
end module MakeGridModule