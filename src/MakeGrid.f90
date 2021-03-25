module MakeGridModule
    use NumberKinds
    use Tools
    implicit none
    save
    private 
    public :: New, GridType, Delete

    type GridType
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
    type (GridType) :: self
    integer(KINT), intent(in) :: numberOfPoints
    real(KREAL), intent(in) :: interval(2)
    integer (KINT) :: i

    ! save the number of gridpoints and the interval
    self%numberOfPoints = numberOfPoints
    self%interval = interval

    ! allocate the grid points array and initialize
    call AllocAndInit(self%gridPoints, self%numberOfPoints)

    ! calculate h
    self%h = CalculateH(self)

    self%gridPoints(1) = self%interval(1)

    do i = 2, self%numberOfPoints
        self%gridPoints(i) = self%gridpoints(i-1) + (i-1)*self%h
    enddo
end subroutine

subroutine NewPrivateUserInput(self)
    type(GridType) :: self
    integer(KINT) :: i

    ! ask for user input
    print *, "Put in number of gridpoints"
    read *, self%numberOfPoints
    print *, "Put in the interval"
    read *, self%interval(:)

    ! allocate the grid points array and initialize
    call AllocAndInit(self%gridPoints, self%numberOfPoints)

    ! calculate h
    self%h = CalculateH(self)

    self%gridPoints(1) = self%interval(1)

    do i = 2, self%numberOfPoints
        self%gridPoints(i) = self%gridpoints(1) + (i-1)*self%h
    enddo

end subroutine

real(KREAL) function CalculateH(self)
    type(GridType) :: self

    CalculateH = (self%interval(2)-self%interval(1))/self%numberOfPoints
end function

subroutine Delete(self)
    type(GridType) :: self

    deallocate(self%gridPoints)
end subroutine
    
end module MakeGridModule