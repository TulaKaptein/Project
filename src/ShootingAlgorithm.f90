!
! Module that uses the shooting algorithm, starting from values calculated from 
! a three-point scheme, to calculate the eigenstates of the 1D Schrödinger equation 
! with a certain potential
!
module ShootingAlgorithmModule
    use MakeGridModule
    use ThreePointSchemeModule
    use NumberKinds
    use integration_module
    use PotentialsModule
    use Tools
    implicit none
    save

    private 
    public :: Shooting

    type doubleArray
        real(KREAL), allocatable :: in(:)
        real(KREAL), allocatable :: out(:)
    end type
    
contains

! the main subroutine
subroutine Shooting(potential)
    character(len=*) :: potential
    type (ThreePointSchemeType) :: TPS
    type (doubleArray), allocatable :: yAlpha(:)
    real(KREAL), allocatable :: results(:,:)
    integer(KINT) :: dimens, xMatch, alpha, i

    ! start with the three-point scheme
    call NewTPS(TPS, potential)
    call RunTPS(TPS)

    dimens = GetDimension(GetGrid(TPS))

    ! allocate yAlpha and results
    allocate(yAlpha(dimens))
    call AllocAndInit(results, dimens)

    ! for each eigenvalue
    do alpha = 1, size(yAlpha)

        ! run the shooting algorithm
        xMatch = dimens/2
        call AllocAndInit(yAlpha(alpha)%out, dimens)
        call AllocAndInit(yAlpha(alpha)%in, dimens)
        call ShootingLoopAlpha(GetGrid(TPS), yAlpha(alpha), GetCertainVector(TPS, alpha), GetCertainValue(TPS, alpha), potential)

        ! save the solutions in a matrix
        results(:xMatch, alpha) = yAlpha(alpha)%out(:xMatch)
        results(xMatch + 1:, alpha) = yAlpha(alpha)%in(xMatch + 1:)

        ! normalize the solutions
        call Normalize(results(:,alpha))

    enddo

    ! write the results to a file
    call WriteToFile(results, "Results")

    ! deallocate the two arrays of each yAlpha(i)
    do i = 1, size(yAlpha)
        deallocate(yAlpha(i)%in, yAlpha(i)%out)
    enddo

    ! deallocate yAlpha
    deallocate(yAlpha)

    ! delete three-point scheme ADT
    call DeleteTPS(TPS)
end subroutine

! the iterative part of the Shooting Algorithm
subroutine ShootingLoopAlpha(grid, y, eigenVector, lambda, potential)
    type(GridType) :: grid
    type(doubleArray) :: y
    real(KREAL) :: eigenVector(:)
    real(KREAL) :: lambda
    character(len=*) :: potential
    real(KREAL) :: correction

    ! start the loop
    do  
        y%in = 0
        y%out = 0

        ! calculate the outwards and inwards integrations of the differential equation
        call Outwards(grid, eigenVector, lambda, y%out, potential)
        call Inwards(grid, eigenVector, lambda, y%in, potential)

        ! normalize the output
        call Normalize(y%in)
        call Normalize(y%out)

        ! correct the eigenvalue
        correction = deltaLambda(grid, y)
        print *, correction
        lambda = lambda - correction

        ! check for convergence
        if (abs(correction) < 1e-15) then
            ! end the loop
            print *, "CONVERGED"
            exit
        endif
    enddo
end subroutine

! integrates the differential equation from the left
subroutine Outwards(grid, eigenVector, lambda, yAlpha, potentialName)
    type(GridType) :: grid
    real(KREAL) :: eigenVector(:)
    real(KREAL) :: lambda
    real(KREAL) :: yAlpha(:)
    character(len=*) :: potentialName
    integer(KINT) :: xMatch, i
    real(KREAL) :: potential, h, l, x, v0, alpha

    ! get xm and h
    xMatch = GetDimension(grid)/2
    h = GetH(grid)

    ! calculate the values of yAlpha from the beginning to xm + 1
    yAlpha(1) = eigenVector(1)
    yAlpha(2) = eigenVector(2)

    do i = 3, xMatch + 1
        x = GetGridPoint(grid, i-1)
        if (potentialName == "ParticleInBox") then
            l = abs(GetLowBound(grid))*2
            potential = ParticleInBox(x, l)
        else if (potentialName == "GaussianPotWell") then
            v0 = 3
            alpha = 0.1
            potential = GaussianPotWell(x, v0, alpha)   
        endif
        yAlpha(i) = -yAlpha(i-2) + 2*h**2*(potential - lambda + 1/h**2)*yAlpha(i-1)  
    enddo

end subroutine

! integrates the differential equation from the right
subroutine Inwards(grid, eigenVector, lambda, yAlpha, potentialName)
    type(GridType) :: grid
    real(KREAL) :: eigenVector(:)
    real(KREAL) :: lambda
    real(KREAL) :: yAlpha(:)
    character(len=*) :: potentialName
    integer(KINT) :: xMatch, numberOfPoints, i
    real(KREAL) :: potential, h, l, x, v0, alpha

    ! get the total number of gridpoints, xm and h
    numberOfPoints = GetDimension(grid)
    xMatch = numberOfPoints/2
    h = GetH(grid)

    ! calculate the values of yAlpha from the end to xm - 1
    yAlpha(numberOfPoints) = eigenVector(numberOfPoints)
    yAlpha(numberOfPoints - 1) = eigenVector(numberOfPoints-1)

    do i = numberOfPoints - 2, xMatch - 1, -1
        x = GetGridPoint(grid, i+1)
        if (potentialName == "ParticleInBox") then
            l = abs(GetLowBound(grid))*2
            potential = ParticleInBox(x, l)
        else if (potentialName == "GaussianPotWell") then
            v0 = 3
            alpha = 0.1
            potential = GaussianPotWell(x, v0, alpha)
        endif
        yAlpha(i) = -yAlpha(i+2) + 2*h**2*(potential - lambda + 1/h**2)*yAlpha(i+1)
    enddo

end subroutine

! calculates the correction to the eigenvalue
real(KREAL) function deltaLambda(grid, yAlpha)
    type (GridType) :: grid
    type (doubleArray) :: yAlpha
    integer(KINT) :: xMatch, numberOfPoints
    real(KREAL) :: dyAlphaIn, dyAlphaOut, partA, partB, partC, h 

    ! get the total number of gridpoints, xm and h
    h = GetH(grid)
    numberOfPoints = GetDimension(grid)
    xMatch = numberOfPoints/2

    ! calculate derivative in and out
    dyAlphaIn = (yAlpha%in(xMatch + 1) - yAlpha%in(xMatch - 1))/(2*h)
    dyAlphaOut = (yAlpha%out(xMatch + 1) - yAlpha%out(xMatch - 1))/(2*h)

    ! calculate part A, B and C of the equation
    partA = 0.5*(dyAlphaIn/yAlpha%in(xMatch) - dyAlphaOut/yAlpha%out(xMatch))
    partB = 0
    call Newton_cotes((yAlpha%out)**2, h, 0, xMatch, partB)
    partB = partB/(yAlpha%out(xMatch)**2)
    partC = 0
    call Newton_cotes((yAlpha%in)**2, h, xMatch, numberOfPoints, partC)
    partC = partC/(yAlpha%in(xMatch)**2)

    ! calculate the correction
    deltaLambda = partA * 1/(partB + partC)

end function

end module ShootingAlgorithmModule