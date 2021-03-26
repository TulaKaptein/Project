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

subroutine Shooting(potential)
    character(len=*) :: potential
    type(ThreePointSchemeType) :: TPS
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

    ! for each eigenvalue/alpha
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

    call WriteToFile(results, "Results")

    ! deallocate the two arrays of each yAlpha(i)
    do i = 1, size(yAlpha)
        deallocate(yAlpha(i)%in, yAlpha(i)%out)
    enddo

    ! deallocate everything
    deallocate(yAlpha)
    ! delete three-point scheme ADT
    call DeleteTPS(TPS)
end subroutine

subroutine ShootingLoopAlpha(grid, y, eigenVector, lambda, potential)
    type(GridType) :: grid
    type(doubleArray) :: y
    real(KREAL) :: eigenVector(:)
    real(KREAL) :: lambda
    character(len=*) :: potential
    real(KREAL) :: correction

    do  
        y%in = 0
        y%out = 0
        call Outwards(grid, eigenVector, lambda, y%out, potential)
        call Inwards(grid, eigenVector, lambda, y%in, potential)

        call Normalize(y%in)
        call Normalize(y%out)

        correction = deltaLambda(grid, y)
        print *, correction

        lambda = lambda - correction

        if (abs(correction) < 1e-10) then
            print *, "CONVERGED"
            exit
        endif
    enddo
end subroutine
  
subroutine Outwards(grid, eigenVector, lambda, yAlpha, potentialName)
    type(GridType) :: grid
    real(KREAL) :: eigenVector(:)
    real(KREAL) :: lambda
    real(KREAL) :: yAlpha(:)
    character(len=*) :: potentialName
    integer(KINT) :: xMatch, i
    real(KREAL) :: potential, h, l, x, v0, alpha, springConstant

    xMatch = GetDimension(grid)/2
    h = GetH(grid)

    yAlpha(1) = eigenVector(1)
    ! yAlpha(1) = 0
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
        else if (potentialName == "HarmOsc") then
            springConstant = 0.5
            potential = HarmOsc(x, springConstant, l)
        endif    
        yAlpha(i) = -yAlpha(i-2) + 2*h**2*(potential - lambda + 1/h**2)*yAlpha(i-1)  
    enddo

end subroutine

subroutine Inwards(grid, eigenVector, lambda, yAlpha, potentialName)
    type(GridType) :: grid
    real(KREAL) :: eigenVector(:)
    real(KREAL) :: lambda
    real(KREAL) :: yAlpha(:)
    character(len=*) :: potentialName
    integer(KINT) :: xMatch, numberOfPoints, i
    real(KREAL) :: potential, h, l, x, v0, alpha, springConstant

    numberOfPoints = GetDimension(grid)
    xMatch = numberOfPoints/2
    h = GetH(grid)

    yAlpha(numberOfPoints) = eigenVector(numberOfPoints)
    ! yAlpha(numberOfPoints) = 0
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
        else if (potentialName == "HarmOsc") then
            springConstant = 0.5
            potential = HarmOsc(x, springConstant, l)
        endif
        yAlpha(i) = -yAlpha(i+2) + 2*h**2*(potential - lambda + 1/h**2)*yAlpha(i+1)
    enddo

end subroutine

real(KREAL) function deltaLambda(grid, yAlpha)
    type (GridType) :: grid
    type (doubleArray) :: yAlpha
    integer(KINT) :: xMatch, numberOfPoints
    real(KREAL) :: yAlphaInDeriv, yAlphaOutDeriv, partA, partB, partC, h 

    h = GetH(grid)
    numberOfPoints = GetDimension(grid)
    xMatch = numberOfPoints/2

    ! calculate derivative in
    yAlphaInDeriv = (yAlpha%in(xMatch + 1) - yAlpha%in(xMatch - 1))/(2*h)
    ! calculate derivative out
    yAlphaOutDeriv = (yAlpha%out(xMatch + 1) - yAlpha%out(xMatch - 1))/(2*h)

    partA = 0.5*(yAlphaInDeriv/yAlpha%in(xMatch) - yAlphaOutDeriv/yAlpha%out(xMatch))

    partB = 0
    call Newton_cotes((yAlpha%out)**2, h, 0, xMatch, partB)
    partB = partB/(yAlpha%out(xMatch)**2)

    partC = 0
    call Newton_cotes((yAlpha%in)**2, h, xMatch, numberOfPoints, partC)
    partC = partC/(yAlpha%in(xMatch)**2)

    deltaLambda = partA * 1/(partB + partC)

end function

end module ShootingAlgorithmModule