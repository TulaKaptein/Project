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

    type tripleArray
        real(KREAL), allocatable :: in(:)
        real(KREAL), allocatable :: out(:)
        real(KREAL), allocatable :: final(:)
    end type
    
contains

subroutine Shooting()
    type(ThreePointSchemeType) :: TPS
    type (tripleArray), allocatable :: yAlpha(:)
    real(KREAL) :: lambda
    integer(KINT) :: xMatch, alpha

    ! start with the three-point scheme
    call NewTPS(TPS, "ParticleInBox")
    call RunTPS(TPS)

    ! allocate yAlpha
    allocate(yAlpha(size(TPS%eigenValues)))

    ! for each eigenvalue/alpha
    do alpha = 1, size(yAlpha)
        ! run the shooting algorithm
        xMatch = TPS%grid%numberOfPoints/2
        call AllocAndInit(yAlpha(alpha)%out, TPS%grid%numberOfPoints)
        call AllocAndInit(yAlpha(alpha)%in, TPS%grid%numberOfPoints)
        call AllocAndInit(yAlpha(alpha)%final, TPS%grid%numberOfPoints)

        lambda = TPS%eigenValues(alpha)

        call ShootingLoopAlpha(TPS%grid, yAlpha(alpha), TPS%eigenVectors(:,alpha),lambda)

        ! save the final solutions
        yAlpha(alpha)%final(:xMatch) = yAlpha(alpha)%out(:xMatch)
        yAlpha(alpha)%final(xMatch + 1:) = yAlpha(alpha)%in(xMatch + 1:)

        ! normalize the solutions
        call Normalize(yAlpha(alpha)%final)

        ! call WriteToFile(yAlpha(alpha)%final, "yAlpha")

    enddo

    call WriteResults(yAlpha)

    ! deallocate everything
    deallocate(yAlpha)
    ! delete three-point scheme ADT
    call DeleteTPS(TPS)
end subroutine

real(KREAL) function deltaLambda(grid, yAlpha, xMatch)
    type (GridType) :: grid
    type (tripleArray) :: yAlpha
    integer(KINT) :: xMatch
    real(KREAL) :: yAlphaInDeriv, yAlphaOutDeriv, partA, partB, partC

    ! calculate derivative in
    yAlphaInDeriv = (yAlpha%in(xMatch + 1) - yAlpha%in(xMatch - 1))/(2*grid%h)
    ! calculate derivative out
    yAlphaOutDeriv = (yAlpha%out(xMatch + 1) - yAlpha%out(xMatch - 1))/(2*grid%h)

    partA = 0.5*(yAlphaInDeriv/yAlpha%in(xMatch) - yAlphaOutDeriv/yAlpha%out(xMatch))

    partB = 0
    call Newton_cotes((yAlpha%out)**2, grid%h, 0, xMatch, partB)
    partB = partB/(yAlpha%out(xMatch)**2)

    partC = 0
    call Newton_cotes((yAlpha%in)**2, grid%h, xMatch, grid%numberOfPoints, partC)
    partC = partC/(yAlpha%in(xMatch)**2)

    deltaLambda = partA * 1/(partB + partC)

end function

subroutine Outwards(grid, eigenVectorLambda, lambda, yAlpha)
    type(GridType) :: grid
    real(KREAL) :: eigenVectorLambda(:)
    real(KREAL) :: lambda
    real(KREAL) :: yAlpha(:)
    integer(KINT) :: xMatch, numberOfPoints, i
    real(KREAL) :: potential = 0

    xMatch = grid%numberOfPoints/2
    numberOfPoints = grid%numberOfPoints

    yAlpha(1) = eigenVectorLambda(1)
    yAlpha(2) = eigenVectorLambda(2)

    do i = 3, xMatch + 1
        yAlpha(i) = -yAlpha(i-2) + 2*grid%h**2*(potential - lambda + 1/(grid%h)**2)*yAlpha(i-1)
    enddo

end subroutine

subroutine Inwards(grid, eigenVectorLambda, lambda, yAlpha)
    type(GridType) :: grid
    real(KREAL) :: eigenVectorLambda(:)
    real(KREAL) :: lambda
    real(KREAL) :: yAlpha(:)
    integer(KINT) :: xMatch, numberOfPoints, i
    real(KREAL) :: potential = 0

    xMatch = grid%numberOfPoints/2
    numberOfPoints = grid%numberOfPoints

    yAlpha(numberOfPoints) = eigenVectorLambda(numberOfPoints)
    yAlpha(numberOfPoints - 1) = eigenVectorLambda(numberOfPoints-1)

    do i = numberOfPoints - 2, xMatch - 1, -1
        yAlpha(i) = -yAlpha(i+2) + 2*grid%h**2*(potential - lambda + 1/(grid%h)**2)*yAlpha(i+1)
    enddo

end subroutine

subroutine ShootingLoopAlpha(grid, y, eigenVector, lambda)
    type(GridType) :: grid
    type(tripleArray) :: y
    real(KREAL) :: eigenVector(:)
    real(KREAL) :: lambda
    real(KREAL) :: correction

    do  
        y%in = 0
        y%out = 0
        call Outwards(grid, eigenVector, lambda, y%out)
        call Inwards(grid, eigenVector, lambda, y%in)

        call Normalize(y%in)
        call Normalize(y%out)

        correction = deltaLambda(grid, y, grid%numberOfPoints/2)
        print *, correction

        lambda = lambda - correction

        if (abs(correction) < 1e-15) then
            print *, "CONVERGED"
            exit
        endif
    enddo
end subroutine

subroutine WriteResults(results)
    type(tripleArray) :: results(:)
    real(KREAL) :: matrixResults(size(results),size(results))
    integer(KINT) :: i

    do i = 1, size(results)
        matrixResults(:,i) = results(i)%final(:)
    enddo

    call WriteToFile(matrixResults, "Results")
end subroutine
    
end module ShootingAlgorithmModule