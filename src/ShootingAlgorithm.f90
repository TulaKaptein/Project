module ShootingAlgorithmModule
    use MakeGridModule
    use ThreePointSchemeModule
    use NumberKinds
    use integration_module
    implicit none
    save
    private 
    public :: Shooting

    type doubleArray
        real(KREAL), allocatable :: in(:)
        real(KREAL), allocatable :: out(:)
    end type
    
contains

subroutine Shooting()
    type (Grid) :: gridCalc
    type (doubleArray) :: yAlpha
    real(KREAL), allocatable :: eigenVectors(:,:)
    real(KREAL), allocatable :: eigenValues(:)
    real(KREAL) :: alpha, norm
    integer(KINT) :: xMatch, i

    ! make grid
    call New(gridCalc)
    ! calculate three point scheme
    call ThreePointScheme(gridCalc, eigenVectors, eigenValues)

    ! for one value of lambda
    xMatch = gridCalc%numberOfPoints/2
    allocate(yAlpha%out(gridCalc%numberOfPoints))
    yAlpha%out = 0
    allocate(yAlpha%in(gridCalc%numberOfPoints))
    yAlpha%in = 0

    alpha = eigenValues(1)

    call Outwards(gridCalc, eigenVectors(:,1), alpha, yAlpha%out)
    call Inwards(gridCalc, eigenVectors(:,1), alpha, yAlpha%in)
    
    ! normalize the solutions
    norm = norm2(yAlpha%out)
    yAlpha%out = yAlpha%out/norm 

    norm = norm2(yAlpha%in)
    yAlpha%in = yAlpha%in/norm

    open(15, file="yAlpha.txt", action="write")
    do i = 1, size(yAlpha%in)
        write(15, *) yAlpha%out(i), yAlpha%in(i)
    enddo
    close(15)
    
    call deltaLambda(gridCalc, yAlpha, xMatch)

    deallocate(eigenVectors, eigenValues)
    call Delete(gridCalc)

end subroutine

subroutine deltaLambda(grid1, yAlpha, xMatch)
    type (Grid) :: grid1
    type (doubleArray) :: yAlpha
    integer(KINT) :: xMatch
    real(KREAL) :: yAlphaInDeriv, yAlphaOutDeriv, partA, partB

    ! calculate derivative in
    yAlphaInDeriv = (yAlpha%in(xMatch + 1) - yAlpha%in(xMatch - 1))/2*grid1%h
    ! calculate derivative out
    yAlphaOutDeriv = (yAlpha%out(xMatch + 1) - yAlpha%out(xMatch - 1))/2*grid1%h
    print *, yAlphaInDeriv, yAlphaOutDeriv, yAlpha%in(xMatch), yAlpha%out(xMatch)

    partA = 0.5*(yAlphaInDeriv/yAlpha%in(xMatch) - yAlphaOutDeriv/yAlpha%out(xMatch))

    ! moet nog ergens een kwadraat in de integraal
    ! call Newton_cotes(yAlpha%out, grid1%h, 0, xMatch, partB)
    ! partB = partB

    print *, partA

end subroutine

subroutine Outwards(grid1, eigenVectorLambda, lambda, yAlpha)
    type(Grid) :: grid1
    real(KREAL) :: eigenVectorLambda(:)
    real(KREAL) :: lambda
    real(KREAL) :: yAlpha(:)
    real(KREAL) :: yOld, yCurrent, yNew
    integer(KINT) :: xMatch, numberOfPoints, i

    yOld = eigenVectorLambda(1)
    yCurrent = eigenVectorLambda(2)
    xMatch = grid1%numberOfPoints/2
    numberOfPoints = grid1%numberOfPoints

    yAlpha(1) = yOld
    yAlpha(2) = yCurrent

    do i = 3, xMatch + 2
        yNew = CalcStep(yOld, yCurrent, lambda, grid1%h)
        yOld = yCurrent
        yCurrent = yNew
        yAlpha(i) = yNew
    enddo

end subroutine

subroutine Inwards(grid1, eigenVectorLambda, lambda, yAlpha)
    type(Grid) :: grid1
    real(KREAL) :: eigenVectorLambda(:)
    real(KREAL) :: lambda
    real(KREAL) :: yAlpha(:)
    real(KREAL) :: yOld, yCurrent, yNew
    integer(KINT) :: xMatch, numberOfPoints, i

    numberOfPoints = grid1%numberOfPoints
    xMatch = grid1%numberOfPoints/2

    yOld = eigenVectorLambda(numberOfPoints)
    yCurrent = eigenVectorLambda(numberOfPoints-1)

    yAlpha(numberOfPoints) = yOld 
    yAlpha(numberOfPoints - 1) = yCurrent

    do i = numberOfPoints - 2, xMatch - 2, -1
        yNew = CalcStep(yOld, yCurrent, lambda, grid1%h)
        yOld = yCurrent
        yCurrent = yNew
        yAlpha(i) = yNew
    enddo

end subroutine

real(KREAL) function CalcStep(yOld, yCurrent, lambda, h) result(yNew)
    real(KREAL) :: yOld, yCurrent, lambda, h
    ! PotentiAL staat er nog niet in!!!
    yNew = -1*yOld + (h**2)*(lambda + 2/(h**2))*yCurrent
end function
    
end module ShootingAlgorithmModule