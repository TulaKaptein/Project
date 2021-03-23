module ShootingAlgorithmModule
    use MakeGridModule
    use ThreePointSchemeModule
    use NumberKinds
    use integration_module
    use PotentialsModule
    implicit none
    save
    private 
    public :: Shooting

    type doubleArray
        real(KREAL), allocatable :: in(:)
        real(KREAL), allocatable :: out(:)
        real(KREAL), allocatable :: final(:)
    end type
    
contains

subroutine Shooting()
    type (Grid) :: gridCalc
    type (doubleArray), allocatable :: yAlpha(:)
    real(KREAL), allocatable :: eigenVectors(:,:)
    real(KREAL), allocatable :: eigenValues(:)
    real(KREAL) :: norm, correction, lambda
    integer(KINT) :: xMatch, i, ii, alpha

    ! make grid
    call New(gridCalc)
    ! calculate three point scheme
    call ThreePointScheme(gridCalc, eigenVectors, eigenValues)

    ! allocate yAlpha array
    allocate(yAlpha(size(eigenValues)))

! for each value of lambda
do i = 1, size(eigenValues)

    xMatch = gridCalc%numberOfPoints/2
    allocate(yAlpha(i)%out(gridCalc%numberOfPoints))
    yAlpha(i)%out = 0
    allocate(yAlpha(i)%in(gridCalc%numberOfPoints))
    yAlpha(i)%in = 0
    allocate(yAlpha(i)%final(gridCalc%numberOfPoints))
    yAlpha(i)%final = 0

    alpha = i
    lambda = eigenValues(alpha)

    do
        yAlpha(i)%in = 0
        yAlpha(i)%out = 0
        call Outwards(gridCalc, eigenVectors(:,alpha), lambda, yAlpha(i)%out)
        call Inwards(gridCalc, eigenVectors(:,alpha), lambda, yAlpha(i)%in)
        
        ! normalize the solutions
        norm = norm2(yAlpha(i)%out)
        yAlpha(i)%out = yAlpha(i)%out/norm 

        norm = norm2(yAlpha(i)%in)
        yAlpha(i)%in = yAlpha(i)%in/norm
    
        correction = deltaLambda(gridCalc, yAlpha(i), xMatch)
        print *, correction

        lambda = lambda - correction

        ! check for convergence
        if (abs(correction) < 1E-10) then 
            print *, "CONVERGED"
            exit
        endif
    enddo

    ! save the final solutions
    yAlpha(i)%final(:xMatch) = yAlpha(i)%out(0:xMatch)
    yAlpha(i)%final(xMatch:) = yAlpha(i)%in(xMatch:)

    ! normalize the solutions
    norm = norm2(yAlpha(i)%final)
    yAlpha(i)%final = yAlpha(i)%final/norm

enddo

    ! print out the values for all eigenvalues
    open(15, file="yAlpha.txt", action="write")
        ! do ii = 1, size(yAlpha(1)%final)
        !     write(15, *) yAlpha(1:)%final(ii)
        ! enddo
    close(15)


    deallocate(eigenVectors, eigenValues)
    call Delete(gridCalc)

end subroutine

real(KREAL) function deltaLambda(grid1, yAlpha, xMatch)
    type (Grid) :: grid1
    type (doubleArray) :: yAlpha
    integer(KINT) :: xMatch
    real(KREAL) :: yAlphaInDeriv, yAlphaOutDeriv, partA, partB, partC

    ! calculate derivative in
    yAlphaInDeriv = (yAlpha%in(xMatch + 1) - yAlpha%in(xMatch - 1))/(2*grid1%h)
    ! calculate derivative out
    yAlphaOutDeriv = (yAlpha%out(xMatch + 1) - yAlpha%out(xMatch - 1))/(2*grid1%h)

    partA = 0.5*(yAlphaInDeriv/yAlpha%in(xMatch) - yAlphaOutDeriv/yAlpha%out(xMatch))

    partB = 0
    call Newton_cotes((yAlpha%out)**2, grid1%h, 0, xMatch, partB)
    partB = partB/(yAlpha%out(xMatch)**2)

    partC = 0
    call Newton_cotes((yAlpha%in)**2, grid1%h, xMatch, grid1%numberOfPoints, partC)
    partC = partC/(yAlpha%in(xMatch)**2)

    deltaLambda = partA * 1/(partB + partC)

end function

subroutine Outwards(grid1, eigenVectorLambda, lambda, yAlpha)
    type(Grid) :: grid1
    real(KREAL) :: eigenVectorLambda(:)
    real(KREAL) :: lambda
    real(KREAL) :: yAlpha(:)
    integer(KINT) :: xMatch, numberOfPoints, i
    real(KREAL) :: potential = 0

    xMatch = grid1%numberOfPoints/2
    numberOfPoints = grid1%numberOfPoints

    yAlpha(1) = eigenVectorLambda(1)
    yAlpha(2) = eigenVectorLambda(2)

    do i = 3, numberOfPoints
        yAlpha(i) = -yAlpha(i-2) + 2*grid1%h**2*(potential - lambda + 1/(grid1%h)**2)*yAlpha(i-1)
    enddo

end subroutine

subroutine Inwards(grid1, eigenVectorLambda, lambda, yAlpha)
    type(Grid) :: grid1
    real(KREAL) :: eigenVectorLambda(:)
    real(KREAL) :: lambda
    real(KREAL) :: yAlpha(:)
    integer(KINT) :: xMatch, numberOfPoints, i
    real(KREAL) :: potential = 0

    xMatch = grid1%numberOfPoints/2
    numberOfPoints = grid1%numberOfPoints

    yAlpha(numberOfPoints) = eigenVectorLambda(numberOfPoints)
    yAlpha(numberOfPoints - 1) = eigenVectorLambda(numberOfPoints-1)

    do i = numberOfPoints - 2, 1, -1
        yAlpha(i) = -yAlpha(i+2) + 2*grid1%h**2*(potential - lambda + 1/(grid1%h)**2)*yAlpha(i+1)
    enddo

end subroutine

! real(KREAL) function CalcStep(yOld, yCurrent, lambda, h) result(yNew)
!     real(KREAL) :: yOld, yCurrent, lambda, h
!     ! PotentiAL staat er nog niet in!!!
!     yNew = -1*yOld + (h**2)*(lambda + 2/(h**2))*yCurrent
! end function
    
end module ShootingAlgorithmModule