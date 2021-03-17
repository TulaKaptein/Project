module ThreePointSchemeModule
    use NumberKinds
    use MakeGridModule
    use Diagonalization
    implicit none
    save
    private 
    public :: ThreePointScheme
    
contains
subroutine ThreePointScheme(grid1)
    type(Grid), intent(in) :: grid1
    real(KREAL), allocatable :: matrixS(:,:), matrixV(:,:), matrixL(:,:)
    real(KREAL), allocatable :: eigenVectors(:,:), eigenValues(:)
    integer(KINT) :: i

    ! allocate and initialize the matrices
    allocate(matrixS(grid1%numberOfPoints, grid1%numberOfPoints))
    allocate(matrixV(grid1%numberOfPoints, grid1%numberOfPoints))
    allocate(matrixL(grid1%numberOfPoints, grid1%numberOfPoints))
    matrixS = 0
    matrixV = 0
    matrixL = 0

    ! fill matrix S 
    forall (i = 1:size(matrixS, 1)) 
        matrixS(i,i) = -2
        matrixS(i, i + 1) = 1
        matrixS(i + 1, i) = 1
    endforall

    ! ! fill matrix V according to the potential chosen
    ! forall (i = 1:size(matrixV,1))
    !   matrixV(i,i) = potential(grid1%gridPoints(i))
    ! end forall

    ! do i = 1, grid1%numberOfPoints
    !     print *, matrixS(i,:)
    ! enddo

    ! calculate matrix L
    matrixL = (-1/(grid1%h)**2)*matrixS + matrixV
    print *, grid1%h

    ! allocate and initialize the eigenvectors and -values
    allocate(eigenVectors(grid1%numberOfPoints, grid1%numberOfPoints))
    allocate(eigenValues(grid1%numberOfPoints))
    eigenVectors = 0
    eigenValues = 0

    ! diagonalize matrix L
    call diagonalize(matrixL, eigenVectors, eigenValues)

    print *
    print *, "eigenValues"
    do i = 1, size(eigenValues)
        print *, eigenValues(i)
    enddo
    print *

    print *, "eigenVectors"
    do i = 1, size(eigenVectors, 1)
        print *, eigenVectors(i, :)
    enddo
    print *

    ! ! deallocate the matrices
    ! deallocate(matrixS(grid1%numberOfPoints, grid1%numberOfPoints))
    ! deallocate(matrixV(grid1%numberOfPoints, grid1%numberOfPoints))
    ! deallocate(matrixL(grid1%numberOfPoints, grid1%numberOfPoints))

end subroutine
    
end module ThreePointSchemeModule