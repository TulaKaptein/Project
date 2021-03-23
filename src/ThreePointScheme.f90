module ThreePointSchemeModule
    use NumberKinds
    use MakeGridModule
    use Diagonalization
    use PotentialsModule
    implicit none
    save
    private 
    public :: ThreePointScheme
    
contains
subroutine ThreePointScheme(grid1, eigenVectors, eigenValues)
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

    ! fill matrix V
    do i = 1, size(matrixV, 1)
        matrixV(i,i) = ParticleInBox(grid1%gridPoints(i), (grid1%interval(1))*2)
    enddo

    ! calculate matrix L
    matrixL = (-1/(2*(grid1%h)**2))*matrixS + matrixV

    ! allocate and initialize the eigenvectors and -values
    allocate(eigenVectors(grid1%numberOfPoints, grid1%numberOfPoints))
    allocate(eigenValues(grid1%numberOfPoints))
    eigenVectors = 0
    eigenValues = 0

    ! diagonalize matrix L
    call diagonalize(matrixL, eigenVectors, eigenValues)

    ! checks whether the eigenvectors are negative
    if (eigenVectors(1,1) > eigenVectors(1,2)) then
        eigenVectors = -1*eigenVectors
    endif

    open(13, file="eigenValues.txt", action="write")
    write(13, *) "eigenValues"
    do i = 1, size(eigenValues)
        write(13, *) eigenValues(i)
    enddo
    close(13)

    open(14, file="eigenVectors.txt", action="write")
    write(14,*) "eigenVectors"
    do i = 1, size(eigenVectors, 1)
        write(14, *) eigenVectors(i, 1)
    enddo
    close(14)

    ! deallocate the matrices
    deallocate(matrixS, matrixV, matrixL)

end subroutine
    
end module ThreePointSchemeModule