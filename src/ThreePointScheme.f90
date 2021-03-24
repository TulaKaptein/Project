module ThreePointSchemeModule
    use NumberKinds
    use MakeGridModule
    use Diagonalization
    use PotentialsModule
    use Tools
    implicit none
    save
    private 
    public :: NewTPS, RunTPS, DeleteTPS, ThreePointSchemeType

    type ThreePointSchemeType
        type(GridType) :: grid
        real(KREAL), allocatable :: eigenVectors(:,:)
        real(KREAL), allocatable :: eigenValues(:)
        character(50) :: potential
    end type
    
contains

! initialize a new three-point scheme 
subroutine NewTPS(self, potential)
    type(ThreePointSchemeType) :: self
    character(len=*) :: potential

    ! initialize the grid
    call New(self%grid)

    ! initialize the TPS
    call AllocAndInit(self%eigenVectors, self%grid%numberOfPoints)
    call AllocAndInit(self%eigenValues, self%grid%numberOfPoints)
    self%potential = potential

end subroutine

! run the three-point scheme and write out the results to files
subroutine RunTPS(self)
    type(ThreePointSchemeType) :: self
    real(KREAL), allocatable :: matrixS(:,:), matrixV(:,:), matrixL(:,:)
    integer(KINT) :: i

    ! allocate and initialize the matrices
    call AllocAndInit(matrixS, self%grid%numberOfPoints)
    call AllocAndInit(matrixV, self%grid%numberOfPoints)
    call AllocAndInit(matrixL, self%grid%numberOfPoints)

    ! fill matrix S 
    forall (i = 1:size(matrixS, 1)) 
        matrixS(i,i) = -2
        matrixS(i, i + 1) = 1
        matrixS(i + 1, i) = 1
    endforall

    ! fill matrix V with the correct potential
    if (self%potential == "ParticleInBox") then 
        do i = 1, size(matrixV, 1)
            matrixV(i,i) = ParticleInBox(self%grid%gridPoints(i), (self%grid%interval(1))*2)
        enddo
    endif

    ! calculate matrix L and diagonalize
    matrixL = (-1/(2*(self%grid%h)**2))*matrixS + matrixV
    call diagonalize(matrixL, self%eigenVectors, self%eigenValues)

    ! checks whether the eigenvectors are negative and corrects
    if (self%eigenVectors(1,1) > self%eigenVectors(2,1)) then
        self%eigenVectors = -1*self%eigenVectors
    endif

    ! write the results to files
    call WriteToFile(self%eigenValues, "Eigenvalues")
    call WriteToFile(self%eigenVectors, "Eigenvectors")

    ! deallocate the matrices
    deallocate(matrixS, matrixV, matrixL)

end subroutine

! delete the three-point scheme
subroutine DeleteTPS(self)
    type(ThreePointSchemeType) :: self

    call Delete(self%grid)
    deallocate(self%eigenValues, self%eigenVectors)
end subroutine
  
end module ThreePointSchemeModule