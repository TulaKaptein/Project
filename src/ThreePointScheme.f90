!
! Module that uses the three-point scheme to calculate eigenvectors and eigenvalues
! for a 1D Schrödinger equation with a certain potential
!
module ThreePointSchemeModule
    use NumberKinds
    use MakeGridModule
    use Diagonalization
    use PotentialsModule
    use Tools
    implicit none
    save

    private 
    public :: ThreePointSchemeType
    public :: NewTPS, RunTPS, DeleteTPS
    public :: GetGrid, GetCertainVector, GetCertainValue, GetPotential

    type ThreePointSchemeType
        type(GridType) :: grid
        real(KREAL), allocatable :: eigenVectors(:,:)
        real(KREAL), allocatable :: eigenValues(:)
        character(50) :: potential
    end type
    
contains

!
! Initialize a new three-point scheme and delete
!
subroutine NewTPS(self, potential)
    type(ThreePointSchemeType) :: self
    character(len=*) :: potential
    integer(KINT) :: dimens

    ! initialize the grid
    call New(self%grid)

    dimens = GetDimension(self%grid)

    ! initialize the TPS
    call AllocAndInit(self%eigenVectors, dimens)
    call AllocAndInit(self%eigenValues, dimens)
    self%potential = potential

end subroutine

subroutine DeleteTPS(self)
    type(ThreePointSchemeType) :: self

    call Delete(self%grid)
    deallocate(self%eigenValues, self%eigenVectors)
end subroutine

!
! Accessors
!
function GetGrid(self)
    type(ThreePointSchemeType) :: self
    type(GridType) :: GetGrid 

    GetGrid = self%grid 
end function

function GetCertainVector(self, index)
    type(ThreePointSchemeType) :: self
    integer(KINT) :: index
    real(KREAL) :: GetCertainVector(self%grid%numberOfPoints)

    GetCertainVector = self%eigenVectors(:,index)
end function

function GetCertainValue(self, index)
    type(ThreePointSchemeType) :: self
    integer(KINT) :: index
    real(KREAL) :: GetCertainValue

    GetCertainValue = self%eigenValues(index)
end function

function GetPotential(self)
    type(ThreePointSchemeType) :: self
    character(50) :: GetPotential

    GetPotential = self%potential 
end function

!
! Run the three-point scheme and write out the results to files
!
subroutine RunTPS(self)
    type(ThreePointSchemeType) :: self
    real(KREAL), allocatable :: matrixS(:,:), matrixV(:,:), matrixL(:,:)
    integer(KINT) :: i, dimens
    real(KREAL) :: h, l, x, v0, alpha
    real(KREAL) :: potentialGraph(size(self%grid%gridPoints))

    ! get the dimension of the matrices and h
    dimens = GetDimension(self%grid)
    h = GetH(self%grid)

    ! allocate and initialize the matrices
    call AllocAndInit(matrixS, dimens)
    call AllocAndInit(matrixV, dimens)
    call AllocAndInit(matrixL, dimens)
    potentialGraph = 0

    ! fill matrix S 
    forall (i = 1:size(matrixS, 1)) 
        matrixS(i,i) = -2
    endforall
    do i = 1, size(matrixS, 1) - 1
        matrixS(i, i + 1) = 1
        matrixS(i + 1, i) = 1
    enddo

    ! fill matrix V and write the potential to a file
    if (self%potential == "ParticleInBox") then
        print *
        print *, "Particle in a box is used"
        print *
        l = abs(GetLowBound(self%grid))*2
        do i = 1, size(matrixV, 1)
            x = GetGridPoint(self%grid, i)
            matrixV(i,i) = ParticleInBox(x, l)
        enddo
    else if(self%potential == "GaussianPotWell") then
        print *
        print *, "Gaussian potential well is used"
        print *
        do i = 1, size(matrixV, 1)
            x = GetGridPoint(self%grid, i)
            v0 = 3
            alpha = 0.1
            matrixV(i,i) = GaussianPotWell(x, v0, alpha)
            potentialGraph(i) = GaussianPotWell(x,v0, alpha)
        enddo
        call WriteToFile(potentialGraph, "GaussianPotWell")
    endif

    ! calculate matrix L and diagonalize
    matrixL = -matrixS/(2*h**2) + matrixV
    call diagonalize(matrixL, self%eigenVectors, self%eigenValues)

    ! check whether the eigenvectors are negative and correct
    if (self%eigenVectors(1,1) > self%eigenVectors(2,1)) then
        self%eigenVectors = -1*self%eigenVectors
    endif

    ! write the results to files
    call WriteToFile(self%eigenValues, "Eigenvalues")
    call WriteToFile(self%eigenVectors, "Eigenvectors")

    ! deallocate the matrices
    deallocate(matrixS, matrixV, matrixL)

end subroutine
  
end module ThreePointSchemeModule