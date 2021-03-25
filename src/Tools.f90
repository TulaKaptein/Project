module Tools
    use NumberKinds
    implicit none
    save
    private 
    public :: AllocAndInit, WriteToFile, Normalize

    interface AllocAndInit
        module procedure AllocAndInitArray
        module procedure AllocAndInitMatrix
    end interface

    interface WriteToFile
        module procedure WriteArray
        module procedure WriteMatrix
    end interface
    
contains

subroutine AllocAndInitArray(array, dimension)
    real(KREAL), allocatable :: array(:)
    integer(KINT) :: dimension

    allocate(array(dimension))
    array = 0

end subroutine

subroutine AllocAndInitMatrix(matrix, dimension)
    real(KREAL), allocatable :: matrix(:,:)
    integer(KINT) :: dimension

    allocate(matrix(dimension, dimension))
    matrix = 0
end subroutine

subroutine WriteArray(array, name)
    real(KREAL) :: array(:)
    character(len=*) :: name
    integer(KINT) :: i, iu

    open(newunit=iu, file=trim(name)//".txt", action="write")

    ! write output with a certain formatstring
    write(iu, *) name
    do i = 1, size(array)
        write(iu, "(f8.4)") array(i)
    enddo

    close(iu)

end subroutine

subroutine WriteMatrix(matrix, name)
    real(KREAL) :: matrix(:,:)
    character(len=*) :: name
    character(len=50) :: formatString 
    integer(KINT) :: i, iu 

    open(newunit=iu, file=trim(name)//".txt", action="write")

    write(iu,*) name

    write(formatString, "(i10)") size(matrix,2)
    formatString = "("//trim(formatString)//"f8.4)"

    do i = 1, size(matrix, 1)
        write(iu, formatString) matrix(i,:)
    enddo

    close(iu)
end subroutine

subroutine Normalize(array)
    real(KREAL) :: array(:)

    array = array/norm2(array)
end subroutine

end module Tools