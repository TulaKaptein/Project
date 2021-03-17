program main
    use Diagonalization
    use MakeGridModule
    use ThreePointSchemeModule
    implicit none
    type(Grid)           :: grid1

    ! initialize a grid based on user input
    call New(grid1)

    ! implement three-point scheme
    call ThreePointScheme(grid1)
    ! implement shooting algorithm
    ! write out output

    ! delete the grid
    call Delete(grid1)

end program
  
