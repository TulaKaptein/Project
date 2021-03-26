!
! Program that implements the ShootingAlgorithm to calculate the 1D Schrodinger equation
! with a certain potential. 
!
program main
    use ShootingAlgorithmModule
    implicit none

    ! possibly put in which potential to use
    ! ParticleInBox : only need boundaries
    ! GaussianPotWell : you need v0 and alpha

    call Shooting("GaussianPotWell")

end program
  
