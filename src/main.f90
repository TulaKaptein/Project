!
! Program that implements the ShootingAlgorithm to calculate the 1D Schrodinger equation
! with a certain potential. 
!
! ParticleInBox     : the boundaries are equal to the interval
! GaussianPotWell   : v0 = 3.0, alpha = 0.1
!
program main
    use ShootingAlgorithmModule
    implicit none

    ! start the shooting algorithm with a certain potential
    call Shooting("ParticleInBox")

end program
  
