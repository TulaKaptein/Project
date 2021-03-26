program main
    use ShootingAlgorithmModule
    implicit none

    ! possibly put in which potential to use
    ! ParticleInBox : only need boundaries
    ! GaussianPotWell : you need v0 and alpha
    ! HarmOsc : you need springConstant and boundaries

    !!!! HARMOSC NOT WORKING YET
    call Shooting("GaussianPotWell")

end program
  
