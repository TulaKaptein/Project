!
! Module with two different potentials: Particle in a Box and Gaussian Potential Well
!
module PotentialsModule
    use NumberKinds
    implicit none
    save
    private 
    public :: ParticleInBox, GaussianPotWell
    
contains

! function that calculates the the Particle in a Box potential for a certain x-value
real(KREAL) function ParticleInBox(x, L) result(y)
    real(KREAL), intent(in) :: x, L

    ! for each value of x within the interval, the potential is zero
    if (abs(x) <= abs(L/2)) then
        y = 0
    else
        !ERROR
        y = 10000000
    endif

end function

! function that calculates the Gaussion potential well potential for a certain x-value
real(KREAL) function GaussianPotWell(x, V0, alpha) result(y)
    real (KREAL), intent(in) :: x, V0, alpha

    y = -V0*exp(-alpha*x**2)

end function
    
end module PotentialsModule