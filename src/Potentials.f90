!
! Module with two different potentials
!
module PotentialsModule
    use NumberKinds
    implicit none
    save
    private 
    public :: ParticleInBox, GaussianPotWell, HarmOsc
    
contains

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

real(KREAL) function GaussianPotWell(x, V0, alpha) result(y)
    real (KREAL), intent(in) :: x, V0, alpha

    y = -V0*exp(-alpha*x**2)

end function

real(KREAL) function HarmOsc(x, springConstant, L) result(y)
    real(KREAL), intent(in) :: x, springConstant, L

    if (abs(x) <= abs(L/2)) then
        y = 0.5*springConstant*x**2
    else
        !ERROR
        y = 1000000
    endif
end function
    
end module PotentialsModule