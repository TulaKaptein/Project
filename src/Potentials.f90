module PotentialsModule
    use NumberKinds
    implicit none
    save
    private 
    public :: ParticleInBox
    
contains

real function ParticleInBox(x, L) result(y)
    real(KREAL), intent(in) :: x
    real(KREAL), intent(in) :: L

    ! for each value of x within the interval, the potential is zero
    if (abs(x) <= abs(L/2)) then
        y = 0
    else
        !ERROR
        y = 10000000
    endif

end function
    
end module PotentialsModule