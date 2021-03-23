module integration_module

   implicit none 
   save
   private

   public :: Newton_cotes 

   contains

   subroutine Newton_cotes(f,h,a,b,v_int)

   ! evaluates an integral over a 1-dimensional function f, tabulated on n equidistant grid 
   ! using composite Newton-Cotes formulas. The subroutine checks the number of grid-points and 
   ! decides which routine is used, giving priority to higher-order Newton-Cotes formulas, starting
   ! from a 7-point rule.

   real(8), intent(in)  :: f(:)              ! integrand
   real(8), intent(in)  :: h                 ! mesh-spacing 
   integer, intent(in)  :: a,b               ! start-/end-point of integration 
   real(8), intent(out) :: v_int             ! definite integral 

   real(8) :: prefactor
   integer :: i,N

   N = b - a + 1 

   if (mod(N-1,6).eq.0) then 

!      write(*,*) '7-point'

      prefactor = h/140.0d0
      do i = a+3,b-3,6
         v_int = v_int +  41.0d0 * f(i-3) +  41.0d0 * f(i+3) &
                       + 216.0d0 * f(i-2) + 216.0d0 * f(i+2) &
                       +  27.0d0 * f(i-1) +  27.0d0 * f(i+1) &
                       + 272.0d0 * f(i)
      enddo

   else if (mod(N-1,5).eq.0) then 

!      write(*,*) '6-point'

      prefactor = 5.0d0 * h/288.0d0
      do i = a+3,b-2,5
         v_int = v_int + 19.0d0 * f(i-3) + 19.0d0 * f(i+2) &
                       + 75.0d0 * f(i-2) + 75.0d0 * f(i+1) &
                       + 50.0d0 * f(i-1) + 50.0d0 * f(i) 
      enddo

   else if (mod(N-1,4).eq.0) then 

!      write(*,*) '5-point'

      prefactor = 2.0d0 * h/45.0d0
      do i = a+2,b-2,4
         v_int = v_int +  7.0d0 * f(i-2) +  7.0d0 * f(i+2) &
                       + 32.0d0 * f(i-1) + 32.0d0 * f(i+1) &
                                         + 12.0d0 * f(i) 
      enddo

   else if (mod(N-1,3).eq.0) then

!      write(*,*) '4-point'

      prefactor = 3.0d0 * h/8.0d0
      do i = a+1,b-2,3
         v_int = v_int + f(i-1) +  3.0d0 * f(i) + 3.0d0 * f(i+1) + f(i+2)
      enddo

   else if (mod(N-1,2).eq.0) then

!      write(*,*) '3-point'

      prefactor = 1.0d0 * h/3.0d0
      do i = a+1,b-1,2
         v_int = v_int + f(i-1) +  4.0d0 * f(i) + f(i+1)
      enddo

   else

!      write(*,*) '2-point'

      prefactor = h/2.0d0
      do i = a+1,b-1
         v_int = v_int + 2.0d0 * f(i)
      enddo 
     v_int = v_int + f(a) + f(b)

   endif

   v_int = v_int * prefactor

   end subroutine

end module

!program integration
!   use integration_module
!   implicit none 

!   real(8) :: res2,h
!   integer, parameter :: Npoints =103, a = 1,b = 103
!   integer :: i
!   real(8), save,allocatable :: function(:) 
!   real(8) :: length

!   allocate(function(Npoints))

!   length = 10.0d0 

!  h = length/(real(Npoints,8)-1) 

!   function = 0.0d0

!   do i = 1,Npoints
!      function(i) = (i-1)*h + 5
!   enddo

!   call Newton_cotes(function,h,a,b,res2)
!  write(*,*) 'res (six-point) = ', res2

!   deallocate(function)
!end program

