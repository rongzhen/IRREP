!BOP
! !ROUTINE: r3vint
! !INTERFACE:
logical function r3vint(x)
! !INPUT/OUTPUT PARAMETERS:
!   x : input vector (in,real(3))
! !DESCRIPTION:
!   Returns {\it true} if the real 3-vector ${\bf x}$ is 
!   an integer 3-vector with and accuracy of $10^{-8}$, 
!   otherwise {\it false}.
!
! !REVISION HISTORY:
!   Created  September 2007 (Clas Persson)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: x(3)
! local variables
real(8) df
!
r3vint=.false.
df = abs( dble(nint(x(1))) - x(1) )&
   + abs( dble(nint(x(2))) - x(2) )&
   + abs( dble(nint(x(3))) - x(3) )
if(df.le.1.d-8) r3vint=.true.
!
return
end function
!EOC


