!BOP
! !ROUTINE: su2
! !INTERFACE:
subroutine  su2(x,s)
! !INPUT/OUTPUT PARAMETERS:
!   x : input  angles        (in, real(3))
!   s : output SU(2) matrix  (out,complex(2,2))
!
! !DESCRIPTION:
!   Calculates SU(2) matrix from the Euler's angles. Matrix is 
!   defined to be positive: $\rm{Re}|\rm{s}_{11}| \ge 0$.
!   See routine {\tt euler} for definition of rotation angles.
!
! !REVISION HISTORY:
!   Created  September 2007 (Clas Persson)
!EOP
!BOC
implicit none
! arguments
real(8),    intent(in)  :: x(3)
complex(8), intent(out) :: s(2,2)
! local variables
complex(8), parameter   :: zi=(0.d0,1.d0)
complex(8)  a(4)
real(8)     b(2)
!
s(1,1)= dcos(x(2)*0.5)*exp( zi*(x(3)+x(1))*0.5)
s(1,2)= dsin(x(2)*0.5)*exp( zi*(x(3)-x(1))*0.5)
s(2,1)=-dsin(x(2)*0.5)*exp(-zi*(x(3)-x(1))*0.5)
s(2,2)= dcos(x(2)*0.5)*exp(-zi*(x(3)+x(1))*0.5)
!
! use positive matrices
if(real(s(1,1)).lt.0.d0) s(1:2,1:2) = -s(1:2,1:2)
!
return
end subroutine
!EOC

