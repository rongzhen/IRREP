!BOP
! !ROUTINE: irmm
! !INTERFACE:
subroutine  irmm(a,b,c)
! !INPUT/OUTPUT PARAMETERS:
!   a : input  matrix                 (in,real(3,3))
!   b : input  matrix                 (in,real(3,3))
!   c : output matrix                 (out,real(3,3))
! !DESCRIPTION:
!   Multiplies $3\times 3$ matrices: ${\bf ab}$ = ${\bf c}$.
!
! !REVISION HISTORY:
!   Created  September 2007 (Clas Persson)
!EOP
!BOC
implicit none
! arguments
real(8),      intent(in)   :: a(3,3)
real(8),      intent(in)   :: b(3,3)
real(8),      intent(out)  :: c(3,3)
! local variables
real(8)       tm(3,3)
!
tm(1,1) = a(1,1)*b(1,1)+a(1,2)*b(2,1)+a(1,3)*b(3,1)
tm(1,2) = a(1,1)*b(1,2)+a(1,2)*b(2,2)+a(1,3)*b(3,2)
tm(1,3) = a(1,1)*b(1,3)+a(1,2)*b(2,3)+a(1,3)*b(3,3)
tm(2,1) = a(2,1)*b(1,1)+a(2,2)*b(2,1)+a(2,3)*b(3,1)
tm(2,2) = a(2,1)*b(1,2)+a(2,2)*b(2,2)+a(2,3)*b(3,2)
tm(2,3) = a(2,1)*b(1,3)+a(2,2)*b(2,3)+a(2,3)*b(3,3)
tm(3,1) = a(3,1)*b(1,1)+a(3,2)*b(2,1)+a(3,3)*b(3,1)
tm(3,2) = a(3,1)*b(1,2)+a(3,2)*b(2,2)+a(3,3)*b(3,2)
tm(3,3) = a(3,1)*b(1,3)+a(3,2)*b(2,3)+a(3,3)*b(3,3)
c = tm  
!
return
end subroutine
!EOC
