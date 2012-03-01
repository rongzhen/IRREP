!BOP
! !ROUTINE: findmmi
! !INTERFACE:
subroutine  findmmi(a,b,n,i,d)
! !INPUT/OUTPUT PARAMETERS:
!   a : input  matrix,            (in,  real(3,3))
!   b : input  matrix array,      (in,  real(3,3,n))
!   n : input  size of b          (in,  integer)
!   i : output matrix number of b (out, integer)
!   d : output sum of |y(i)-x|    (out, real)
!
! !DESCRIPTION:
!   Finds $3\times 3$ real matrix ${\bf b}(i)$  $i=1,\ldots,n$ 
!   that minimizes
!   $d=\sum_{jj'}|{\bf b}_{jj'}(i) - {\bf a}_{jj'}|$. 
!   Returns $i$ and $d$, where $d=0$ if ${\bf b}(i)\equiv{\bf a}$.
!     
! !REVISION HISTORY:
!   Created  September 2007 (Clas Persson)
!EOP
!BOC
implicit none 
! arguments
integer,    intent(in)  :: n
integer,    intent(out) :: i
real(8),    intent(out) :: d
real(8),    intent(in)  :: a(3,3),b(3,3,n)
!local variables
real(8)     ds
integer     is
!
d=1.d+34 
do is=1,n
  ds = abs(b(1,1,is)-a(1,1))&
     + abs(b(1,2,is)-a(1,2))&
     + abs(b(1,3,is)-a(1,3))&
     + abs(b(2,1,is)-a(2,1))&
     + abs(b(2,2,is)-a(2,2))&
     + abs(b(2,3,is)-a(2,3))&
     + abs(b(3,1,is)-a(3,1))&
     + abs(b(3,2,is)-a(3,2))&
     + abs(b(3,3,is)-a(3,3))
  if(ds.lt.d) then 
    d=ds
    i=is
  endif
enddo 
return
end subroutine
!EOC

