!BOP
! !ROUTINE: irlkgrp 
! !INTERFACE:
subroutine  irlkgrp(v,n,r,nk,lk) 
! !INPUT/OUTPUT PARAMETERS: 
!   v  : input  k-vector input                      (in, real(3)) 
!   n  : input  number of symmetry operations       (in, integer) 
!   r  : input  rotation operators                  (in, integer(3,3,48)) 
!   nk : output number of symmetry operations       (out,integer) 
!   lk : output list of k-group symmetry operations (out,integer(48)) 
! 
! !DESCRIPTION: 
!   Determines the space group $\mathcal{G}({\bf k})$ of the allowed 
!   {\bf k}-vector, see Sect.~\ref{s:ALLOW}. 
!   Returns the list of crystallographic symmetry operations belonging   
!   to $\mathcal{G}({\bf k})$.   
! 
! !REVISION HISTORY: 
!   Created  September 2007 (Clas Persson) 
!EOP 
!BOC 
implicit none 
! external functions 
logical    irvint 
external   irvint 
! arguments 
real(8),      intent(in)  :: v(3) 
integer,      intent(in)  :: n,r(3,3,48)  
integer,      intent(out) :: nk,lk(48) 
! local variables 
integer    is,i 
real(8)    det,ri(3,3),vr(3) 
! 
nk=0 
lk(1:48)=0 
do is=1,n 
  ri(1:3,1:3)=dble(r(1:3,1:3,is)) 
  call irminv(ri,det,ri) 
  vr(1) = ri(1,1)*v(1) +ri(2,1)*v(2) +ri(3,1)*v(3)
  vr(2) = ri(1,2)*v(1) +ri(2,2)*v(2) +ri(3,2)*v(3)
  vr(3) = ri(1,3)*v(1) +ri(2,3)*v(2) +ri(3,3)*v(3)

  do i=1,3 
    vr(i)=vr(i)-v(i) 
  enddo 
  if(irvint(vr)) then 
    nk=nk+1 
    lk(nk)=is 
  endif 
end do 
! 
return 
end subroutine 
!EOC 
