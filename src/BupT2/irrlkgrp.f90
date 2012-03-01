!BOP
! !ROUTINE: irrlkgrp 
! !INTERFACE:
subroutine irrlkgrp(v,n,r,nk,lk) 
! !INPUT/OUTPUT PARAMETERS: 
!   v  : input  k-vector input                      (in, real(3)) 
!   n  : input  number of symmetry operations       (in, integer) 
!   r  : input  rotation operators                  (in, integer(3,3,48)) 
!   nk : output number of symmetry operations       (out,integer) 
!   lk : output list of k-group symmetry operations (out,integer(48)) 
! 
! !DESCRIPTION: 
!   Determines the space group $\mathcal{G}({\bf k})$ of the allowed 
!   {\bf k}-vector, i.e., those ${\bf R}_g$ having the property 
!   ${\bf k}$$\cdot$${\bf R}_g^{-1}$ = ${\bf k}$+${\bf K}$. 
!   Returns the list of crystallographic symmetry operations belonging   
!   to $\mathcal{G}({\bf k})$.   
! 
! !REVISION HISTORY: 
!   Created  September 2007 (Clas Persson) 
!EOP 
!BOC 
implicit none 
! external functions 
logical    r3vint 
external   r3vint 
! arguments 
real(8)    v(3) 
integer    r(3,3,48),nk,lk(48) 
! local variables 
integer    is,i,n 
real(8)    ri(3,3),vr(3) 
! 
nk=0 
lk(1:48)=0 
do is=1,n 
  ri(1:3,1:3)=dble(r(1:3,1:3,is)) 
  call r3minv(ri,ri) 
  call r3mtv(ri,v,vr) 
  do i=1,3 
    vr(i)=vr(i)-v(i) 
  enddo 
  if(r3vint(vr)) then 
    nk=nk+1 
    lk(nk)=is 
  endif 
end do 
! 
return 
end subroutine 
!EOC 
