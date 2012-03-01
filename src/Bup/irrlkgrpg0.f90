!BOP
! !ROUTINE: irrlkgrpg0
! !INTERFACE: 
subroutine irrlkgrpg0(v,l,n,r,t,f) 
! 
! !INPUT/OUTPUT PARAMETERS: 
!   v  : input  k-vector input                      (in, real(3)) 
!   n  : input  number of symmetry operations       (in, integer) 
!   r  : input  rotation operators                  (in, integer(3,3,48)) 
!   t  : input  translation operators               (in, real(3,48)) 
!   f  : output flag,                               (out,logical) 
! 
! !DESCRIPTION:  
!   Returns $f = true$ if the IRs of the space group $\mathcal{G}({\bf k})$
!   of the allowed {\bf k}-vector can be represented by the IRs of the
!   corresponding point group $\mathcal{G}_0({\bf k})$, otherwise $false$.   
!   See Sect.~\ref{s:ALLOW}. 
! 
! !REVISION HISTORY: 
!   Created  September 2007 (Clas Persson) 
!EOP 
!BOC 
implicit none 
! arguments 
real(8)    v(3),t(3,48) 
integer    n,r(3,3,48),l(48) 
logical    f 
! local variables 
integer    i,is,it 
real(8)    ri(3,3),vr(3),arg,dtest,ti(3),tr(3) 
real(8),   parameter :: twopi=6.2831853071795864769d0
! 
f =.true. 
do is=1,n 
  do it=1,n 
    ri(1:3,1:3)=dble(r(1:3,1:3,is)) 
    ti(1:3)=t(1:3,it) 
    tr(1)  =ri(1,1)*ti(1) +ri(1,2)*ti(2) +ri(1,3)*ti(3)
    tr(2)  =ri(2,1)*ti(1) +ri(2,2)*ti(2) +ri(2,3)*ti(3)
    tr(3)  =ri(3,1)*ti(1) +ri(3,2)*ti(2) +ri(3,3)*ti(3)
    do i=1,3 
      tr(i) = tr(i)-ti(1) 
    enddo 
    arg=-twopi*(v(1)*tr(1) +v(2)*tr(2) +v(3)*tr(3)) 
    dtest=(cos(arg)-1.0)**2 +(sin(arg))**2 
    if(dtest.gt.1d-12) f=.false. 
  enddo  
enddo  
return 
end subroutine 
!EOC  
