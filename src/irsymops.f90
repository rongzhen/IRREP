!BOP
! !ROUTINE: irsymops
! !INTERFACE:
subroutine  irsymops(fl,avec,ainv,nsym,symrot,symtrn,symaxc, &
                                       symsu2,iopii,iopjij)
!
! !INPUT/OUTPUT PARAMETERS: 
!   avec   : lattice vectors                          (in ,real(3,3))
!   ainv   : inverse of lattice vector                (in ,real(3,3))
!   nsym   : number of space-group operations         (in, integer) 
!   symrot : crystallographic symmetry rotations R    (in ,integer(3,3,48)) 
!   symtrn : crystallographic symmetry translations t (in ,real(3,48)) 
!   symaxc : rotations properties, where              (out,real(3,48)) 
!                (1,i) = theta,  Cartesian coordinates 
!                (2,i) = phi,    Cartesian coordinates 
!                (3,i) = 1 for clockwise,  
!                       -1 for counterclockwise, and  
!                        0 for unit and 180 degree rotations 
!   symsu2 : SU(2) matrices                           (out,complex(2,2,48))
!   iopii  : n, for which Pn=Pi*inv(Pi)               (out,integer(48)) 
!   iopjij : n, for which Pn=Pj*Pi*inv(Pj)            (out,integer(48,48)) 
!                
! !DESCRIPTION: 
!   Determines the properties of the crystallographic symmetry 
!   operations ${\bf P}_g$ = $\{{\bf R}_g|{\bf t}_g\}$ and 
!   generates ${\bf U}({\bf R}_g^p)$. 
!    
! !REVISION HISTORY: 
!   Created  September 2007 (Clas Persson) 
!EOP 
!BOC 
implicit none 
! external functions 
logical        irvint 
external       irvint 
! arguments 
logical,       intent(inout) :: fl(10)
real(8),       intent(in)    :: ainv(3,3), avec(3,3)
integer,       intent(in)    :: symrot(3,3,48),nsym 
real(8),       intent(in)    :: symtrn(3,48)  
real(8),       intent(out)   :: symaxc(3,48) 
complex(8),    intent(out)   :: symsu2(2,2,48) 
integer,       intent(out)   :: iopii(48),iopjij(48,48) 
! local variables 
integer        n,is,it 
real(8)        tau(3),tav(3) 
real(8)        ds,det,ang(3),rot(3,3),ros(3,3),roi(3,3) 
! 
fl(8)=.false. 
fl(9)=.false. 
do is=1,nsym 
!---------------------------------------------------------------! 
!          Determine crystal symmetry                           ! 
!---------------------------------------------------------------! 
  rot(1:3,1:3)=dble(symrot(1:3,1:3,is)) 
  tau(1:3)    =symtrn(1:3,is) 
  if(dabs(tau(1))+dabs(tau(2))+dabs(tau(3)).gt.5e-6) fl(8)=.true.  
  if(abs(rot(1,1)+1)+abs(rot(1,2))  +abs(rot(1,3)  ) &
    +abs(rot(2,1))  +abs(rot(2,2)+1)+abs(rot(2,3)  ) &
    +abs(rot(3,1))  +abs(rot(3,2))  +abs(rot(3,3)+1) &
                                           .gt.5e-6) fl(9)=.true.
!---------------------------------------------------------------! 
!          Determine proper or improper rotation                ! 
!---------------------------------------------------------------! 
  rot(1:3,1:3)=dble(symrot(1:3,1:3,is)) 
  call  irmm(rot,ainv,rot) 
  call  irmm(avec,rot,rot) 
  call  irminv(rot,det,ros)
  if(abs(abs(det)-1.0).ge.1d-6) then 
    write(*,*) 
    write(*,'("Error(irsymops): |Rg|=",F8.5)') det 
    write(*,*) 
    stop   
  endif 
  if(det.lt.-0.5) rot(1:3,1:3)=-rot(1:3,1:3)  
  call irrotaxis(rot,symaxc(1,is)) 
!---------------------------------------------------------------! 
!           SU(2) matrices of the proper rotations              ! 
!---------------------------------------------------------------! 
  call irsu2(rot,ang,symsu2(1,1,is)) 
!---------------------------------------------------------------! 
!           Find iopii(i=is)=n for which Rn=Pi*Pi               ! 
!---------------------------------------------------------------! 
  rot(1:3,1:3)=dble(symrot(1:3,1:3,is)) 
  tau(1:3)    =symtrn(1:3,is) 
!write(34,*) " ttt = ",is
!write(34,'(4F8.4)')  (rot(1,it),it=1,3),tau(1)
!write(34,'(4F8.4)')  (rot(2,it),it=1,3),tau(2)
!write(34,'(4F8.4)')  (rot(3,it),it=1,3),tau(3)

  call irmm(rot,rot,ros) 
  call irfindmmi(ros,dble(symrot),nsym,n,ds) 
  tav(1) = rot(1,1)*tau(1) +rot(1,2)*tau(2) +rot(1,3)*tau(3)
  tav(2) = rot(2,1)*tau(1) +rot(2,2)*tau(2) +rot(2,3)*tau(3)
  tav(3) = rot(3,1)*tau(1) +rot(3,2)*tau(2) +rot(3,3)*tau(3)
!
  tav(1) = (tav(1)-tau(1)) -symtrn(1,n) 
  tav(2) = (tav(2)-tau(2)) -symtrn(2,n) 
  tav(3) = (tav(3)-tau(3)) -symtrn(3,n) 
  if(ds.lt.1e-6.and.irvint(tav)) then 
    iopii(is)=n 
  else 
    write(*,*) 
    write(*,'("Error(irsymops): cannot find Pn=Pj*Pj")') 
    write(*,*) 
    stop 
  end if 
!---------------------------------------------------------------! 
!     Find iopjij(i=is,j=it)=n for which Rn=Rj*Ri*inv(Rj)       ! 
!---------------------------------------------------------------! 
  rot(1:3,1:3)=dble(symrot(1:3,1:3,is)) 
  do it=1,nsym 
    ros(1:3,1:3)=dble(symrot(1:3,1:3,it)) 
    call irminv(ros,det,roi) 
    call irmm(rot,roi,roi) 
    call irmm(ros,roi,roi) 
    call irfindmmi(roi,dble(symrot),nsym,n,ds) 
    if(ds.lt.1d-6) then 
      iopjij(is,it)=n 
    else 
      write(*,*) 
      write(*,'("Error(irsymops): no Pn=Pj*Pi*inv(Rj)")') 
      write(*,*) 
      stop 
    end if 
  end do 
!
end do 
return 
end subroutine 
!EOC 
