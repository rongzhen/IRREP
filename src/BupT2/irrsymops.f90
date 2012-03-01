!BOP
! !ROUTINE: irrsymops
! !INTERFACE:
subroutine  irrsymops(avec,ainv,nsym,symrot,symtrn,symaxc, &
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
real(8)       r3mdet 
logical       r3vint 
external      r3mdet,r3vint 
! arguments 
real(8),      intent(in)  :: ainv(3,3), avec(3,3)
integer,      intent(in)  :: symrot(3,3,48),nsym 
real(8),      intent(in)  :: symtrn(3,48)  
real(8),      intent(out) :: symaxc(3,48) 
complex(8),   intent(out) :: symsu2(2,2,48) 
integer,      intent(out) :: iopii(48),iopjij(48,48) 
! local variables 
integer       n,is,it 
real(8)       tau(3),tav(2) 
real(8)       ds,det,ang(3),rot(3,3),ros(3,3),roi(3,3) 
! 
do is=1,nsym 
!---------------------------------------------------------------! 
!          Determine proper or improper rotation                ! 
!---------------------------------------------------------------! 
  rot(1:3,1:3)=dble(symrot(1:3,1:3,is)) 
  call  r3mm(rot,ainv,rot) 
  call  r3mm(avec,rot,rot) 
  det = r3mdet(rot) 
  if(abs(abs(det)-1.0).ge.1d-6) then 
    write(*,*) 
    write(*,'("Error(irrsymops): det(R)=",F8.5)') det 
    write(*,*) 
    stop   
  endif 
  if(det.lt.-0.5) rot(1:3,1:3)=-rot(1:3,1:3)  
  call irrrotaxis(rot,symaxc(1,is)) 
!---------------------------------------------------------------! 
!           SU(2) matrices of the proper rotations              ! 
!---------------------------------------------------------------! 
  call euler(rot(1,1),ang(1)) 
  call su2(ang(1:3),symsu2(1,1,is)) 
!---------------------------------------------------------------! 
!           Find iopii(i=is)=n for which Rn=Pi*Pi               ! 
!---------------------------------------------------------------! 
  rot(1:3,1:3)=dble(symrot(1:3,1:3,is)) 
  tau(1:3)    =symtrn(1:3,is) 
  call r3mm(rot,rot,ros) 
  call findmmi(ros,dble(symrot),nsym,n,ds) 
  call r3mv(rot,tau,tav) 
  tav(1) = (tav(1)-tau(1)) -symtrn(1,n) 
  tav(2) = (tav(2)-tau(2)) -symtrn(2,n) 
  tav(3) = (tav(3)-tau(3)) -symtrn(3,n) 
  if(ds.lt.1e-6.and.r3vint(tav)) then 
    iopii(is)=n 
  else 
    write(*,*) 
    write(*,'("Error(irrsymops): cannot find Pn=Pj*Pj")') 
    write(*,*) 
    stop 
  end if 
!---------------------------------------------------------------! 
!     Find iopjij(i=is,j=it)=n for which Rn=Rj*Ri*inv(Rj)       ! 
!---------------------------------------------------------------! 
  rot(1:3,1:3)=dble(symrot(1:3,1:3,is)) 
  do it=1,nsym 
    ros(1:3,1:3)=dble(symrot(1:3,1:3,it)) 
    call r3minv(ros,roi) 
    call r3mm(rot,roi,roi) 
    call r3mm(ros,roi,roi) 
    call findmmi(roi,dble(symrot),nsym,n,ds) 
    if(ds.lt.1d-6) then 
      iopjij(is,it)=n 
    else 
      write(*,*) 
      write(*,'("Error(irrsymops): no Pn=Pj*Pi*inv(Rj)")') 
      write(*,*) 
      stop 
    end if 
  end do 
!
end do 
return 
end subroutine 
!EOC 
