!BOP 
! !ROUTINE: irrchar 
! !INTERFACE: 
subroutine irrchar(evalfv, & 
                   maxdeg,v,ttir,ctir,ztir,ntir,nkgrp,lkgrp, &
                   kgrpr,kgrpt,kgrps,ik,& 
                   kgrpclass,irlab,evecfv,evecsv,toldg) 
! !USES: 
use modmain  
! !DESCRIPTION:  
!   Calculates the characters $\chi^{\mathcal{G}}_{\alpha}({\bf P}_g)$ 
!   of the representation matrices
!   for each energy states. Degeneracy of eigenstates is determined
!   by the tolerance variable {\tt toldg} which can be changed in     
!   the routine {\tt irrep}. 
! 
! !REVISION HISTORY: 
!   Created  September 2007 (Clas Persson) 
!EOP 
!BOC 
implicit none 
character(30),intent(out) :: irlab(nstfv)  
real(8),      intent(in)  :: evalfv(nstfv,nspnfv)
complex(8),   intent(in)  :: evecfv(nmatmax,nstfv,nspnfv)
complex(8),   intent(in)  :: evecsv(nstsv,nstsv)
! external functions 
real(8)   r3mdet 
logical   r3vint 
external  r3mdet,r3vint 
! arguments 
character(80) ctir(48),ttir 
complex(8)    ztir(48,48),xi,xr(48),zir 
integer       maxdeg,ntir(2),lkgrp(48),nkgrp,kgrpr(3,3,48),ik 
real(8)       kgrpt(3,48),v(3),toldg  
complex(8)    kgrps(2,2,48),csum 
character(6)  kgrpclass(48) 
! 
integer       n,nsp,r(3,3,48) 
logical       f,flg1 
! initialise universal variables 
! local variables 
integer       ig,in,i,ir,is,it,i1,i2,iu,j,j1,j2,idi,ip 
integer       nvec,isp,ie,ndeg 
integer       iwrt,nj1i(100),i3 
real(8)       ri(3,3),vg(3),vr(3),arg,dtest,ti(3),tr(3) 
real(8)       dph,sum,ds 
! 
! allocatable arrays 
complex(8),   allocatable :: ph(:,:) 
integer,      allocatable :: lv(:,:) 
complex(8),   allocatable :: pvecv(:,:,:) 
complex(8),   allocatable :: rvecv(:,:,:) 
complex(8),   allocatable :: xm(:,:) 
! allocate local allocatable variables 
nsp = 1
isp = 1
write(*,*) "WARNING SPIN = ",nsp
nvec=ngk(ik,isp) 
write(*,*) "WARNING SPIN = ",nsp

allocate(ph(nkgrp,nvec)) 
allocate(lv(nkgrp,nvec)) 
allocate(pvecv(nvec,nstfv,2)) 
allocate(rvecv(nvec,nstfv,2)) 
allocate(xm(48,nstfv)) 
xm(1:48,1:nstfv)=cmplx(9999.0, 9999.0) 
! 
!----------------------------------------------------------! 
!           find G' such that (k+G)inv(R)=k+G'             ! 
!----------------------------------------------------------! 
do is=1,nkgrp 
  ri(1:3,1:3)=dble(symcryss(1:3,1:3,lkgrp(is))) 
  call r3minv(ri,ri) 
  do ip=1,nvec 
    vg(1:3) = vgkl(1:3,ip,ik,isp)  
    call r3mtv(ri,vg,vr) 
    in=1 
    do while(( abs(vr(1)-vgkl(1,in,ik,isp)) & 
              +abs(vr(2)-vgkl(2,in,ik,isp)) & 
              +abs(vr(3)-vgkl(3,in,ik,isp)) ).gt.1d-8) 
      in=in+1 
      if(in.eq.nvec+1)then 
        write(*,*) 
        write(*,'("Error(irrchar): no (k+G)inv(R)")') 
        write(*,*) 
        stop 
      endif 
    enddo 
    arg=-twopi*( (vgkl(1,in,ik,isp)-v(1))*symtran(1,lkgrp(is)) & 
                +(vgkl(2,in,ik,isp)-v(2))*symtran(2,lkgrp(is)) & 
                +(vgkl(3,in,ik,isp)-v(3))*symtran(3,lkgrp(is)) )                 
    ph(is,ip)=cmplx(cos(arg),sin(arg)) 
    lv(is,ip)=in 
  enddo 
enddo 
!----------------------------------------------------------! 
!              normalize the plane waves                   ! 
!----------------------------------------------------------! 
pvecv(1:nvec,1:nstfv,1:2)=cmplx(0.d0,0.d0)
rvecv(1:nvec,1:nstfv,1:2)=cmplx(0.d0,0.d0) 


pvecv(:,:,1) =evecfv(:,:,1) 

!clas
pvecv(ip,j1,1) =evecfv(ip,j1,1) 
! 
do j=1,nstfv
  sum=0.d0 
  do ip=1,nvec 
    do isp=1,nsp  
      sum=sum +abs(pvecv(ip,j,isp))**2 
    enddo 
  enddo 
  do ip=1,nvec 
    do isp=1,nsp   
     pvecv(ip,j,isp)=(pvecv(ip,j,isp))/sqrt(sum) 
   enddo 
 enddo 
enddo 
!----------------------------------------------------------! 
!     characters for all energies and symm. ops.           ! 
!----------------------------------------------------------! 
write(*,*) "nst nsp nstfv ",nstfv, nsp, nstfv 
write(*,*) "char ",(evalfv(ie,1),ie=1,10)

!Clas0
isp = 1
nsp = 1
write(*,*) " WARNING jjjj"
!Clas1



ie=1 
do while(ie.le.nstfv) 

!Clas0
isp = 1
nsp = 1
!Clas1

  ndeg=1 
  do while((ie+ndeg).le.(nstfv-6).and. & 
           (evalfv(ie+ndeg,isp)-evalfv(ie,isp)).lt.toldg) 
    ndeg=ndeg+1 
  enddo 
  if(ndeg.gt.maxdeg) then 
!    write(*,*)  
!    write(*,*) "y1 ", ik,isp, ie, ndeg 
!    write(*,*) "y2 ", evalfv(ie+ndeg,isp), evalfv(ie,isp)
!    write(*,'("Error(irrchar): ndeg > maxdeg")') 
!    write(*,'("                increase maxdeg in irrep.f90")') 
!    write(*,*)  
!    stop 
     ie=nstfv+1
  else 
! 
! for all symmetry operations of little k-group 

!clas
    do is=1,nkgrp  
      j2=0 
      do j=ie,ie+ndeg-1 
        j2=j2+1 
        do ip=1,nvec 
          if(nspinor.eq.1) then 
            rvecv(ip,j2,1:2)=pvecv(ip,j,1:2)*ph(is,ip)   
          else 
write(*,*) "WARNING not SO yet"
            rvecv(ip,j2,1)=pvecv(ip,j,1)*symsu2(1,1,lkgrp(is))& 
                          +pvecv(ip,j,2)*symsu2(1,2,lkgrp(is)) 
            rvecv(ip,j2,2)=pvecv(ip,j,1)*symsu2(2,1,lkgrp(is))& 
                          +pvecv(ip,j,2)*symsu2(2,2,lkgrp(is)) 
            rvecv(ip,j2,1:2)=rvecv(ip,j2,1:2)*ph(is,ip)  
          endif 
        enddo 
      enddo 
! 
! character 
      xm(is,ie)=cmplx(0.d0,0.d0) 
      do i1=1,ndeg 
        do i2=1,ndeg  
          csum=cmplx(0.d0,0.d0)
          do i=1,nvec 
            do isp=1,nsp  
              csum=csum+conjg(pvecv(lv(is,i),ie+i2-1,isp)) & 
                             *rvecv(i,i1,isp)  
            enddo 
          enddo 
          if(i1.eq.i2) xm(is,ie)=xm(is,ie)+csum 
        enddo 
      enddo 
      xm(is,ie:ie+ndeg-1)=xm(is,ie) 
    enddo 
!
!----------------------------------------------------------! 
!     identify the irreducible representations,            ! 
!     for iwrt:th element in the ir:th class               ! 
!----------------------------------------------------------! 
    xr(1:48)=cmplx(0.d0,0.d0) 
    do ir=1,ntir(2) 
      iwrt=0 
      do is=1,nkgrp 
        if(kgrpclass(is)(1:6).eq.ttir(6*(ir-1)+10:6*ir+9)) then 
          iwrt=iwrt+1     
          xi=xm(is,ie) 
          if(iwrt.eq.1) then 
            xr(ir)=xi 
          else 
            if(abs(xi-xr(ir)).gt.1e-6) then 
!              write(*,*)  
              write(*,'("Error(irrchar): X(Ri) not X(Rj)")') 
!              write(*,*)  
!              stop 
            endif 
          endif 
        endif 
      enddo 
    enddo 
! 
    flg1=.true. 
    i1=0 
    do while(i1.lt.maxdeg.and.flg1) 
      i1=i1+1 
      nj1i(1:maxdeg)=0 
      i2=-1 
      do while(i2.lt.(ntir(1)**i1-1).and.flg1) 
        i2=i2+1 
        nj1i(1)=mod(i2,ntir(1))+1 
        do i3=2,i1 
          if(mod(i2,ntir(1)**(i3-1)).eq.0) then 
            nj1i(i3)=nj1i(i3)+1 
            if(nj1i(i3).gt.ntir(1)) nj1i(i3)=1 
          endif 
        enddo 
        ds=0.d0 
        do ir=1,ntir(2) 
          zir=cmplx(0.d0,0.d0) 
          do i3=1,i1 
            zir=zir+ ztir(nj1i(i3),ir) 
          enddo 
          ds=ds+real(abs(xr(ir)-zir)) 
        enddo 
        if(ds.lt.dble(ntir(1))*1.0d-3) flg1=.false. 
      enddo 
    enddo 
!
!----------------------------------------------------------! 
!                  label for output                        ! 
!----------------------------------------------------------! 
    irlab(ie)(:)=' ' 
    if(flg1.or.(ie+ndeg).gt.nstfv) then 
      write(irlab(ie)(1:2),'("--")')  
    else 
      in=1 
      do i3=i1,1,-1 
        if(i3.ne.i1) then 
          write(irlab(ie)(in:in+1),'("+",$)')  
          in=in+1 
        endif 
        write(irlab(ie)(in:in+3),'(A3,$)') CTIR(nj1i(i3))(1:3) 
        in=in+3 
      enddo 
    endif 
    ie=ie+ndeg 
  endif 
enddo 
! 
deallocate(ph,lv,pvecv,rvecv,xm) 
return 
end subroutine 
!EOC 
