!BOP
! !ROUTINE: irrwrite
! !INTERFACE:
subroutine  irrwrite(iout,maxdeg,ik,irlab,kgrpnam,nkgrp,nkgrpclass,kgrpclass,kgrpt,& 
                     v,nsym, nstfv,evalfv, nspnfv, toldg, &
                     kgrpr,ttir,ctir,ztir,ntir,lkgrp,trlab) 
! !DESCRIPTION: 
!   Outputs the relevant IRs $\Gamma^{\mathcal{G}_0({\bf k})}_{\alpha}$ and  
!   $e^{-i{\bf k}({\bf \tau}_g+{\bf t}_n )}$ to {\tt IRREP.OUT}.
! 
! !REVISION HISTORY: 
!   Created  September 2007 (Clas Persson) 
!EOP 
!BOC 
implicit none 
! arguments 
real(8),      intent(in)  :: v(3),toldg
integer,      intent(in)  :: nsym,nstfv,nspnfv
real(8),      intent(in)  :: evalfv(nstfv,nspnfv)
integer,      intent(in)  :: iout
character(30),intent(in)  :: irlab(nstfv)  
integer       maxdeg,ik 
! local variables 
integer       i,j,ie,ndeg,ntir(2),it 
integer       ispn,is,ia,ias,iv,ist,l,iclas 
integer       kcase,msrep(48),msneg(48),mspro(2) 
integer       nkgrp,nkgrpclass,kgrpr(3,3,48),lkgrp(48),ltest(48) 
real(8)       kgrpt(3,48),dplph,arg 
complex(8)    kgrps(2,2,48),ztir(48,48) 
character(10) kgrpnam(2) 
character(6)  kgrpclass(48) 
character(80) ctir(48),ttir,trlab 
logical       flg 
real(8),      parameter :: twopi=6.2831853071795864769d0
! 
write(iout,'(////,2X,"k-point:",i5,5X,"(",3f7.4,")") ') ik,(v(i),i=1,3) 
write(iout,'(2X,"point group:  ",A3," ",A5)') kgrpnam(1)(1:3),kgrpnam(2)(1:5) 
write(iout,'(2X,I2," symmetry operations in ",I2," classes",/)') & 
             nkgrp, nkgrpclass 
write(iout,'(2X,"classes",6X,"symmetry operation number",9X,"exp(-i*k*tau)",$)') 
! 
ltest(1:48)=0 
do iclas=1,nkgrpclass 
  i=1 
  do while(i.le.nkgrp.and.ltest(i).ne.0)  
    i=i+1 
  enddo 
  write(iout,'(/,2X,A6,5X," ",$)') kgrpclass(i) 
  it=0 
  do is=i,nsym
    if(kgrpclass(is)(1:6).eq.kgrpclass(i)(1:6)) then 
      write(iout,'(I3,$)') lkgrp(is) 
      ltest(is)=iclas 
      it=it+1 
    endif 
  enddo 
  arg=-twopi*( v(1)*kgrpt(1,i) & 
              +v(2)*kgrpt(2,i) & 
              +v(3)*kgrpt(3,i) ) 
  do j=1,11-it 
    write(iout,'("   ",$)')  
  enddo   
  write(iout,'(X,"(",2F6.3,")",$)') cos(arg),sin(arg) 
enddo 
! 
write(iout,'(///,5X,A70)') ttir(1:70)  
do is=1,ntir(1) 
  if(is.eq.ntir(2)+1) write(iout,'(4x,81a)') ('-',it=1,8+ntir(2)*6) 
  write(iout,'(4x,a4,x,a10,$)') ctir(is)(1:4),ctir(is)(5:14) 
  write(iout,'(11a6)') (ctir(is)(8+it*6+1:8+it*6+7),it=1,ntir(2)-1) 
enddo 
if(abs(dplph-(twopi/8.d0)).le.1e-6) write(iout,'(4x,"d=exp(2pi*i/8)")')  
if(abs(dplph-(twopi/3.d0)).le.1e-6) write(iout,'(4x,"e=exp(2pi*i/3)")')  
 
if(trlab(1:4).ne."    ") write(iout,'(                     & 
    /,4X,"Degeneration due to time-reversal symmetry: ",   & 
    /,4X,A60)') trlab(1:60) 
! 
write(iout,'(//,2X,"band degen   eigval     irred repr")') 
ie=1 
do while(ie.le.nstfv) 
  ndeg=1 
  do while((ie+ndeg).le.nstfv.and. & 
           (evalfv(ie+ndeg,1)-evalfv(ie,1)).lt.toldg) 
    ndeg=ndeg+1 
  enddo 
  if((ie+ndeg-1).lt.nstfv) then
    write(iout,'(2I5,3X,F10.6,7X,A30)') ie,ndeg,evalfv(ie,1),irlab(ie)(1:30) 
  endif
  ie=ie+ndeg 
enddo  
write(iout,'(/)') 
! 
return 
end subroutine 
!EOC 
