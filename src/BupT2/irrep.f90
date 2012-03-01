!BOP
! !ROUTINE: irrep
! !INTERFACE:
subroutine irrep(ik,evecfv,evecsv,evalfv)
! !USES:
use modmain
! !DESCRIPTION:
!   Main subroutine for determining the irreducible
!   representations $\Gamma_{\alpha}^{\mathcal{G}_0({\bf k})}$ of  
!   the electron eigenstates at each {\bf k}-point.  
!
! !REVISION HISTORY: 
!   Created  September 2007 (Clas Persson) 
!EOP 
!BOC 
implicit none 
integer,      intent(in)  :: ik 
real(8),      intent(in)  :: evalfv(nstfv,nspnfv)
complex(8),   intent(in)  :: evecfv(nmatmax,nstfv,nspnfv)
complex(8),   intent(in)  :: evecsv(nstsv,nstsv)
!
! local variables 
character(10)  cdate,ctime,czone
integer        icvalues(8)
integer       i,j,ie,ndeg,ntir(2),it 
integer       ispn,is,ia,ias,iv,ist,l,iclas 
integer       kcase,msrep(48),msneg(48),mspro(2) 
integer       nkgrp,nkgrpclass,kgrpr(3,3,48),lkgrp(48),ltest(48) 
integer       iopii(48),iopjij(48,48) 
real(8)       kgrpt(3,48),dplph,arg 
real(8)       symcrysaxc(3,48) 
!
real(8)       symaxc(3,48) 
integer       symrot(3,3,48),nsym 
real(8)       symtrn(3,48)
character(6)  symcls(48)
integer       ncls
integer       iout
!
complex(8)    kgrps(2,2,48),ztir(48,48) 
character(10) kgrpnam(2) 
character(6)  kgrpclass(48)
character(10) pgrpname(2) 
character(80) ctir(48),ttir,trlab 
logical       flg 
integer,      parameter :: maxdeg=6 
real(8),      parameter :: toldg = 1.d-6 
! 
! allocatable arrays 
character(30), allocatable :: irlab(:) 
!
!---------------------------------------------------------------! 
!                  Allocate local variable                      ! 
!---------------------------------------------------------------! 
! allocate labels of irreducible representations 
allocate(irlab(nstfv))  
iout=34  
!
!---------------------------------------------------------------! 
!    Open IRREP.OUT for symmetry ops and irreducible reprs      ! 
!---------------------------------------------------------------! 
if(ik.eq.1) then
  open(iout,file='IRREP.OUT',action='WRITE',form='FORMATTED') 
  call writetime(iout)
endif
!
!---------------------------------------------------------------! 
!   Basic info about the crystallographic symmetry operations   ! 
!---------------------------------------------------------------! 
nsymcryss=nsymcrys  
nsym=nsymcrys  
do is=1,nsym
!Clastest
  symcryss(:,:,is)=symlat(:,:,lsplsymc(is)) 
  symtran(:,is) = vtlsymc(:,is)
  symrot(:,:,is)= symlat(:,:,lsplsymc(is)) 
  symtrn(:,is)  = vtlsymc(:,is)
  lkgrp(is)     = is
enddo
!  
! properties of the symmetry operations 
call irrsymops(avec,ainv,nsym,symrot,symtrn,symaxc, &
                                 symsu2,iopii,iopjij)
!CLAStest
symcrysaxc=symaxc
!
! find crystallographic point group 
call irrpgrpname(nsym,lkgrp,symrot,symtrn,pgrpname) 
! 
! find classes of the point group symmetry operations 
call irrclasses(nsym,lkgrp,symrot,ncls,symcls) 

! write symmetry operations to IRREP.OUT 
if(ik.eq.1) & 
    call irrwritesym(iout,avec,ainv,pgrpname,nsym,symrot, &
                           symtrn,symaxc,symsu2,ncls,symcls) 

symcryss = symrot
nsymcrysclass=ncls
symcrysclass=symcls
!---------------------------------------------------------------! 
!    Calculate the eigenfunctions, and determined the           ! 
!    irreducible representations.                               ! 
!---------------------------------------------------------------! 
!do ik=1,nkpt 
! reset variables 
kcase = 0 
mspro = 0 
lkgrp = 0 
! 
! properties of the k-point 
call irrkpt(ik,kcase,msrep,msneg,mspro) 
! 
! determine the space group G(k) of the allowed k-vector 
call irrlkgrp(vkl(1:3,ik),nsymcryss,symcryss,nkgrp,lkgrp) 
write(iout,*) "ut irrlkgrp ", nkgrp,lkgrp(1)  
write(iout,*) "ut irrlkgrp ", ((symcryss(i,j,lkgrp(1)),j=1,3),i=1,3)
call flush(iout)
! 
! check the IRs of the k-point can be classified with pure  
! point group notations of Go(k) 
! call irrlkgrpg0(vkl(1:3,ik),nkgrp,lkgrp,kgrpr,kgrpt,flg) 
! 
! find point group Go(k) 
call irrpgrpname(nkgrp,lkgrp,kgrpr,kgrpt,kgrpnam) 

write(iout,*) "ut irrgrpname ", kgrpnam 
write(iout,*) "ut irrgrpname ", nkgrp 
write(iout,*) "ut irrgrpname ", ((kgrpr(i,j,lkgrp(1)),j=1,3),i=1,3)
call flush(iout)
! 
call irrclasses(nsym,lkgrp,symrot,nkgrpclass,kgrpclass)
write(iout,*) "yyyyyyyyyyyy ut irrclasses"
call flush(iout)
! find classes of the point group symmetry operations 
call irrclasses(nkgrp,lkgrp,symrot,nkgrpclass,kgrpclass) 
write(iout,*) "ut irrclasses"
call flush(iout)

!call irrclasses(nkgrp,lkgrp,kgrpr,nkgrpclass,kgrpclass) 
!write(iout,*) "ut irrclasses"
call flush(iout)
! 
if(flg) then 
! classify according to Go(k) 
! 
! get character table 
  call irrchartab(kgrpnam,ttir,ctir,ztir,ntir,dplph) 
! 
! check degeneracies due to time-reversal symmetry 
  call irrtimerev(vkl(1:3,ik),ztir,ctir,ttir,ntir,kgrpclass,nkgrp,lkgrp,&  
                         kcase,msrep,msneg,mspro,kgrpr,kgrpt,iopii,iopjij,trlab) 
! 
! calculate characters and find irreducible representation 
  call irrchar( evalfv, &
                maxdeg,vkl(1:3,ik),ttir,ctir,ztir,ntir,nkgrp,lkgrp,kgrpr,& 
                kgrpt,kgrps,ik,kgrpclass,irlab,evecfv,evecsv,toldg) 
! 
! write irreducible representations to IRREP.OUT 
write(*,*) "in irrwrite "
  call irrwrite(iout,maxdeg,ik,irlab,kgrpnam,nkgrp,nkgrpclass,kgrpclass,kgrpt,& 
                kgrpr,ttir,ctir,ztir,ntir,lkgrp,trlab) 
write(*,*) "ut irrwrite "
else 
! cannot classify the irreducible representations
  write(*,*) 
  write(*,   '("Warning(irrep): cannot classify this k-point with point group.")') 
  write(iout,'("Warning(irrep): cannot classify this k-point with point group.")') 
  write(*,*)
endif 
!---------------------------------------------------------------! 
!                Close IRREP.OUT and deallocate                 ! 
!---------------------------------------------------------------! 
deallocate(irlab) 
if(ik.eq.nkpt) then
  close(iout) 
  write(*,*) 
  write(*,'("Info(irrep):")') 
  write(*,'(" irreducible representations written to IRREP.OUT")') 
  write(*,*) 
endif
!
return 
end subroutine 
!EOC 
