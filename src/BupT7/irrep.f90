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
real(8),      allocatable ::  ener(:,:)
complex(8),   allocatable ::  evec(:,:,:)
complex(8),   allocatable ::  evey(:,:)

!
real(8),      intent(in)  :: evalfv(nstfv,nspnfv)
complex(8),   intent(in)  :: evecfv(nmatmax,nstfv,nspnfv)
complex(8),   intent(in)  :: evecsv(nstsv,nstsv)
!
! local variables 
integer       i,j,ie,ndeg,ntir(2),it 
integer       ispn,is,ia,ias,iv,ist,l,iclas 
integer       kcase,msrep(48),msneg(48),mspro(2) 
integer       nkgrp,nkgrpclass,kgrpr(3,3,48),lkgrp(48),ltest(48) 
integer       iopii(48),iopjij(48,48) 
real(8)       kgrpt(3,48),dplph,arg 
!
real(8)       symaxc(3,48) 
integer       symrot(3,3,48),nsym 
real(8)       symtrn(3,48),                test(3,3)
complex(8)    symsu2(2,2,48)
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
logical       fl(10)
character(80) title
integer,      parameter :: maxdeg= 6 
real(8),      parameter :: toldg = 9.d-6 
! 
! allocatable arrays 
character(30), allocatable :: irlab(:) 
!
!---------------------------------------------------------------! 
!                  Allocate local variable                      ! 
!---------------------------------------------------------------! 
! allocate labels of irreducible representations 
allocate(irlab(nstfv))  
allocate(ener(nstfv,nspnfv) )
allocate(evec(nmatmax,nstfv,nspnfv) ) 
allocate(evey(nstsv,nstsv) ) 

ener(:,:)   = evalfv(:,:)
evec(:,:,:) = evecfv(:,:,:)
evey(:,:)   = evecsv(:,:)  

write(*,*) "av1 ",(avec(1,is),is=1,3)
write(*,*) "av2 ",(avec(2,is),is=1,3)
write(*,*) "av3 ",(avec(3,is),is=1,3)

write(*,*) "ai1 ",(ainv(1,is),is=1,3)
write(*,*) "ai2 ",(ainv(2,is),is=1,3)
write(*,*) "ai3 ",(ainv(3,is),is=1,3)

call irmm(ainv,avec,test)
write(*,*) "te1 ",(test(1,is),is=1,3)
write(*,*) "te2 ",(test(2,is),is=1,3)
write(*,*) "te3 ",(test(3,is),is=1,3)

! call irralloc(nsym,symrot,symtrn)
!
!---------------------------------------------------------------! 
!    Open IRREP.OUT for symmetry ops and irreducible reprs      ! 
!---------------------------------------------------------------! 
iout=34
if(ik.eq.1) then
  open(iout,file='IRREP.OUT',action='WRITE',form='FORMATTED') 
  call irwritedate(iout)
endif
!
!---------------------------------------------------------------! 
!   Basic info about the crystallographic symmetry operations   ! 
!---------------------------------------------------------------!  
nsym=nsymcrys
symrot(:,:,1:nsym)= symlat(:,:,lsplsymc(1:nsym))
symtrn(:,1:nsym)  = vtlsymc(:,1:nsym)
lkgrp             = (/(i,i=1,48)/)

title(1:80) = '.'
fl(1:10) = .true.
! call irinit(fl,title,avec,ainv,nsym,symrot,symtrn)
!
! properties of the symmetry operations
call irrsymops(fl,avec,ainv,nsym,symrot,symtrn,symaxc, &
                                 symsu2,iopii,iopjij)
!
! find crystallographic point group
call irrpgrpname(nsym,lkgrp,symrot,pgrpname)
!
! find classes of the point group symmetry operations
call irrclasses(nsym,lkgrp,symrot,ncls,symcls)

! write symmetry operations to IRREP.OUT
if(ik.eq.1) &
    call irrwritesym(iout,fl,title,avec,ainv,pgrpname,nsym,symrot, &
                           symtrn,symaxc,symsu2,ncls,symcls)
!---------------------------------------------------------------! 
!    Calculate the eigenfunctions, and determined the           ! 
!    irreducible representations.                               ! 
!---------------------------------------------------------------! 
!do ik=1,nkpt 
! reset variables 
ispn  = 1
nspinor=1
kcase = 0 
mspro = 0 
lkgrp = 0 
kgrpr = 0
nkgrp = 0 
irlab(1:nstfv) = ' '
! 
! properties of the k-point 
call irrkpt(iout,vkl(:,ik),vkc(:,ik),kcase,msrep,msneg,mspro, &
            nsym,symrot,symtrn) 
! 
! determine the space group G(k) of the allowed k-vector 
call irrlkgrp(vkl(1:3,ik),nsym,symrot,nkgrp,lkgrp) 
! 
! check the IRs of the k-point can be classified with pure point 
! group notations of Go(k) 
call irrlkgrpg0(vkl(1:3,ik),nkgrp,lkgrp,kgrpr,kgrpt,flg) 
! 
! find point group Go(k) 
call irrpgrpname(nkgrp,lkgrp,symrot,kgrpnam) 
! 
! find classes of the point group symmetry operations 
call irrclasses(nkgrp,lkgrp,symrot,nkgrpclass,kgrpclass) 
! 
if(flg) then 
! classify according to Go(k) 
! 
! get character table 
  call irrchartab(kgrpnam,ttir,ctir,ztir,ntir,dplph) 
! 
! check degeneracies due to time-reversal symmetry 
  call irrtimerev(vkl(1:3,ik),ztir,ctir,ttir,ntir,kgrpclass,nkgrp,lkgrp,&  
                         nsym, symrot, symtrn, & 
                         kcase,msrep,msneg,mspro,kgrpr,kgrpt,iopii,iopjij,trlab) 
! 
! calculate characters and find irreducible representation 
  call irrchar( ngkmax,nmatmax,nstfv,nspnfv,nspinor, &
                vgkl(:,:,ik,ispn), ngk(ik,ispn), ispn, symrot, symtrn, symsu2,  &
                maxdeg,vkl(1:3,ik),ttir,ctir,ztir,ntir,nkgrp,lkgrp,kgrpr,& 
                kgrpt,kgrps,ik,kgrpclass,irlab,ener,evec,toldg) 
! 
! write irreducible representations to IRREP.OUT 
  call irrwrite(iout,maxdeg,ik,irlab,kgrpnam,nkgrp,nkgrpclass,kgrpclass,kgrpt,& 
                vkl(1:3,ik),nsym,nstfv, ener, nspnfv, toldg, &
                kgrpr,ttir,ctir,ztir,ntir,lkgrp,trlab) 
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
deallocate(irlab,ener,evec,evey) 
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
