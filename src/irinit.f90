!BOP
! !ROUTINE: irinit 
! !INTERFACE:
subroutine  irinit(fl,title,avec,ainv,nsym,symrot,symtrn,nnlo,&
                   br1,br2,db1,dr2,di1)
! !INPUT/OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!   Reads in crystal structure, the symmetry operations, and 
!   the number of local orbitals. 
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
character(80), intent(inout) :: title
real(8),       intent(out)   :: avec(3,3)
real(8),       intent(out)   :: ainv(3,3)
integer,       intent(out)   :: nsym
integer,       intent(out)   :: symrot(3,3,48)
real(8),       intent(out)   :: symtrn(3,48)
integer,       intent(out)   :: nnlo  
! local variables
integer,       allocatable   :: natom(:)
real(8),       allocatable   :: elo(:,:)

integer        i,j,ntyp,is,it,lmax2,lomax,nloat
real(8)        det,alpha(3),aa,bb,cc,tm1(3,3),tm2(3,3)
real(8)        br1(3,3),br2(3,3),dr1(3,3),dr2(3,3)
real(8)        dum,db1(3,3),di1(3,3)
character*4    latt
logical        ortho 
!---------------------------------------------------------------! 
!        Allocate elo and emist, not using param.inc            !
!---------------------------------------------------------------!
! Note, lomax should be  consistent with lapw1 and lapwso
! Note, elo(1:lomax+1), and not elo(0:lomax)
lomax=3
lmax2=10 
nloat=3 
allocate(elo(lomax+1,nloat)) 
!---------------------------------------------------------------! 
!        Input from case.struct                                 !
!---------------------------------------------------------------!
read(20,'(a80)')       title
read(20,'(a4,23x,i3)') latt,ntyp
read(20,*) 
read(20,*) aa,bb,cc,(alpha(i),i=1,3)
allocate(natom(ntyp))
do it=1,ntyp
  read(20,*) 
  read(20,'(15x,i2)') natom(it)
  do i=2,natom(it)   
    read(20,*) 
  enddo
  read(20,*) 
  read(20,*) 
  read(20,*) 
  read(20,*) 
enddo
read(20,*) nsym
do is=1,nsym
  read(20,'(2(3i2,f10.5,/),(3i2,f10.5))')((symrot(i,j,is),j=1,3),&
                                             symtrn(i,is),i=1,3)  
  read(20,*) 
enddo
!---------------------------------------------------------------! 
!        Generates lattice vectors, NOTE the units in latgen    !
!---------------------------------------------------------------!
call latgen2(latt,alpha,br1,br2,ortho)
call irminv(br1,det,dr1)
call irminv(br2,det,dr2)
call irmm(dr2,br1,db1)
call irminv(db1,det,di1)
!
! column-wise in units of bohr
avec(1,1:3) = db1(1:3,1)*aa 
avec(2,1:3) = db1(1:3,2)*bb
avec(3,1:3) = db1(1:3,3)*cc 
call irminv(avec,det,ainv(:,:))
!
! symmetry operations into lattice coordinates
do is=1,nsym
  do i=1,3
    tm1(i,1) = db1(i,1)*symtrn(1,is)*0.d0 &
             + db1(i,2)*symtrn(2,is)*0.d0 &
             + db1(i,3)*symtrn(3,is)*0.d0
  enddo
  symtrn(:,is)=tm1(:,1)
!
  call irmm(dble(symrot(:,:,is)),db1,tm1)
  call irmm(di1,tm1,tm2)
  if(irvint(tm2(1,:)).and.irvint(tm2(2,:)).and.&
                          irvint(tm2(3,:))) then
    symrot(:,:,is)=nint(tm2)     
  else
    write(*,*) 
    write(*,'("Error(irinit): R not integer matrix")') 
    write(*,*) 
  endif
enddo
fl(2)=.false. 
write(*,*) "fixa fl(2)"
!
!write(6,*) is
!write(6,'(3I4,2X,F8.4)') symrot(1,1:3,is), symtrn(1,is)
!write(6,'(3I4,2X,F8.4)') symrot(2,1:3,is), symtrn(2,is)
!write(6,'(3I4,2X,F8.4)') symrot(3,1:3,is), symtrn(3,is)
!---------------------------------------------------------------! 
!        Read header of the vector file,                        !
!        Determine number of local orbitals = nnlo              !
!---------------------------------------------------------------!
nnlo=0
do it=1,ntyp
  read(10) dum
  read(10) elo
  if(fl(2)) read(9) dum
  if(fl(2)) read(9) elo
  do j=0,lomax
    do i=1,nloat
      if(elo(j+1,i).lt.(995.d+0)) &
        nnlo=nnlo+((2*j+1))*natom(it)    
    enddo
  enddo
enddo 
!---------------------------------------------------------------! 
!        Deallocate ALL allocatable local variables             !
!---------------------------------------------------------------!
deallocate(natom,elo)
!
return
end subroutine
!EOC

