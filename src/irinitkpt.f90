!BOP
! !ROUTINE: irinitkpt 
! !INTERFACE:
subroutine  irinitkpt(fl,ik,nnlo,avec,ainv,ngk,nne,vl,vc,&
                       nstfv,nspnfv,nmatmax,kkl,en,ev, &
                       br1,br2,db1,dr2,di1)

                                               


! !INPUT/OUTPUT PARAMETERS:
!   fl   : input  flag                              (in, logical(10))
!   ik   : input  k-point number                    (in,integer)
!   nnlo : input  number of local orbital           (in,integer)
!   avec : input  crystal lattice vectors           (in,real(3,3))
!   ainv : input  inverted vectors                  (in,real(3,3))
!   vl   : output k in lattice coordinates          (out, real(3))
!   vc   : output k in Cartesian coordinates        (out,real(3))
!   kkl  : output k+K in lattice coordinates        (out, real(3,:))
!   ener : output eigenvalues                       (out,real(:,:))
!   evec : output eigenvectors (plane-wave coeff.)  (out,complex(:,:,:))
!
! !DESCRIPTION:
!   Reads in crystal k-vector, eigenvalues, and eigenvectors. 
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
logical,       intent(in)    ::  fl(10)
integer,       intent(in)    ::  nnlo,      nstfv,nspnfv,nmatmax
real(8),       intent(in)    ::  avec(3,3)
real(8),       intent(in)    ::  ainv(3,3)
integer,       intent(inout) ::  ik
integer,       intent(inout) ::  ngk(100,nspnfv),nne
real(8),       intent(inout) ::  vl(3,100)
real(8),       intent(inout) ::  vc(3,100)
real(8),       intent(inout) ::  kkl(3,nmatmax,100,nspnfv)
real(8),       intent(inout) ::  en(nstfv,nspnfv)
complex(8),    intent(inout) ::  ev(nmatmax,nstfv,nspnfv)

real(8)        br1(3,3),br2(3,3),dr1(3,3),dr2(3,3)
real(8)        dum,db1(3,3),di1(3,3)


!real(8),       allocatable   ::  kkl(:,:)
!real(8),       allocatable   ::  en(:)
!complex(8),    allocatable   ::  ev(:,:)
!
! local variables
real(8),       allocatable    ::  a(:)
integer,       allocatable    ::  ikv(:,:)
integer        i,in,j,nv,ne,num,ispn
real(8)        det,sk(3),wght,tm(3)
character(10)  kname
!---------------------------------------------------------------! 
!        Input from case.vectors[up/dn][so]                     !
!---------------------------------------------------------------!
nne=0
read(10,end=999)          (sk(i),i=1,3),kname,nv,nne,wght

write(*,*)  (sk(i),i=1,3),kname,nv,nne,wght

ispn=1

! if(fl(2)) read(9,end=999) (sk(i),i=1,3),kname,nv,nne,wght
!if(allocated(kkl)) deallocate(kkl)
!if(allocated(ev))  deallocate(ev)
!if(allocated(en))  deallocate(en)
!allocate(kkl(3,nv),ev(nv,nne),en(nne))
 
  
allocate(a(nv),ikv(3,nv))
!
read(10) ((ikv(j,in),j=1,3),in=1,nv)
do in=1,nv
  kkl(1:3,in,ik,1) = dble(ikv(1:3,in))
enddo
if(fl(2)) read(9) ((ikv(j,in),j=1,3),in=1,nv)
!
do j=1,nne
  read(10) num,en(j,ispn)
!  if(fl(2)) read(9) num,en(j,ispn)
!  do in=1,nv
!    b(i,nne)=cmplx(0.d0,0.d0)
!  enddo
  if(fl(1)) then
    read(10) (ev(in,j,ispn),in=1,nv)
    if(fl(2)) then
      write(*,*) "irinitkpt: not so"
      stop
!      read(9)  (b(in,j),in=1,nv)
    endif
  else
    read(10) (a(in),in=1,nv)
    do in=1,nv
      ev(in,j,ispn)=cmplx(a(in),0.d0)
    enddo
!    if(fl(2)) then
!      read(9) (a(in),in=1,nv)
!      do  in=1,nv
!              B(In,NNE)=CMPLX(APA(In),0.D0)
!      enddo  
!    endif
  endif 
  write(*,*) "j,ik= ",j,ik,nv,nnlo,ngk(ik,ispn)
enddo

! -------- OOBBSSSS, anvaend bara avec, ej ainv ------------
!

!---------------------------------------------------------------! 
!        Change to primitive coordinates for kkl and vl         !
!        Also, K => (K+k)                                       !
!---------------------------------------------------------------!
do j=1,3
  tm(j)=db1(j,1)*sk(1)+db1(j,2)*sk(2)+db1(j,3)*sk(3) 
enddo
vl(1:3,ik)=tm(:)

write(*,*) "######### ", (vl(j,ik),j=1,3)

do in=1,nv
  do j=1,3
    tm(j)=db1(j,1)*kkl(1,in,ik,1) & 
         +db1(j,2)*kkl(2,in,ik,1) & 
         +db1(j,3)*kkl(3,in,ik,1) 
  enddo
  kkl(1:3,in,ik,1)=tm(1:3)  
  do j=1,3
    kkl(j,in,ik,1) = kkl(j,in,ik,1)+vl(j,ik)
  enddo
enddo
!
! Cartesian coordinates for vc
do j=1,3
  vc(j,ik)=di1(j,1)*vl(1,ik)+di1(j,2)*vl(2,ik)&
                            +di1(j,3)*vl(3,ik)
enddo
! do not use coefficients of the local orbitals  
ngk(ik,ispn)=nv-nnlo
!
deallocate(a,ikv)
write(*,*) "######### vl ", (vl(j,ik),j=1,3)
write(*,*) "######### vc ", (vc(j,ik),j=1,3)

!do i=1,10
!write(*,*) "k ",nv, ngk(ik,ispn), (kkl(j,i,ik,ispn),j=1,3)  
!enddo


!
return
999 ik=-ik
end subroutine
!EOC

