!BOP
! !ROUTINE: irrwritesym
! !INTERFACE:
subroutine  irrwritesym(iout,fl,title,avec,ainv,pgrpname,nsym,symrot, &
                              symtrn,symaxc,symsu2,ncls,symcls) 
!
! !DESCRIPTION: 
!   Outputs the properties of the symmetry operations 
!   $\{{\bf R}_g|{\bf t}_g\}$ and ${\bf U}({\bf R}_g^p)$ to 
!   {\tt IRREP.OUT}. 
! 
! !REVISION HISTORY: 
!   Created  September 2007 (Clas Persson) 
!EOP 
!BOC 
implicit none 
! arguments 
integer,      intent(in)    :: iout 
logical,      intent(inout) :: fl(10)
character(80),intent(in)    :: title 
real(8),      intent(in)    :: ainv(3,3), avec(3,3)
character(10),intent(in)    :: pgrpname(2)
integer,      intent(in)    :: symrot(3,3,48),nsym 
real(8),      intent(in)    :: symtrn(3,48)  
real(8),      intent(in)    :: symaxc(3,48) 
complex(8),   intent(in)    :: symsu2(2,2,48)
character(6), intent(in)    :: symcls(48)
integer,      intent(in)    :: ncls
! local variables 
integer       i,j,is,ia,it,nc 
integer       iclas, ltest(48) 
real(8)       ds,rot(3,3),axc(3),axl(3),tra(3),tav(3) 
character(6)  cnold  
! 
!---------------------------------------------------------------! 
!        Crystal symmetry operators to IRREP.OUT                ! 
!---------------------------------------------------------------! 
write(iout,'(a80)') title
if(fl(9)) fl(3)=.true.
if(fl(2)) fl(1)=.true.
if(     fl(8)) write(iout,'(21a)',advance='NO') " Non-symmorphic crystal"
if(.not.fl(8)) write(iout,'(19a)',advance='NO') " Symmorphic crystal"
if(     fl(3)) write(iout,'(24a)',advance='NO') " with inversion symmetry,"
if(.not.fl(3)) write(iout,'(28a)',advance='NO') " without inversion symmetry,"
if(     fl(1)) write(iout,'(30a)')              " complex eigenfunctions       "
if(.not.fl(1)) write(iout,'(30a)')              " real eigenfunctions          "
if(     fl(4)) write(iout,'(16a)',advance='NO') " Spin-polarized,"
if(.not.fl(4)) write(iout,'(20a)',advance='NO') " Not spin-polarized,"
if(     fl(2)) write(iout,'(30a)')              " with spin-orbit coupling     "
if(.not.fl(2)) write(iout,'(30a)')              " without spin-orbit coupling  "
!
write(iout,'(" Point group: ",A3," ",A5)')pgrpname(1)(1:3),pgrpname(2)(1:3) 
write(iout,'(/," classes",5X,I2," symmetry operations {R|t}u(R)",$)') nsym 
ltest(1:48)=0 
j=0 
do iclas=1,ncls 
  i=1 
  do while(i.le.nsym.and.ltest(i).ne.0)  
    i=i+1 
  enddo 
  j=0 
  do is=i,nsym
    if(symcls(is).eq.symcls(i)) then 
      j=j+1 
      ltest(is)=iclas 
    endif 
  enddo 
  cnold(1:6)=symcls(i)(1:6) 
  select case(cnold(2:3))   
    case(' E')           
      write(iout,'(/,X,A6,6X," 1 unit rotation",$)') cnold 
    case(' I')          
      write(iout,'(/,X,A6,6X," 1 inversion",$)') cnold 
    case('IC') 
      read(cnold(4:4),'(I1)') it 
      write(iout,'(/,X,A6,6X,I2,I4,"-degree rotation times inversion",$)')& 
                                  cnold,j,360/it 
    case default  
      read(cnold(4:4),'(I1)') it 
      write(iout,'(/,X,A6,6X,I2,I4,"-degree rotation",$)') cnold,j,360/it       
  end select 
enddo  
! 
write(iout,'(////" n    R     t    axis           R        t   axis",15X,"u(R)")') 
write(iout,'(    "   lattice coordinates",7X,"Cartesian coordinates",8X,"spin coordinates")')  
do is=1,nsym   
  rot(1:3,1:3) = dble(symrot(1:3,1:3,is)) 
  call  irmm(rot,ainv,rot) 
  call  irmm(avec,rot,rot) 
  tra(:) =symtrn(:,is) 
  tav(1) = avec(1,1)*tra(1) +avec(1,2)*tra(2) +avec(1,3)*tra(3)
  tav(2) = avec(2,1)*tra(1) +avec(2,2)*tra(2) +avec(2,3)*tra(3)
  tav(3) = avec(3,1)*tra(1) +avec(3,2)*tra(2) +avec(3,3)*tra(3)
  tra(:) = tav(:) 
!
  axc(1:3)=0.0 
  axl(1:3)=0.0 
  if(it.gt.1) then 
    axc(1)=sin(symaxc(1,is))*cos(symaxc(2,is)) 
    axc(2)=sin(symaxc(1,is))*sin(symaxc(2,is)) 
    axc(3)=cos(symaxc(1,is)) 
    axl(1)=ainv(1,1)*axc(1)+ainv(1,2)*axc(2)+ainv(1,3)*axc(3)
    axl(2)=ainv(2,1)*axc(1)+ainv(2,2)*axc(2)+ainv(2,3)*axc(3)
    axl(3)=ainv(3,1)*axc(1)+ainv(3,2)*axc(2)+ainv(3,3)*axc(3)
    ds=sqrt(axl(1)**2+axl(2)**2+axl(3)**2) 
    axl(1:3)=axl(1:3)/ds  
  endif 
  write(iout,'(/,2X,3I2,F6.3,X,F6.2,5X,3F4.1,F6.3,X,F6.2)')  & 
       (symrot(1,j,is),j=1,3),symtrn(1,is),axl(1),(rot(1,j),j=1,3),& 
        tra(1),axc(1)  
  do i=2,3 
  write(iout,'(2X,3I2,F6.3,X,F6.2,5X,3F4.1,F6.3,X,F6.2,5X,"(",2F4.1,")(",2F4.1,")")')& 
       (symrot(i,j,is),j=1,3),symtrn(i,is),axl(i),(rot(i,j),j=1,3),& 
        tra(i),axc(i), & 
        symsu2(i-1,1,is),symsu2(i-1,2,is) 
  
  end do
!
  write(iout,'(I3,"  class:", A6,8X,$)') is, symcls(is)(1:6)
  cnold(1:6)=symcls(is)(1:6) 
  select case(cnold(2:3))   
    case(' E')           
      write(iout,'(3X," unit rotation")') 
    case(' I')          
      write(iout,'(3X," inversion")') 
    case('IC') 
      read(cnold(4:4),'(I1)') it 
      if(symaxc(3,is).eq.1) then
        write(iout,'(3X,I4,"-degree clockwise rotation, times inversion")') 360/it 
      elseif(symaxc(3,is).eq.1) then
        write(iout,'(3X,I4,"-degree rotation, times inversion")') 360/it
      else
        write(iout,'(3X,I4,"-degree counterclockwise rotation, times inversion")') 360/it  
      endif 
   case default  
      read(cnold(4:4),'(I1)') it 
      if(symaxc(3,is).eq.1) then
        write(iout,'(3X,I4,"-degree clockwise rotation")') 360/it 
      elseif(symaxc(3,is).eq.1) then
        write(iout,'(3X,I4,"-degree rotation")') 360/it 
      else
        write(iout,'(3X,I4,"-degree counterclockwise rotation")') 360/it  
      endif       
  end select 
enddo 
!
return 
end subroutine 
!EOC  
