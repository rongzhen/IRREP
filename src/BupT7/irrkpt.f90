!BOP 
! !ROUTINE: irrkpt 
! !INTERFACE: 
subroutine irrkpt(iout,v,vc,kcase,msrep,msneg,mspro, &
                  nsym,symrot,symtrn  )  
!
! !INPUT/OUTPUT PARAMETERS: 
!   ik     : input  k-vector                        (in,integer) 
! !DESCRIPTION: 
!   Determines properties of the {\bf k}-point. 
! 
! !REVISION HISTORY: 
!   Created  September 2007 (Clas Persson) 
!EOP 
!BOC 
! 
implicit none 
integer,      intent(in)  :: iout
real(8),      intent(in)  :: v(3),vc(3)
integer,      intent(in)  :: symrot(3,3,48),nsym
real(8),      intent(in)  :: symtrn(3,48)
real(8)       avkl,dtest,dmini,roti(3,3),rkl(3) 
integer       kequ(3) 
integer       msrep(48),msneg(48),mspro(2) 
integer       is,it,iu,itest,nik 
integer       i1,ii1,ii2,i,ix1,iy1,iz1,ik,ik2 
integer       kcase,ikn(3) 
logical       fg1,fg2,fg3 
real(8)       kstl(3,48) 
! 
msrep(1:48) = 0 
msneg(1:48) = 0 
mspro(1:2)  = 0 
kcase = 0 
nik   = 1 
kstl(1:3,nik) = v(1:3) 
avkl  = sqrt(v(1)**2 +v(2)**2 +v(3)**2) 
fg1=.false. 
fg2=.false. 
ii1=0 
ii2=0 
do i=1,nsym 
  call irminv(dble(symrot(1:3,1:3,i)),dtest,roti) 
  dtest = abs(roti(1,1)-1)+abs(roti(1,2))  +abs(roti(1,3)) & 
        + abs(roti(2,1))  +abs(roti(2,2)-1)+abs(roti(2,3)) & 
        + abs(roti(3,1))  +abs(roti(3,2))  +abs(roti(3,3)-1) 
  if(dtest.le.1e-6) msrep(1)=i 
!
! k*inv(Ri) 
  rkl(1) = roti(1,1)*v(1) +roti(2,1)*v(2) +roti(3,1)*v(3)
  rkl(2) = roti(1,2)*v(1) +roti(2,2)*v(2) +roti(3,2)*v(3)
  rkl(3) = roti(1,3)*v(1) +roti(2,3)*v(2) +roti(3,3)*v(3)
!
! test for k*inv(Ri)=k+K  
  dtest=abs( dble(nint(rkl(1)-v(1)))-(rkl(1)-v(1)) )& 
       +abs( dble(nint(rkl(2)-v(2)))-(rkl(2)-v(2)) )& 
       +abs( dble(nint(rkl(3)-v(3)))-(rkl(3)-v(3)) )  
  if(dtest.le.(avkl*1e-6)) kcase=1  
!
! test for k*inv(Ri)=-k+K  
  dtest=abs( dble(nint(rkl(1)+v(1)))-(rkl(1)+v(1)) )& 
       +abs( dble(nint(rkl(2)+v(2)))-(rkl(2)+v(2)) )& 
       +abs( dble(nint(rkl(3)+v(3)))-(rkl(3)+v(3)) )  
  if(dtest.le.(avkl*1e-6)) then 
    ii2=ii2+1 
    msneg(ii2)=I 
    fg2=.true.            
  endif   
!
! check if k-star point 
  fg3=.true. 
  do it=1,nik 
    dtest=abs( dble(nint(rkl(1)-kstl(1,it)))-(rkl(1)-kstl(1,it)) )& 
         +abs( dble(nint(rkl(2)-kstl(2,it)))-(rkl(2)-kstl(2,it)) )& 
         +abs( dble(nint(rkl(3)-kstl(3,it)))-(rkl(3)-kstl(3,it)) ) 
    if(dtest.le.(avkl*1e-6)) fg3=.false. 
  end do 
  if(fg3) then 
    nik=nik+1 
    msrep(nik)=i 
    kstl(1:3,nik) = rkl(1:3)  
  endif 
end do 
if(msrep(1).eq.0) STOP 'kstar: cannot find E' 
if(fg2) kcase = 2 
mspro(1) = nik 
mspro(2) = ii2 
! check if -k is equivalent to k 
dtest=abs( dble(nint(v(1)+v(1)))-(v(1)+v(1)) )& 
     +abs( dble(nint(v(2)+v(2)))-(v(2)+v(2)) )& 
     +abs( dble(nint(v(3)+v(3)))-(v(3)+v(3)) ) 
if(dtest.le.(avkl*1e-6)) then 
  kcase   = 3 
  kequ(1) = nint(v(1)+v(1)) 
  kequ(2) = nint(v(2)+v(2)) 
  kequ(3) = nint(v(3)+v(3)) 
endif 
! 
! output 
write(iout,'(//," k-point number",i6)') ik 
write(iout,'(" k=      (",3f7.4,")")') (v(it),it=1,3) 
write(iout,'("         (",3f7.4,") in Cartesian coordinates")') & 
               (vc(it),it=1,3) 
write(iout,'(" k-star: ",$)') 
if(kcase.ne.0) then 
  do is=1,nik   
    write(iout,'("(",3f7.4,") ",$)') (kstl(it,is),it=1,3) 
    if(mod(is,2).eq.0.and.is.ne.nik) write(iout,'(/,8x," ",$)')  
  end do 
endif 
select case(kcase) 
case(1) 
  write(iout,'(/,10x,"-k does not belong to star of k")')    
case(2)  
  write(iout,'(/,10x,"-k belong to star of k, using sym.op. no. ",i2)') msneg(1)  
case(3)  
  write(iout,'(/,10x,"-k is equivalent to k, with K=(",3I2,")") ') (kequ(it),it=1,3) 
case default 
  write(iout,'(/,10x,"k is a general k-point")') 
end select 
write(iout,*)  
return 
end subroutine 
!EOC  
