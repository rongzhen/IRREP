!BOP 
! !ROUTINE: irrpgrname 
! !INTERFACE:  
subroutine irrpgrpname(n,lk,r,cn) 
! !INPUT/OUTPUT PARAMETERS: 
!   n  : number of crystal symmetry operations      (in, integer) 
!   lk : output list of k-group symmetry operations (out,integer(48))
!   r  : symmetry rotations Rg                      (in ,integer(3,3,48)) 
!   cn : point group name where columns             (out,character(2))  
!            1 = Schoenflis notation 
!            2 = International notation  
!    
! !DESCRIPTION: 
!   Finds the point group name. 
! 
! !REVISION HISTORY: 
!   Created  September 2007 (Clas Persson) 
!EOP 
!BOC 
implicit none 
! external functions 
! arguments 
integer        n,lk(48) 
integer        r(3,3,48) 
character(10)  cn(2) 
! local variables 
integer        is,it,ni,nrot(10) 
real(8)        det,ds,rot(3,3),ros(3,3) 
character(11)  bte 
! 
nrot(1:10)=0 
do is=1,n  
  rot(1:3,1:3)=dble(r(1:3,1:3,lk(is))) 
  call irminv(rot,det,ros)
  ros(1:3,1:3)=0.0 
  ros(1,1)    =1.0 
  ros(2,2)    =1.0 
  ros(3,3)    =1.0 
  ds=9.0 
  it=0 
  do while(ds.ge.1d-8.and.it.le.7) 
    it=it+1 
    call irmm(rot,ros,ros) 
    ds=abs(ros(1,1)-1.)+abs(ros(1,2)   )+abs(ros(1,3)   )& 
      +abs(ros(2,1)   )+abs(ros(2,2)-1.)+abs(ros(2,3)   )& 
      +abs(ros(3,1)   )+abs(ros(3,2)   )+abs(ros(3,3)-1.) 
  end do  
  ni=abs(it)+5*nint( 0.5*(1-abs(det)/det) )    
  if(abs(it).eq.6) ni=ni-1 
  nrot(ni)=nrot(ni)+1 
end do 
write(bte,'(5i1,X,5i1)') (nrot(it),it=1,10) 
! 
!-----------------------------------------------------------------! 
!     Identify the point group: for example bte='13000 02020'     ! 
!     means the identity operator, three 180-degree rotatis,      ! 
!     two 180-degree rotations times inversion, two 90-degree     ! 
!     rotations times invesion.                                   ! 
!-----------------------------------------------------------------! 
cn(1:2)(1:10)=' ' 
select case(bte) 
case('10000 00000')  
  cn(1:2)(1:5) = (/'C1   ','1    '/) 
case('10000 10000')   
  cn(1:2)(1:5) = (/'Ci  ' ,'-1   '/) 
case('11000 00000')   
  cn(1:2)(1:5) = (/'C2   ','2    '/) 
case('10000 01000')   
  cn(1:2)(1:5) = (/'Cs   ','m    '/) 
case('11000 11000')   
  cn(1:2)(1:5) = (/'C2h  ','2/m  '/)   
case('13000 00000')   
  cn(1:2)(1:5) = (/'D2   ','222  '/) 
case('11000 02000')   
  cn(1:2)(1:5) = (/'C2v  ','mm2  '/)      
case('13000 13000')   
  cn(1:2)(1:5) = (/'D2h  ','mmm  '/) 
case('11020 00000')   
  cn(1:2)(1:5) = (/'C4   ','4    '/) 
case('11000 00020')   
  cn(1:2)(1:5) = (/'S4   ','-4   '/) 
case('11020 11020')   
  cn(1:2)(1:5) = (/'C4h  ','4/m  '/) 
case('15020 00000')   
  cn(1:2)(1:5) = (/'D4   ','422  '/) 
case('11020 04000')   
  cn(1:2)(1:5) = (/'C4v  ','4mm  '/) 
case('13000 02020')   
  cn(1:2)(1:5) = (/'D2d  ','-42m '/) 
case('15020 15020')   
  cn(1:2)(1:5) = (/'D4h  ','4/mmm'/) 
case('10200 00000')   
  cn(1:2)(1:5) = (/'C3   ','3    '/) 
case('10200 10200')   
  cn(1:2)(1:5) = (/'C3i  ','-3   '/) 
case('13200 00000')   
  cn(1:2)(1:5) = (/'D3   ','32   '/) 
case('10200 03000')   
  cn(1:2)(1:5) = (/'C3v  ','3m   '/) 
case('13200 13200')   
  cn(1:2)(1:5) = (/'D3d ','-3m   '/) 
case('11202 00000')   
  cn(1:2)(1:5) = (/'C6   ','6    '/) 
case('10200 01002')   
  cn(1:2)(1:5) = (/'C3h  ','-6   '/) 
case('11202 11202')   
  cn(1:2)(1:5) = (/'C6h  ','6/m  '/) 
case('17202 00000')   
  cn(1:2)(1:5) = (/'D6   ','622  '/) 
case('11202 06000')   
  cn(1:2)(1:5) = (/'C6v  ','6mm  '/) 
case('13200 04002')   
  cn(1:2)(1:5) = (/'D3h  ','-6m2 '/) 
case('17202 17202')   
  cn(1:2)(1:5) = (/'D6h  ','6/mmm'/) 
case('13800 00000')   
  cn(1:2)(1:5) = (/'T    ','23   '/) 
case('13800 13800')   
  cn(1:2)(1:5) = (/'Th   ','m-3  '/) 
case('19860 00000')   
  cn(1:2)(1:5) = (/'O    ','432  '/) 
case('13800 06060')   
  cn(1:2)(1:5) = (/'Td   ','-43m '/) 
case('19860 19860')   
  cn(1:2)(1:5) = (/'Oh   ','m-3m '/) 
case default 
  write(*,*) 
  write(*,'("Error(irrpgrpname): not found point group")') 
  write(*,*) 
  stop  
end select   
return 
end subroutine 
!EOC 
