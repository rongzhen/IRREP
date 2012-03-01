!BOP 
! !ROUTINE: irrrotaxis 
! !INTERFACE:  
subroutine irrrotaxis(a,x) 
! !INPUT/OUTPUT PARAMETERS: 
!   a : input  matrix (in, real(3,3)) 
!   x : output angles (out,real(3)) 
!      
! !DESCRIPTION: 
!   Finds the normalized rotation axis $\{\theta,\phi\}$ of 
!   the proper rotation matrix {\bf a}. Angles are defined as 
!   $0 \le x_1$$\equiv$$\theta \le\pi$ and 
!   $0 \le x_2$$\equiv$$\phi\le2\pi$. 
!   $x_3$ equals 1 for clockwise rotation, -1 for 
!   counterclockwise rotation, and 0 for unity rotation 
!   or 180 degree rotation. 
! 
! !REVISION HISTORY: 
!   Created  September 2007 (Clas Persson) 
!EOP 
!BOC 
implicit none 
! external functions 
real(8)   r3dist,r3mdet 
external  r3dist,r3mdet 
! arguments 
real(8), intent(in)  :: a(3,3) 
real(8), intent(out) :: x(3) 
! local variables 
integer, parameter   :: itmx=1000 
real(8), parameter   :: pi=3.1415926535897932385d0 
integer  i,j,it 
real(8)  b(3,3),an(2),bn(2),v(3) 
real(8)  fc,ds,dmin,r(3),ro(3) 
!  
if(abs(r3mdet(a)-1.0).gt.1e-6) then  
  write(*,*) 
  write(*,'("Error(irrrotaxis): has to be proper rotation")') 
  write(*,*) 
  stop 
endif 
! find rotation axis (can be done more elegantly) 
! for all but 0 degree rotations 
x(1:3)=0.0 
ds=abs(a(1,1)-1.)+abs(a(1,2)   )+abs(a(1,3)   )& 
  +abs(a(2,1)   )+abs(a(2,2)-1.)+abs(a(2,3)   )& 
  +abs(a(3,1)   )+abs(a(3,2)   )+abs(a(3,3)-1.) 
if(ds.gt.1d-6) then 
  x(1)=pi*0.25 
  x(2)=pi 
  dmin=98.7654321 
  it  =0 
  fc  =1.0 
  do while(it.le.itmx.and.dmin.ge.1d-9) 
    it=it+1 
    an(1)=x(1) 
    an(2)=x(2) 
    do i =-2,2 
      do j =-2,2 
        bn(1)=an(1)+0.25*pi*fc*dble(i)/2.0 
        bn(2)=an(2)+     pi*fc*dble(j)/2.0 
        v(1) =sin(bn(1))*cos(bn(2)) 
        v(2) =sin(bn(1))*sin(bn(2)) 
        v(3) =cos(bn(1)) 
        call r3mv(a,v,r) 
        ds=r3dist(v,r) 
        if(ds.lt.dmin) then 
          dmin=ds 
          x(1)=bn(1) 
          x(2)=bn(2) 
        endif 
      end do 
    end do 
    fc=fc/2.0 
  end do 
  if(abs(x(1)).le.1d-7) x(2)=0.0 
  if(it.ge.itmx) then  
    write(*,*) 
    write(*,'("Error(irrrotaxis): cannot find angles")') 
    write(*,*) 
    stop 
  endif 
! 
! clockwise or counterclockwise, 
! all but for 0 and 180 degree rotations 
  call r3mm(a,a,b) 
  ds=abs(b(1,1)-1.)+abs(b(1,2)   )+abs(b(1,3)   )& 
    +abs(b(2,1)   )+abs(b(2,2)-1.)+abs(b(2,3)   )& 
    +abs(b(3,1)   )+abs(b(3,2)   )+abs(b(3,3)-1.) 
  if(ds.gt.1d-6) then 
    ro(1)=sin(x(1))*cos(x(2)) 
    ro(2)=sin(x(1))*sin(x(2)) 
    ro(3)=cos(x(1)) 
    dmin=0.0 
    do it=1,3 
      v(1:3)=0.0 
      v(it) =1.0 
      call r3mv(a,v,r) 
      call r3cross(v,r,r) 
      dmin=dmin+ro(1)*r(1)+ro(2)*r(2)+ro(3)*r(3) 
    end do 
    if(dmin.gt.1e-6) then 
      x(3)= 1.0 
    elseif(dmin.lt.-1e-6) then 
      x(3)=-1.0 
    else 
      write(*,*) 
      write(*,'("Error(irrrotaxis): cannot find rotation direction")') 
      write(*,*)  
!      stop 
    endif 
  endif 
endif 
return 
end subroutine 
!EOC 
