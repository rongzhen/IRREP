!BOP
! !ROUTINE: irsu2
! !INTERFACE:
subroutine  irsu2(r,a,s)
! !INPUT/OUTPUT PARAMETERS:
!   r : input  rotation matrix  (in, real(3,3))
!   a : Euler's angles          (out, real(3))
!   s : output SU(2) matrix     (out,complex(2,2))
!
! !DESCRIPTION:
!   Calculates Euler's angles and the SU(2) matrix. The matrix 
!   is defined to be positive: $\rm{Re}|\rm{s}_{11}| \ge 0$.
!   See Sect.~\ref{s:DG} for definition of rotation angles.
!
! !REVISION HISTORY:
!   Created  September 2007 (Clas Persson)
!EOP
!BOC
implicit none
! arguments
real(8),    intent(in)  :: r(3,3)
real(8),    intent(out) :: a(3)
complex(8), intent(out) :: s(2,2)
! local variables
integer     i1,i2
real(8)     pi,amax
complex(8), parameter   :: zi=(0.d0,1.d0)
!
! Euler's angles
pi=dacos(-1.d0)
a(1)=dacos((r(3,3)))
do while(a(1).gt.pi.or.a(1).lt.0)  
   a(1)=a(1)-sign(pi,a(1))
enddo
!
if(abs(dabs(r(3,3))-1.d0).lt.5e-6) then
  a(2)=dacos( r(3,3)*r(1,1) )
  if(r(1,2).lt.0.d0)              a(2)=2*pi-a(2)
  a(3)=0.d0
else
  a(2)=dacos( -r(1,3)/dsin(a(1)) )
  if((r(2,3)/dsin(a(1))).lt.0.d0) a(2)=2*pi-a(2)
  a(3)=dacos( -r(3,1)/dsin(a(1)) )
  if((r(3,2)/dsin(a(1))).lt.0.d0) a(3)=2*pi-a(3)
endif
!
!.....0<= a(1) <=pi,   0<= a(2), a(3) <=2*pi     
do i1=1,3
  do i2=-48,48
    if(dabs(  12*a(i1)-i2*pi     ).lt.5e-6) a(i1)=i2*pi/12.d0
  enddo
  amax=pi
  if(i1.eq.2) amax=2*pi
  do while(a(i1).gt.(amax-5e-6).or.a(i1).lt.0.d0)
    a(i1)=a(i1)-sign(amax,a(i1))
  enddo
  if(dabs(a(i1)-0.d0).lt.5e-6) a(i1)=0.d0
  if(dabs(a(i1)-amax).lt.5e-6) a(i1)=amax
enddo   
!
! SU(2)
s(1,1)= dcos(a(1)*0.5)*exp( zi*(a(2)+a(3))*0.5)
s(1,2)= dsin(a(1)*0.5)*exp( zi*(a(2)-a(3))*0.5)
s(2,1)=-dsin(a(1)*0.5)*exp(-zi*(a(2)-a(3))*0.5)
s(2,2)= dcos(a(1)*0.5)*exp(-zi*(a(2)+a(3))*0.5)
!
! use positive matrices
if(real(s(1,1)).lt.0.d0) s(1:2,1:2) = -s(1:2,1:2)
!
return
end subroutine
!EOC

