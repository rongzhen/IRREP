!BOP
! !ROUTINE: irminv
! !INTERFACE:
subroutine  irminv(m,d,im)
! !INPUT/OUTPUT PARAMETERS:
!   m  : input  matrix                    (in, real(3,3))
!   d  : output determinant of m          (out,real)
!   im : output inverted matrix           (out,real(3,3))
! 
! !DESCRIPTION:
!   Determinant and inverse of the $3\times3$ matrix. 
!
! !REVISION HISTORY:
!   Created  September 2007 (Clas Persson)
!EOP
!BOC
implicit none
! arguments
real(8),      intent(in)   :: m(3,3)
real(8),      intent(out)  :: d 
real(8),      intent(out)  :: im(3,3)
! local variables
real(8)       tm(3,3)
integer       i,j
!
d = m(1,1)*m(2,2)*m(3,3) + m(1,2)*m(2,3)*m(3,1) &
  + m(1,3)*m(2,1)*m(3,2) - m(1,1)*m(2,3)*m(3,2) &
  - m(1,2)*m(2,1)*m(3,3) - m(1,3)*m(2,2)*m(3,1)
if(dabs(d).lt.1.0E-16) then
  write(*,*)
  write(*,'("Error(irminv): |m| = 0")')
  write(*,*)
  stop
endif
!
tm(1,1)= (m(2,2)*m(3,3)-m(2,3)*m(3,2))
tm(1,2)=-(m(2,1)*m(3,3)-m(2,3)*m(3,1))
tm(1,3)= (m(2,1)*m(3,2)-m(2,2)*m(3,1))
tm(2,1)=-(m(1,2)*m(3,3)-m(1,3)*m(3,2))
tm(2,2)= (m(1,1)*m(3,3)-m(1,3)*m(3,1))
tm(2,3)=-(m(1,1)*m(3,2)-m(1,2)*m(3,1))
tm(3,1)= (m(1,2)*m(2,3)-m(1,3)*m(2,2))
tm(3,2)=-(m(1,1)*m(2,3)-m(1,3)*m(2,1))
tm(3,3)= (m(1,1)*m(2,2)-m(1,2)*m(2,1))
!
do i=1,3
do j=1,3     
  im(i,j)=tm(j,i)/d  
enddo
enddo
return
end subroutine
!EOC
