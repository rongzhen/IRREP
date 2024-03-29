
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: r3mtv
! !INTERFACE:
subroutine r3mtv(a,x,y)
! !INPUT/OUTPUT PARAMETERS:
!   a : input matrix (in,real(3,3))
!   x : input vector (in,real(3))
!   y : output vector (out,real(3))
! !DESCRIPTION:
!   Multiplies the transpose of a real $3\times 3$ matrix with a vector. The
!   output vector can also be the input vector.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: a(3,3)
real(8), intent(in) :: x(3)
real(8), intent(out) :: y(3)
! local variables
real(8) z(3)
z(1)=a(1,1)*x(1)+a(2,1)*x(2)+a(3,1)*x(3)
z(2)=a(1,2)*x(1)+a(2,2)*x(2)+a(3,2)*x(3)
z(3)=a(1,3)*x(1)+a(2,3)*x(2)+a(3,3)*x(3)
! copy to output vector
y(:)=z(:)
return
end subroutine
!EOC
