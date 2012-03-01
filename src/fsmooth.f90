
! Copyright (C) 2002-2005 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: fsmooth
subroutine fsmooth(m,n,ld,f)
! !INPUT/OUTPUT PARAMETERS:
!   m  : number of 3-point running averages to perform (in,integer)
!   n  : number of point (in,integer)
!   ld : leading dimension (in,integer)
!   f  : function array (inout,real(ld,n))
! !DESCRIPTION:
!   Removes numerical noise from a function by performing $m$ successive
!   3-point running averages on the data. The endpoints are kept fixed.
!
! !REVISION HISTORY:
!   Created December 2005 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: m
integer, intent(in) :: n
integer, intent(in) :: ld
real(8), intent(inout) :: f(ld,n)
! local variables
integer i,j
do i=1,m
  do j=2,n-1
    f(1,j)=0.3333333333333333333d0*(f(1,j-1)+f(1,j)+f(1,j+1))
  end do
end do
return
end subroutine
!EOC

