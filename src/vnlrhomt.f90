
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: vnlrhomt
! !INTERFACE:
subroutine vnlrhomt(is,wfmt1,wfmt2,zrhomt)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   is     : species number (in,integer)
!   wfmt1  : muffin-tin part of wavefunction 1 in spherical coordinates
!            (in,complex(lmmaxvr,nrcmtmax,natmtot,nspinor))
!   wfmt2  : muffin-tin part of wavefunction 2 in spherical coordinates
!            (in,complex(lmmaxvr,nrcmtmax,natmtot,nspinor))
!   zrhomt : muffin-tin charge density in spherical harmonics
!            (out,complex(lmmaxvr,nrcmtmax))
! !DESCRIPTION:
!   Calculates the complex overlap density in a single muffin-tin from two input
!   wavefunctions. See routine {\tt vnlrho}.
!
! !REVISION HISTORY:
!   Created November 2004 (Sharma)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: is
complex(8), intent(in) :: wfmt1(lmmaxvr,nrcmtmax)
complex(8), intent(in) :: wfmt2(lmmaxvr,nrcmtmax)
complex(8), intent(out) :: zrhomt(lmmaxvr,nrcmtmax)
! local variables
integer irc
! automatic arrays
complex(8) zfmt(lmmaxvr,nrcmtmax)
! generate complex density
do irc=1,nrcmt(is)
  zfmt(:,irc)=conjg(wfmt1(:,irc))*wfmt2(:,irc)
end do
! convert to spherical harmonics
call zgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,zone,zfshtvr,lmmaxvr,zfmt, &
 lmmaxvr,zzero,zrhomt,lmmaxvr)
return
end subroutine
!EOC

