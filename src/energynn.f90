
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine energynn
use modmain
implicit none
! local variables
integer is,ia,ias,lmax
real(8) t1
complex(8) zrho0
! allocatable arrays
real(8), allocatable :: jlgr(:,:,:)
complex(8), allocatable :: zpchg(:)
complex(8), allocatable :: zrhomt(:,:,:)
complex(8), allocatable :: zrhoir(:)
complex(8), allocatable :: zvclmt(:,:,:)
complex(8), allocatable :: zvclir(:)
allocate(jlgr(0:lmaxvr+npsden+1,ngvec,nspecies))
allocate(zpchg(natmtot))
allocate(zrhomt(lmmaxvr,nrmtmax,natmtot))
allocate(zrhoir(ngrtot))
allocate(zvclmt(lmmaxvr,nrmtmax,natmtot))
allocate(zvclir(ngrtot))
! set the density to zero
zrhomt(:,:,:)=0.d0
zrhoir(:)=0.d0
! set up the point charge array
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    zpchg(ias)=spzn(is)
  end do
end do
! compute the required spherical Bessel functions
lmax=lmaxvr+npsden+1
call genjlgpr(lmax,gc,jlgr)
! solve the complex Poisson's equation
call zpotcoul(nrmt,nrmtmax,spnrmax,spr,1,gc,jlgr,ylmg,sfacg,zpchg,zrhomt, &
 zrhoir,zvclmt,zvclir,zrho0)
! compute the nuclear-nuclear energy
engynn=0.d0
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    t1=dble(zvclmt(1,1,ias))*y00-spzn(is)/spr(1,is)
    engynn=engynn+spzn(is)*t1
  end do
end do
engynn=0.5d0*engynn
deallocate(jlgr,zpchg,zrhomt,zrhoir,zvclmt,zvclir)
return
end subroutine

