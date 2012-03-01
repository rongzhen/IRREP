!BOP
! !ROUTINE: irrtimerev
! !INTERFACE:
subroutine irrtimerev(v,ztir,ctir,ttir,ntir,cnam,nkgrp,lkgrp,& 
                      kcase,msrep,msneg,mspro,irot,tran,iopii,iopjij,trlab) 
! !USES: 
use modmain 
! !INPUT/OUTPUT PARAMETERS: 
!    
! !DESCRIPTION: 
!   Checks if the irreducible representations are degenerate 
!   due to time reversal symmetry. 
! 
! !REVISION HISTORY: 
!   Created  September 2007 (Clas Persson) 
!EOP 
!BOC 
! 
implicit none 
! external functions 
real(8)   r3mdet 
logical   r3vint 
external  r3mdet,r3vint 
! arguments 
real(8)       tran(3,48), v(3) 
integer       ntir(2),nkgrp,irot(3,3,48),lkgrp(48) 
integer       kcase,msrep(48),msneg(48),mspro(2),iopii(48),iopjij(48,48) 
character(6)  cnam(48) 
character(80) ctir(48),ttir,trlab 
complex(8)    ztir(48,48) 
! 
! 
complex(8)    phs(48),csum,php1,phq1 
real(8)       sumt, sumi,arg,dtest,tav(3),tau(3),rot(3,3) 
real(8)       ros(3,3),ds 
integer       iwre, iwrt,iuj,iui,iordk,ir,ip,ip2 
integer       ircas(48) 
integer       i,j,i1,i2,icc,is,it,iu,itest,nik,n 
integer       idgn,ig,ir2,ip1,ig2,iq1,ip4 
logical       fg1 
! 
!--------------------------------------------------------------------------! 
!     Time-reversal symmetry can cause IRp and (IRp)* to belong to the     ! 
!     same eigenvaue, although they span different eigenspaces.            ! 
!                                                                          ! 
!     (a) If IRp  ~ (IRp)*  ~ IRq_real:  IRp is potentially real           ! 
!     (b) If IRp !~ (IRp)*            :  IRp is essentially complex        ! 
!     (c) If IRp  ~ (IRp)* !~ IRq_real:  IRp is pseudoreal                 ! 
!--------------------------------------------------------------------------! 
trlab(:)=' '
iordk = nint(dble(nsymcryss)/dble(mspro(1))) 
do ir=1,ntir(2) 
!=================================================================OBSSS 
  ircas(ir)=-1000 
  csum=cmplx(0.d0,0.d0)
  do ip=1,mspro(2) 
    ip2=iopii( msneg(ip) ) 
    ip4=0 
    do i=1,nkgrp 
      if(lkgrp(i).eq.ip2) ip4=i 
    enddo 
    if(ip4.eq.0) then  
      write(*,*)  
      write(*,'("Error(irrtimerev): not in G(k)")') 
      write(*,*)  
      stop   
    endif 
    icc=0 
    fg1=.true. 
    do while(fg1.and.icc.lt.ntir(1)) 
      icc=icc+1 
      if( cnam(ip4).eq.ttir(6*(icc-1)+10:6*icc+9)) fg1=.false. 
    end do  
    if(fg1) stop 'not G(k)' 
! 
    do i1=1,3 
      tav(i1)= symcryss(i1,1,msneg(ip))*symtran(1,msneg(ip))& 
              +symcryss(i1,2,msneg(ip))*symtran(2,msneg(ip))& 
              +symcryss(i1,3,msneg(ip))*symtran(3,msneg(ip))+symtran(i1,msneg(ip)) 
    end do 
    arg=-twopi*(v(1)*tav(1) +v(2)*tav(2) +v(3)*tav(3)) 
    csum=csum + ztir(ir,icc)*cmplx(cos(arg),sin(arg)) 
  end do 
  if(imag(csum).gt.1d-6.or.icc.gt.ntir(1)) then 
    write(*,*) 
    write(*,'("Error(irreptrsym): incorrect sum")') 
    write(*,*) 
    stop   
  elseif(abs(real(csum)-dble(iordk)).lt.1d-6)  then 
    ircas(ir) =  iordk 
  elseif(abs(real(csum)).lt.0.01) then 
    ircas(ir) =  0  
  elseif(abs(real(csum)+dble(iordk)).lt.1d-6) then 
    ircas(ir) = -iordk 
  else 
    write(*,*)  
    write(*,'("Error(irrtimerev): cannot find case (a), (b) or (c)")') 
    write(*,*) mspro(1),mspro(2),iordk,dble(nsymcryss)/dble(mspro(1)),csum,arg 
    write(*,*)  
    stop   
  endif 
end do 
! 
! if case equal 2(=b) or 3 (=c), for all pairs of IRi and IRj 
if(kcase.ge.2) then 
  idgn=0 
  iwre=0 
  do i=1,ntir(1) 
    fg1=.false.     
    do j=i+1,ntir(1) 
! if case (b) 
      if(ircas(i).eq.0.and.ircas(j).eq.0) then 
        sumt=0.d0 
        iwre=1 
        do ir=1,ntir(2) 
          iwrt=0 
          do ig=1,nkgrp 
            if(cnam(ig).eq.ttir(6*(ir-1)+10:6*ir+9)) then 
              ip1=lkgrp(ig) 
              if(kcase.eq.2) then 
                iq1=iopjij(ip1,msneg(1)) 
 
                do iui=1,nkgrp 
                  if(lkgrp(iui).eq.iq1) iq1=iui 
                enddo 
                fg1=.true. 
                ir2=0 
                do while(fg1) 
                  ir2=ir2+1 
                  if(ir2.gt.nkgrp) stop 'trsyma: cannot find class of Rneg*R*inv(Rneg)' 
                  if(cnam(iq1).eq.ttir(6*(ir2-1)+10:6*ir2+9)) fg1=.false.  
                end do 
                iq1=lkgrp(iq1) 
              else 
                iq1=ip1 
                ir2=ir 
              endif 
! sum over |X_i(Pg)*-X_j(Pg')| 
              iwrt=iwrt+1 
              arg=-twopi*( v(1)*symtran(1,ip1) & 
                          +v(2)*symtran(2,ip1) & 
                          +v(3)*symtran(3,ip1) ) 
              php1=cmplx(cos(arg),sin(arg)) 
! 
              arg=-twopi*( v(1)*symtran(1,iq1) & 
                          +v(2)*symtran(2,iq1) & 
                          +v(3)*symtran(3,iq1) ) 
              phq1=cmplx(cos(arg),sin(arg)) 
              sumt=sumt +abs(conjg(php1*ztir(i,ir )) & 
                                - (phq1*ztir(j,ir2)) ) 
            endif 
          enddo 
          if(iwrt.gt.iwre) iwre=iwrt   
        enddo 
        if(abs(sumt).lt.(nkgrp*1.0d-6)) then 
! 
          idgn=idgn+1 
          if(idgn.eq.1) write(*,'(                                  & 
             /,4X,"Degeneration due to time-reversal symmetry: ",   & 
             /,4X," ",$)') 
          iui=1 
          do while(ctir(I)(iui+1:iui+1).ne.' '.and.iui.le.7) 
            iui=iui+1 
          end do 
          iuj=1 
          do while(ctir(J)(iuj+1:iuj+1).ne.' '.and.iuj.le.7) 
            iuj=iuj+1 
          end do 
          write(*,   '(A,"+",A,";  ",$)') ctir(i)(1:iui),ctir(j)(1:iuj) 
          write(trlab,'(A,"+",A,";  ",$)') ctir(i)(1:iui),ctir(j)(1:iuj) 
        endif 
      endif 
    enddo 
  enddo 
endif 
! 
return 
end subroutine 
!EOC   
