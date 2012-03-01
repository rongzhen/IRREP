!BOP
! !ROUTINE: irrclasses
! !INTERFACE: 
subroutine irrclasses(ns,lr,r,nc,cn) 
!
! !INPUT/OUTPUT PARAMETERS: 
!   ns  : number of symmetry operations      (in, integer) 
!   lr  : list of the symmetry operations    (in ,integer(48)) wwwwwwwwwwwwwwwwwwwwwwwwww
!   r   : symmetry operations                (in ,integer(3,3,48)) 
!   nc  : number of classes                  (out,integer) 
!   cn  : class name, where                  (out,character(48)) 
!             E  = unit operation 
!             I  = inversion 
!            nCx = n (360/x)-degree rotations 
!           nICx = n (360/x)-degree rotations times inversion 
!                  
! !DESCRIPTION: 
!   Generates the classes of a group of $3\times3$ rotation matrices.
!    
! !REVISION HISTORY: 
!   Created  September 2007 (Clas Persson) 
!EOP 
!BOC 
implicit none 
! external functions 
real(8)       r3mdet 
logical       r3vint 
external      r3mdet,r3vint 
! arguments
integer,      intent(in)  :: r(3,3,48),ns 
integer,      intent(in)  :: lr(48)
integer,      intent(out) :: nc
character(6), intent(out) :: cn(48)

!integer      r(3,3,48),ns
!real(8)      t(3,48)
!integer      lr(48)
!integer      nc
!character(6) cn(48)
!
! local variables 
real(8)       irotaxc(3,48) 
integer       clas(2,48),iopjij(48,48) 
character(6)  prfx  
!  
integer       is,iu,it,iv,in,n,idet,nrot(10) 
real(8)       ds,rot(3,3),ros(3,3),roi(3,3) 
real(8)       v(3),u(3) 
integer       lcls(48),eclas(48),gclas(48,48) 
logical       flg,flg2 
real(8),      parameter    :: eps=1.d-12 
! 
!---------------------------------------------------------------! 
! Determine i where (R)^i=E,  and n where Rn=Rj*Ri*inv(Ri)      ! 
!---------------------------------------------------------------! 
do is=1,ns 
! determine proper or improper rotation 
  rot(1:3,1:3)=dble(r(1:3,1:3,lr(is))) 

  write(*,*) "######## iee ",is, lr(is), ((rot(it,n),n=1,3),it=1,3)

  idet=nint(r3mdet(rot))/nint(abs(r3mdet(rot))) 
! 
! determine it, where (R)^i=E 
  ros(1:3,1:3)= 0.d0 
  ros(1,1)    = 1.d0 
  ros(2,2)    = 1.d0 
  ros(3,3)    = 1.d0 
  clas(2,is)=0 
  it=0 
  do while(clas(2,is).eq.0.and.it.le.7) 
    it=it+1 
    call r3mm(rot,ros,ros) 
    ds=abs(ros(1,1)-1)+abs(ros(1,2)  )+abs(ros(1,3)  )& 
      +abs(ros(2,1)  )+abs(ros(2,2)-1)+abs(ros(2,3)  )& 
      +abs(ros(3,1)  )+abs(ros(3,2)  )+abs(ros(3,3)-1) 
    if(ds.le.1d-8) clas(2,is)=it*idet 
  end do 
  rot(1:3,1:3)=dble(r(1:3,1:3,lr(is))) 
  do it=1,ns 
    ros(1:3,1:3)=dble(r(1:3,1:3,lr(it))) 
    call r3minv(ros,roi) 
    call r3mm(rot,roi,roi) 
    call r3mm(ros,roi,roi) 
    call findmmi(roi,dble(r(1:3,1:3,lr(1:ns))),ns,n,ds) 
    if(ds.lt.1d-6) then 
      iopjij(is,it)=n 
    else 
      write(*,*) 
      write(*,'("Error(irrclasses): no Pn=Pj*Pi*inv(Rj)")') 
      write(*,*) 
      stop 
    end if 
  end do 
end do 
!---------------------------------------------------------------! 
!     Number of classes (nc) and class number (lcls)        ! 
!---------------------------------------------------------------! 
nc=0 
lcls(1:48)=0 
eclas(1:48)=0 
gclas(1:48,1:48)=0 
do iu=1,6 
  do iv=1,-1,-2 
    do is=1,ns 
      if(clas(2,is).eq.(iu*iv)) then 
        flg=.false. 
        it=0 
        do while( (it.lt.ns).and.(.not.flg) ) 
          it=it+1 
          n=iopjij(is,it) 
          if(n.lt.is) flg=.true. 
        end do 
        if(flg) then 
          lcls(is)=lcls(n)   
          eclas(lcls(is))=eclas(lcls(is))+1 
          n=1 
          do while(gclas(lcls(is),n).ne.0) 
             n=n+1 
          end do 
          gclas(lcls(is),n)=is 
        else 
          nc=nc+1  
          lcls(is)=nc  
          eclas(nc)=1 
          gclas(nc,1)=is 
        endif 
      endif 
    enddo 
  enddo 
enddo 
!---------------------------------------------------------------! 
!                        Class names                            ! 
!---------------------------------------------------------------! 
do is=1,ns 
  select case(clas(2,is)) 
  case( 1)  
    cn(is)(1:6) ="  E   " 
  case(-1) 
    cn(is)(1:6) ="  I   " 
  case( 2, 3, 4, 6) 
    if(eclas(lcls(is)).eq.1) then 
      write(cn(is)(1:6),'("  C",I1,"  ")') clas(2,is) 
    else 
      write(cn(is)(1:6),'(" ",I1,"C",I1,"  ")') & 
                            eclas(lcls(is)),clas(2,is) 
    endif  
    do it=1,is-1 
      if(cn(is)(1:6).eq.cn(it)(1:6).and.lcls(is).ne.& 
                                              lcls(it)) then 
        if(cn(it)(5:5).eq.' ') then 
          cn(is)(5:5)='`' 
        else 
          cn(is)(5:5)='"' 
        endif 
      endif 
    enddo   
  case(-2,-3,-4,-6) 
    if(eclas(lcls(is)).eq.1) then 
      write(cn(is)(1:6),'(" IC",I1,"  ")') abs(clas(2,is)) 
    else 
      write(cn(is)(1:6),'(I1,"IC",I1,"  ")') & 
                            eclas(lcls(is)),abs(clas(2,is)) 
    endif 
    do it=1,is-1 
      if(cn(is)(1:6).eq.cn(it)(1:6).and.lcls(is).ne.& 
                                              lcls(it)) then 
        if(cn(it)(5:5).eq.' ') then 
          cn(is)(5:5)='`' 
        else 
          cn(is)(5:5)='"' 
        endif 
      endif 
    enddo 
  case default 
    write(*,*) 
    write(*,'("Error(irrclasses): cannot find class name")') 
    write(*,*) 
    stop 
  end select 
enddo  
!
return 
end subroutine 
!EOC 
