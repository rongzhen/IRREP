!BOP 
! !ROUTINE: irrchartab  
! !INTERFACE:
subroutine irrchartab(name,ttir,ctir,ztir,ntir,dph) 
! !INPUT/OUTPUT PARAMETERS: 
!   name : input  point group name               (in, character) 
!   f    : output flag,                          (out,logical) 
! 
! !DESCRIPTION:  
!   Generates the point-group character table for the given
!   point group name. 
! 
! !REVISION HISTORY: 
!   Created  September 2007 (Clas Persson) 
!EOP 
!BOC 
implicit none 
! arguments 
character(10) name(2) 
real(8)       v(3),t(3,48) 
integer       n,r(3,3,48) 
logical       f 
! initialise universal variables 
! local variables 
character(80) ctir(48),ttir 
integer       i,is,it,i1,i2,iu,j1,j2,idi,ntir(2) 
real(8)       ri(3,3),vr(3),arg,dtest,ti(3),tr(3) 
real(8)       dph 
complex(8)    ztir(48,48) 
real(8),      parameter :: twopi=6.2831853071795864769d0 
! 
!----------------------------------------------------------! 
!     identify the point group and the character table     ! 
!----------------------------------------------------------! 
ctir(1:48)(1:80)=' ' 
select case(name(1)(1:3)) 
case("C1 ")    
      ttir     ='           E  '  
      ctir( 1) ='G1   A1    1  ' 
      ctir( 2) ='G2   A1/2  1  ' 
case("Ci ")
! #            
case("C2 ","Cs ") 
      ttir     ='           E     C2 '  
      if(name(1)(1:3).eq."Cs ") & 
      ttir     ='           E    IC2 '  
      ctir( 1) ='G1   A     1     1  ' 
      ctir( 2) ='G2   B     1    -1  ' 
      ctir( 3) ='G3  1E1/2  1     i  ' 
      ctir( 4) ='G4  2E1/2  1    -i  ' 
case("C2h")
! #
case("D2 ","C2v")
      ttir     ='           E     C2    C2`   C2"' 
      if(name(1)(1:3).eq."C2v") & 
      ttir     ='           E     C2   IC2   IC2`'
      ctir( 1) ='G1   A1    1     1     1     1  ' 
      ctir( 2) ='G2   B3    1    -1     1    -1  ' 
      ctir( 3) ='G3   B1    1     1    -1    -1  ' 
      ctir( 4) ='G4   B2    1    -1    -1     1  ' 
      ctir( 5) ='G5   E1/2  2     0     0     0  ' 
case("D2h ")
      ttir     ='           E     C2    C2`   C2"   I    IC2   IC2`  IC2"'                         
      ctir( 1) ='G1+  A1g   1     1     1     1     1     1     1     1  '
      ctir( 2) ='G2+  B3g   1    -1     1    -1     1    -1     1    -1  '
      ctir( 3) ='G3+  B1g   1     1    -1    -1     1     1    -1    -1  '
      ctir( 4) ='G4+  B2g   1    -1    -1     1     1    -1    -1     1  '
      ctir( 5) ='G1-  A1u   1     1     1     1    -1    -1    -1    -1  '
      ctir( 6) ='G2-  B3u   1    -1     1    -1    -1     1    -1     1  '
      ctir( 7) ='G3-  B1u   1     1    -1    -1    -1    -1     1     1  '
      ctir( 8) ='G4-  B2u   1    -1    -1     1    -1     1     1    -1  '
      ctir( 9) ='G5+  E1/2g 2     0     0     0     2     0     0     0  '
      ctir(10) ='G5-  E1/2u 2     0     0     0    -2     0     0     0  '
case("C4 ","S4 ")  
      ttir     ='           E     C4    C2    C4-' 
      if(name(1)(1:3).eq."S4 ") &
      ttir     ='           E    IC4    C2   IC4-'
      ctir( 1) ='G1   A     1     1     1     1  ' 
      ctir( 2) ='G2   B     1    -1     1    -1  ' 
      ctir( 3) ='G3  2E     1     i    -1    -i  ' 
      ctir( 4) ='G4  1E     1    -i    -1     i  ' 
      ctir( 5) ='G5  2E1/2  1     d     i     d* ' 
      ctir( 6) ='G6  1E1/2  1     d*   -i     d  ' 
      ctir( 7) ='G7  2E3/2  1    -d     i    -d* ' 
      ctir( 8) ='G8  1E3/2  1    -d*   -i    -d  ' 
case("C4h")
! #
case("D4 ","C4v","D2d")  
      ttir     ='           E    2C4    C2   2C2   2C2 ' 
      if(name(1)(1:3).eq."C4v") &
      ttir     ='           E    2C4    C2  2IC2` 2IC2"'
      if(name(1)(1:3).eq."D2d") &
      ttir     ='           E   2IC4    C2   2C2  2IC2 '
      ctir( 1) ='G1   A1    1     1     1     1     1  ' 
      ctir( 2) ='G2   A2    1     1     1    -1    -1  ' 
      ctir( 3) ='G3   B1    1    -1     1     1    -1  ' 
      ctir( 4) ='G4   B2    1    -1     1    -1     1  ' 
      ctir( 5) ='G5   E     2     0    -2     0     0  ' 
      ctir( 6) ='G6   E1/2  2     /2    0     0     0  ' 
      ctir( 7) ='G7   E3/2  2    -/2    0     0     0  '
case("D4h")
      ttir     ='           E    2C4    C2   2C2`  2C2"   I   2IC4   IC2  2IC2` 2IC2"'           
      ctir( 1) ='G1+  A1g   1     1     1     1     1     1     1     1     1     1  '
      ctir( 2) ='G2+  A2g   1     1     1    -1    -1     1     1     1    -1    -1  '
      ctir( 3) ='G3+  B1g   1    -1     1     1    -1     1    -1     1     1    -1  '
      ctir( 4) ='G4+  B2g   1    -1     1    -1     1     1    -1     1    -1     1  '
      ctir( 5) ='G5+  Eg    2     0    -2     0     0     2     0    -2     0     0  '
      ctir( 6) ='G1-  A1u   1     1     1     1     1    -1    -1    -1    -1    -1  '
      ctir( 7) ='G2-  A2u   1     1     1    -1    -1    -1    -1    -1     1     1  '
      ctir( 8) ='G3-  B1u   1    -1     1     1    -1    -1     1    -1    -1     1  '
      ctir( 9) ='G4-  B2u   1    -1     1    -1     1    -1     1    -1     1    -1  '
      ctir(10) ='G5-  Eu    2     0    -2     0     0    -2     0     2     0     0  ' 
      ctir(11) ='G6+  E1/2g 2     /2    0     0     0     2     /2    0     0     0  ' 
      ctir(12) ='G7+  E3/2g 2    -/2    0     0     0     2    -/2    0     0     0  '
      ctir(13) ='G6-  E1/2u 2     /2    0     0     0    -2    -/2    0     0     0  '
      ctir(14) ='G7-  E3/2u 2    -/2    0     0     0    -2     /2    0     0     0  '
case("C3 ") 
      ttir     ='           E     C3    C3-' 
      ctir( 1) ='G1   A     1     1     1  ' 
      ctir( 2) ='G2  2E     1     e     e* ' 
      ctir( 3) ='G3  1E     1     e*    e  ' 
      ctir( 4) ='G4  1E1/2  1    -e*   -e  ' 
      ctir( 5) ='G5  2E1/2  1    -e    -e* ' 
      ctir( 6) ='G6   A3/2  1    -1    -1  ' 
case("C3i")
! #
case("D3 ","C3v")  
      ttir     ='           E    2C3   3C2 ' 
      if(name(1)(1:3).eq."C3v") &
      ttir     ='           E    2C3  3IC2 ' 
      ctir( 1) ='G1   A1    1     1     1  ' 
      ctir( 2) ='G2   A2    1     1    -1  ' 
      ctir( 3) ='G3   E     2    -1     0  ' 
      ctir( 4) ='G4   E1/2  2     1     0  ' 
      ctir( 5) ='G5  1E3/2  1    -1     i  ' 
      ctir( 6) ='G6  2E3/2  1    -1    -i  ' 
case("D3d")
      ttir     ='           E    2C3   3C2    I   2IC3  3IC2 '            
      ctir( 1) ='G1+  A1g   1     1     1     1     1     1  '
      ctir( 2) ='G2+  A2g   1     1    -1     1     1    -1  '
      ctir( 3) ='G3+  Eg    2    -1     0     2    -1     0  '
      ctir( 4) ='G1-  A1u   1     1     1    -1    -1    -1  '
      ctir( 5) ='G2-  A2u   1     1    -1    -1    -1     1  '
      ctir( 6) ='G3-  Eu    2    -1     0    -2     1     0  '
      ctir( 7) ='G4+  E1/2g 2     1     0     2     1     0  '
      ctir( 8) ='G5+ 1E3/2g 1    -1     i     1    -1     i  '
      ctir( 9) ='G6+ 2E3/2g 1    -1    -i     1    -1    -i  '
      ctir(10) ='G4-  E1/2u 2     1     0    -2    -1     0  '
      ctir(11) ='G5- 1E3/2u 1    -1     i    -1     1    -i  '
      ctir(12) ='G6- 2E3/2u 1    -1    -i    -1     1     i  '
case("C6 ","C3h")  
      ttir     ='           E     C6    C3    C2    C3-   C6-' 
      if(name(1)(1:3).eq."C3h") &
      ttir     ='           E    IC6    C3   IC2    C3-  IC6-'
      ctir( 1) ='G1   A     1     1     1     1     1     1  ' 
      ctir( 2) ='G2  2E2    1     e*    e     1     e*    e  ' 
      ctir( 3) ='G3  1E2    1     e     e*    1     e     e* ' 
      ctir( 4) ='G4   B     1    -1     1    -1     1    -1  ' 
      ctir( 5) ='G5  2E1    1    -e*    e    -1     e*   -e  ' 
      ctir( 6) ='G5  1E1    1    -e     e*   -1     e    -e* ' 
      ctir( 7) ='G7  1E1/2  1   -ie    -e*    i    -e    ie* ' 
      ctir( 8) ='G8  2E1/2  1    ie*   -e    -i    -e*  -ie  ' 
      ctir( 9) ='G9  2E5/2  1    ie    -e*   -i    -e   -ie* ' 
      ctir(10) ='G10 1E5/2  1   -ie*   -e     i    -e*   ie  ' 
      ctir(11) ='G11 1E3/2  1    -i    -1     i    -1     i  ' 
      ctir(12) ='G12 2E3/2  1     i    -1    -i    -1    -i  ' 
case("C6h")
! #
case("D6 ","C6v","D3h") 
      ttir     ='           E     C2   2C3   2C6   3C2`  3C2"' 
      if(name(1)(1:3).eq."C6v") &
      ttir     ='           E     C2   2C3   2C6  3IC2` 3IC2"'
      if(name(1)(1:3).eq."D3h") &
      ttir     ='           E    IC2   2C3  2IC6   3C2` 3IC2"'
      ctir( 1) ='G1   A1    1     1     1     1     1     1  ' 
      ctir( 2) ='G2   A2    1     1     1     1    -1    -1  ' 
      ctir( 3) ='G3   B1    1    -1     1    -1     1    -1  ' 
      ctir( 4) ='G4   B2    1    -1     1    -1    -1     1  ' 
      ctir( 5) ='G5   E1    2    -2    -1     1     0     0  ' 
      ctir( 6) ='G6   E2    2     2    -1    -1     0     0  ' 
      ctir( 7) ='G7   E1/2  2     0     1     /3    0     0  ' 
      ctir( 8) ='G8   E5/2  2     0     1    -/3    0     0  ' 
      ctir( 9) ='G9   E3/2  2     0    -2     0     0     0  ' 
case("D6h")
      ttir     ='           E     C2   2C3   2C6   3C2`  3C2"   I    IC2  2IC3  2IC6  3IC2` 3IC2"' 
      ctir( 1) ='G1+  A1g   1     1     1     1     1     1     1     1     1     1     1     1  '
      ctir( 2) ='G2+  A2g   1     1     1     1    -1    -1     1     1     1     1    -1    -1  '
      ctir( 3) ='G3+  B1g   1    -1     1    -1     1    -1     1    -1     1    -1     1    -1  '
      ctir( 4) ='G4+  B2g   1    -1     1    -1    -1     1     1    -1     1    -1    -1     1  '
      ctir( 5) ='G5+  E1g   2    -2    -1     1     0     0     2    -2    -1     1     0     0  '
      ctir( 6) ='G6+  E2g   2     2    -1    -1     0     0     2     2    -1    -1     0     0  '
      ctir( 7) ='G1-  A1u   1     1     1     1     1     1    -1    -1    -1    -1    -1    -1  '
      ctir( 8) ='G2-  A2u   1     1     1     1    -1    -1    -1    -1    -1    -1     1     1  '
      ctir( 9) ='G3-  B1u   1    -1     1    -1     1    -1    -1     1    -1     1    -1     1  '
      ctir(10) ='G4-  B2u   1    -1     1    -1    -1     1    -1     1    -1     1     1    -1  '
      ctir(11) ='G5-  E1u   2    -2    -1     1     0     0    -2     2     1    -1     0     0  '
      ctir(12) ='G6-  E2u   2     2    -1    -1     0     0    -2    -2     1     1     0     0  '
      ctir(13) ='G7+  E1/2g 2     0     1     /3    0     0     2     0     1     /3    0     0  ' 
      ctir(14) ='G8+  E5/2g 2     0     1    -/3    0     0     2     0     1    -/3    0     0  '
      ctir(15) ='G9+  E3/2g 2     0    -2     0     0     0     2     0    -2     0     0     0  '
      ctir(16) ='G7-  E1/2u 2     0     1     /3    0     0    -2     0    -1    -/3    0     0  '
      ctir(17) ='G8-  E5/2u 2     0     1    -/3    0     0    -2     0    -1     /3    0     0  '
      ctir(18) ='G9-  E3/2u 2     0    -2     0     0     0    -2     0     2     0     0     0  '
case("T  ") 
      ttir     ='           E    3C2   4C3   4C3-'  
      ctir( 1) ='G1   A     1     1     1     1  ' 
      ctir( 2) ='G2  1E     1     1     e     e* ' 
      ctir( 3) ='G3  2E     1     1     e*    e  ' 
      ctir( 4) ='G4   T     3    -1     0     0  ' 
      ctir( 5) ='G5   E1/2  2     0     1     1  ' 
      ctir( 6) ='G6  1F3/2  2     0     e     e* ' 
      ctir( 7) ='G7  2F3/2  2     0     e*    e  ' 
case("Th ")
! #
case("O  ","Td ")  
      ttir     ='           E    8C3   3C2   6C4   6C2`'  
      if(name(1)(1:3).eq."Td ") &
      ttir     ='           E    8C3   3C2  6IC4  6IC2 '
      ctir( 1) ='G1   A1    1     1     1     1     1  ' 
      ctir( 2) ='G2   A2    1     1     1    -1    -1  ' 
      ctir( 3) ='G3   E     2    -1     2     0     0  ' 
      ctir( 4) ='G4   T1    3     0    -1     1    -1  ' 
      ctir( 5) ='G5   T2    3     0    -1    -1     1  ' 
      ctir( 6) ='G6   E1/2  2     1     0     /2    0  ' 
      ctir( 7) ='G7   E5/2  2     1     0    -/2    0  ' 
      ctir( 8) ='G8   F3/2  4    -1     0     0     0  ' 
case("Oh ")
      ttir     ='           E    8C3   3C2   6C4   6C2`   I   8IC3  3IC2  6IC4  6IC2`'             
      ctir( 1) ='G1+  A1g   1     1     1     1     1     1     1     1     1     1  '
      ctir( 2) ='G2+  A2g   1     1     1    -1    -1     1     1     1    -1    -1  '
      ctir( 3) ='G3+  Eg    2    -1     2     0     0     2    -1     2     0     0  '
      ctir( 4) ='G4+  T1g   3     0    -1     1    -1     3     0    -1     1    -1  '
      ctir( 5) ='G5+  T2g   3     0    -1    -1     1     3     0    -1    -1     1  '
      ctir( 6) ='G1-  A1u   1     1     1     1     1    -1    -1    -1    -1    -1  '
      ctir( 7) ='G2-  A2u   1     1     1    -1    -1    -1    -1    -1     1     1  '
      ctir( 8) ='G3-  Eu    2    -1     2     0     0    -2     1    -2     0     0  '
      ctir( 9) ='G4-  T1u   3     0    -1     1    -1    -3     0     1    -1     1  '
      ctir(10) ='G5-  T2u   3     0    -1    -1     1    -3     0     1     1    -1  '
      ctir(11) ='G6+  E1/2g 2     1     0     /2    0     2     1     0     /2    0  '
      ctir(12) ='G7+  E5/2g 2     1     0    -/2    0     2     1     0    -/2    0  '
      ctir(13) ='G8+  F3/2g 4    -1     0     0     0     4    -1     0     0     0  '
      ctir(14) ='G6-  E1/2u 2     1     0     /2    0    -2    -1     0    -/2    0  '
      ctir(15) ='G7-  E5/2u 2     1     0    -/2    0    -2    -1     0     /2    0  '
      ctir(16) ='G8-  F3/2u 4    -1     0     0     0    -4     1     0     0     0  '
case default 
  write(*,*) 
  write(*,'("Error(irrepchartab): wrong point group")') 
  write(*,*) 
  write(*,*) "======================  stop ==============" 
end select 
!
! find size of character table ctir( 1)  
ntir(1:2)=0 
do while(ntir(1).lt.48.and.ctir(ntir(1)+1)(1:3).ne.'   ') 
  ntir(1)=ntir(1)+1 
  if(ctir(ntir(1))(8:9).ne.'/2') ntir(2)=ntir(2)+1 
enddo  
!  
!----------------------------------------------------------------------! 
!    transform ctir of the character tables to complex numbers ztir    ! 
!----------------------------------------------------------------------! 
dph=0.0 
do i1=1,ntir(1) 
  j1=0 
  do i2=8+1,8+ntir(2)*6,6 
    j1=j1+1  
    ztir(i1,j1)=cmplx(0.0,0.0) 
    select case( ctir(i1)(i2+3:i2+3) ) 
    case('1') 
      ztir(i1,j1)=1.d0 
    case('2')  
      ztir(i1,j1)=2.d0 
    case('3')  
      ztir(i1,j1)=3.d0 
    case('4') 
      ztir(i1,j1)=4.d0 
    case('i') 
      ztir(i1,j1)=cmplx(0.0,1.0) 
    case('/')  
      read(ctir(i1)(i2+4:i2+4),'(i1)') idi 
      ztir(i1,j1)=sqrt(dble(idi)) 
    case('e')  
      dph=twopi/3.d0 
      if(ctir(i1)(i2+4:i2+4).eq.'*') then 
        ztir(i1,j1)=cmplx(cos(dph),-sin(dph)) 
      else 
        ztir(i1,j1)=cmplx(cos(dph), sin(dph)) 
      endif 
      if(ctir(i1)(i2+2:i2+2).eq.'i')  & 
                  ztir(i1,j1)= ztir(i1,j1)*cmplx(0.0,1.0) 
    case('d') 
      dph=twopi/8.d0 
      if(ctir(i1)(i2+4:i2+4).eq.'*') then 
        ztir(i1,j1)=cmplx(cos(dph),-sin(dph)) 
      else 
        ztir(i1,j1)=cmplx(cos(dph), sin(dph)) 
      endif 
    end select 
    if(ctir(i1)(i2+2:i2+2).eq.'-') ztir(i1,j1)=-ztir(i1,j1) 
    if(ctir(i1)(i2+1:i2+1).eq.'-') ztir(i1,j1)=-ztir(i1,j1) 
  enddo 
enddo 
! 
return 
end subroutine 
!EOC 
