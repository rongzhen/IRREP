      SUBROUTINE DPRODCI(CTMP,B1,N1)
      USE BCLASS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      CHARACTER*80     CTMP(MAXIR),B1  
      CHARACTER*72     C1,C2(MAXIR)   
      DIMENSION        N1(4),N2(4)
!***********************************************************************
!
      I1=N1(4)
      I2=N1(3)-N1(4)
      J1=2*I1
      J2=2*I2
      N2(3)=J1+J2
      N2(4)=J1
      LE1=6*N1(4)
      LE2=6*N2(4)
!
      IE=8
      B1(LE1+1+IE:LE1+6+IE)='   I  '
      B1(LE1+7+IE:LE1+7+LE1+IE) = B1(7+IE:LE1+IE)
      DO 2 I=7,LE1-2
        IF(B1(I+IE:I+IE).EQ.' '.AND.(B1(I+1+IE:I+1+IE).NE.' ')) THEN
            IF( B1(I+1+IE:I+1+IE).EQ.'C' ) THEN
               B1(LE1+I+IE:LE1+I+IE)='I'
            ELSE
               B1(LE1+I+IE:LE1+I+IE)= B1(I+1+IE:I+1+IE)
               B1(LE1+I+1+IE:LE1+I+1+IE)='I'
            ENDIF
        ENDIF 
 2    CONTINUE
!
!.....the double group IRs
      DO 14 I=1,I2
        C1(1:72)='  '//CTMP(I1+I)(11:80)        
        C2(J1+I   )(1:LE2)=C1(1:LE1)//C1(1:LE1)
        C2(J1+I+I2)(1:LE2)=C1(1:LE1)//C1(1:LE1)
        DO 18 J=LE1+3,LE2-3,6
          IF(C2(J1+I+I2)(J+1:J+1).NE.'0') THEN
            IF(C2(J1+I+I2)(J:J).EQ.' ') THEN
              C2(J1+I+I2)(J:J)='-'
            ELSE
              C2(J1+I+I2)(J:J)=' '
            ENDIF
          ENDIF
   18     CONTINUE
        CTMP(J1+I)=CTMP(I1+I)(1:10)//C2(J1+I)(3:72)
        CTMP(J1+I+I2)=CTMP(I1+I)(1:10)//C2(J1+I+I2)(3:72)
!
!.......labeling of symmetric and antisymmetric
        IF(CTMP(J1+I)(3:3).EQ.' ') THEN
          CTMP(J1+I)(3:3)='+'
          CTMP(J1+I+I2)(3:3)='-'
        ELSE
          CTMP(J1+I)(4:4)='+'
          CTMP(J1+I+I2)(4:4)='-'
        ENDIF
        IF(CTMP(J1+I)(6:6).EQ.' ') THEN
          CTMP(J1+I)(6:6)='g'
          CTMP(J1+I+I2)(6:6)='u'
        ELSEIF(CTMP(J1+I)(7:7).EQ.' ') THEN
          CTMP(J1+I)(7:7)='g'
          CTMP(J1+I+I2)(7:7)='u'
        ELSEIF(CTMP(J1+I)(8:8).EQ.' ') THEN
          CTMP(J1+I)(8:8)='g'
          CTMP(J1+I+I2)(8:8)='u'
        ELSEIF(CTMP(J1+I)(9:9).EQ.' ') THEN
          CTMP(J1+I)(9:9)='g'
          CTMP(J1+I+I2)(9:9)='u'
        ELSE
          CTMP(J1+I)(10:10)='g'
          CTMP(J1+I+I2)(10:10)='u'
        ENDIF
 14   CONTINUE
!      
!.....the single group IRs
      DO 4 I=1,I1
        C1(1:72)='  '//CTMP(I)(11:80)
        C2(I   )(1:LE2)=C1(1:LE1)//C1(1:LE1)
        C2(I+I1)(1:LE2)=C1(1:LE1)//C1(1:LE1)
        DO 8 J=LE1+3,LE2-3,6
          IF(C2(I+I1)(J+1:J+1).NE.'0') THEN
            IF(C2(I+I1)(J:J).EQ.' ') THEN
              C2(I+I1)(J:J)='-'
            ELSE
              C2(I+I1)(J:J)=' '
            ENDIF
          ENDIF
    8     CONTINUE 
        CTMP(I)=CTMP(I)(1:10)//C2(I)(3:70)
        CTMP(I+I1)=CTMP(I)(1:10)//C2(I+I1)(3:70)
!
!.......labeling of symmetric and antisymmetric
        IF(CTMP(I)(3:3).EQ.' ') THEN
          CTMP(I)(3:3)='+'
          CTMP(I+I1)(3:3)='-'
        ELSE
          CTMP(I)(4:4)='+'
          CTMP(I+I1)(4:4)='-'
        ENDIF
        IF(CTMP(I)(6:6).EQ.' ') THEN
          CTMP(I)(6:6)='g'
          CTMP(I+I1)(6:6)='u'
        ELSEIF(CTMP(I)(7:7).EQ.' ') THEN
          CTMP(I)(7:7)='g'
          CTMP(I+I1)(7:7)='u'
        ELSEIF(CTMP(I)(8:8).EQ.' ') THEN
          CTMP(I)(8:8)='g'
          CTMP(I+I1)(8:8)='u'
        ELSEIF(CTMP(I)(9:9).EQ.' ') THEN
          CTMP(I)(9:9)='g'
          CTMP(I+I1)(9:9)='u'
        ELSE
          CTMP(I)(10:10)='g'
          CTMP(I+I1)(10:10)='u'
        ENDIF
    4 CONTINUE
!

      DO 121 I=3,4
 121    N1(I)=N2(I)
!
      RETURN
      END


