!BOP
! !ROUTINE: writetime
! !INTERFACE:
subroutine  writetime(i)
! !INPUT/OUTPUT PARAMETERS:
!   i : input  unit        (in, integer)
!
! !DESCRIPTION:
!   Outputs date and time to unit $i$, using a call to the 
!   built-in routine {\tt date}_{\tt and}_{\tt time}.
!
! !REVISION HISTORY:
!   Created  September 2007 (Clas Persson)
!EOP
!BOC
implicit none
! arguments
integer,    intent(in)  :: i
! local variables
character(10)  cdate,ctime,czone
integer        icvalues(8)
!
call date_and_time(cdate,ctime,czone,icvalues)
write(i,'("### ",a4,"-",a2,"-",a2,2x,a2,":",a2,":",a2,/)' ) &
                          cdate(1:4),cdate(5:6),cdate(7:8), &
                          ctime(1:2),ctime(3:4),ctime(5:6)
return
end subroutine
!EOC

