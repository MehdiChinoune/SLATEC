!** MPCHK
SUBROUTINE MPCHK(I,J)
  USE MPCOM
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to DQDOTA and DQDOTI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (MPCHK-A)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !  Checks legality of B, T, M, MXR and LUN which should be set
  !  in COMMON. The condition on MXR (the dimension of the EP arrays)
  !  is that  MXR .GE. (I*T + J)
  !
  !***
  ! **See also:**  DQDOTA, DQDOTI, MPBLAS
  !***
  ! **Routines called:**  I1MACH, MPERR
  !***
  ! COMMON BLOCKS    MPCOM

  !* REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   ??????  Modified for use with BLAS.  Blank COMMON changed to named
  !           COMMON.  R given dimension 12.
  !   891009  Removed unreferenced statement label.  (WRB)
  !   891009  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   930124  Increased Array size in MPCON for SUN -r8.  (RWC)

  INTEGER I, I1MACH, ib, J, mx
  !* FIRST EXECUTABLE STATEMENT  MPCHK
  LUN = I1MACH(4)
  ! NOW CHECK LEGALITY OF B, T AND M
  IF ( B<=1 ) THEN
    WRITE (LUN,99001) B
    99001 FORMAT (' *** B =',I10,' ILLEGAL IN CALL TO MPCHK,'/&
      ' PERHAPS NOT SET BEFORE CALL TO AN MP ROUTINE ***')
    CALL MPERR
  END IF
  IF ( T<=1 ) THEN
    WRITE (LUN,99002) T
    99002 FORMAT (' *** T =',I10,' ILLEGAL IN CALL TO MPCHK,'/&
      ' PERHAPS NOT SET BEFORE CALL TO AN MP ROUTINE ***')
    CALL MPERR
  END IF
  IF ( M<=T ) THEN
    WRITE (LUN,99003)
    99003 FORMAT (' *** M .LE. T IN CALL TO MPCHK,'/&
      ' PERHAPS NOT SET BEFORE CALL TO AN MP ROUTINE ***')
    CALL MPERR
  END IF
  ! 8*B*B-1 SHOULD BE REPRESENTABLE, IF NOT WILL OVERFLOW
  ! AND MAY BECOME NEGATIVE, SO CHECK FOR THIS
  ib = 4*B*B - 1
  IF ( (ib<=0).OR.((2*ib+1)<=0) ) THEN
    WRITE (LUN,99004)
    99004 FORMAT (' *** B TOO LARGE IN CALL TO MPCHK ***')
    CALL MPERR
  END IF
  ! CHECK THAT SPACE IN COMMON IS SUFFICIENT
  mx = I*T + J
  IF ( MXR>=mx ) RETURN
  ! HERE COMMON IS TOO SMALL, SO GIVE ERROR MESSAGE.
  WRITE (LUN,99005) I, J, mx, MXR, T
  99005 FORMAT (' *** MXR TOO SMALL OR NOT SET TO DIM(R) BEFORE CALL',&
    ' TO AN MP ROUTINE *** '/' *** MXR SHOULD BE AT LEAST',I3,'*T +',&
    I4,' =',I6,'  ***'/' *** ACTUALLY MXR =',I10,', AND T =',I10,'  ***')
  CALL MPERR
END SUBROUTINE MPCHK