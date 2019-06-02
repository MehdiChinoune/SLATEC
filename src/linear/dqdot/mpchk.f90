!** MPCHK
SUBROUTINE MPCHK(I,J)
  !>
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
  USE MPCOM, ONLY : b_com, lun_com, m_com, t_com, mxr_com
  USE service, ONLY : I1MACH
  INTEGER :: I, J
  INTEGER :: ib, mx
  !* FIRST EXECUTABLE STATEMENT  MPCHK
  lun_com = I1MACH(4)
  ! NOW CHECK LEGALITY OF B, T AND M
  IF ( b_com<=1 ) THEN
    WRITE (lun_com,99001) b_com
    99001 FORMAT (' *** b_com =',I10,' ILLEGAL IN CALL TO MPCHK,'/&
      ' PERHAPS NOT SET BEFORE CALL TO AN MP ROUTINE ***')
    CALL MPERR
  END IF
  IF ( t_com<=1 ) THEN
    WRITE (lun_com,99002) t_com
    99002 FORMAT (' *** t_com =',I10,' ILLEGAL IN CALL TO MPCHK,'/&
      ' PERHAPS NOT SET BEFORE CALL TO AN MP ROUTINE ***')
    CALL MPERR
  END IF
  IF ( m_com<=t_com ) THEN
    WRITE (lun_com,99003)
    99003 FORMAT (' *** m_com .LE. t_com IN CALL TO MPCHK,'/&
      ' PERHAPS NOT SET BEFORE CALL TO AN MP ROUTINE ***')
    CALL MPERR
  END IF
  ! 8*B*B-1 SHOULD BE REPRESENTABLE, IF NOT WILL OVERFLOW
  ! AND MAY BECOME NEGATIVE, SO CHECK FOR THIS
  ib = 4*b_com*b_com - 1
  IF ( (ib<=0).OR.((2*ib+1)<=0) ) THEN
    WRITE (lun_com,99004)
    99004 FORMAT (' *** b_com TOO LARGE IN CALL TO MPCHK ***')
    CALL MPERR
  END IF
  ! CHECK THAT SPACE IN COMMON IS SUFFICIENT
  mx = I*t_com + J
  IF ( mxr_com>=mx ) RETURN
  ! HERE COMMON IS TOO SMALL, SO GIVE ERROR MESSAGE.
  WRITE (lun_com,99005) I, J, mx, mxr_com, t_com
  99005 FORMAT (' *** mxr_com TOO SMALL OR NOT SET TO DIM(r_com) BEFORE CALL',&
    ' TO AN MP ROUTINE *** '/' *** mxr_com SHOULD BE AT LEAST',I3,'*t_com +',&
    I4,' =',I6,'  ***'/' *** ACTUALLY mxr_com =',I10,', AND t_com =',I10,'  ***')
  CALL MPERR
END SUBROUTINE MPCHK
