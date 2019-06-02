!** MPERR
SUBROUTINE MPERR
  !>
  !  Subsidiary to DQDOTA and DQDOTI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (MPERR-A)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !  This routine is called when a fatal error condition is
  !  encountered, and after a message has been written on
  !  logical unit LUN.
  !
  !***
  ! **See also:**  DQDOTA, DQDOTI, MPBLAS
  !***
  ! **Routines called:**  (NONE)
  !***
  ! COMMON BLOCKS    MPCOM

  !* REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   ??????  Modified for use with BLAS.  Blank COMMON changed to named
  !           COMMON.  R given dimension 12.
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
  USE service, ONLY : XERMSG
  !* FIRST EXECUTABLE STATEMENT  MPERR
  CALL XERMSG('MPERR', &
    ' *** EXECUTION TERMINATED BY CALL TO MPERR IN MP VERSION 770217 ***',1,2)
  !
  ! AT PRESENT JUST STOP, BUT COULD DUMP B, T, ETC. HERE.
  ! ACTION COULD EASILY BE CONTROLLED BY A FLAG IN LABELLED COMMON.
  ! ANSI VERSION USES STOP, UNIVAC 1108 VERSION USES
  ! RETURN 0 IN ORDER TO GIVE A TRACE-BACK.
  ! FOR DEBUGGING PURPOSES IT MAY BE USEFUL SIMPLY TO
  ! RETURN HERE.  MOST MP ROUTINES RETURN WITH RESULT
  ! ZERO AFTER CALLING MPERR.
  STOP
END SUBROUTINE MPERR
