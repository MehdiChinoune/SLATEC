!** MPUNFL
SUBROUTINE MPUNFL(X)
  !>
  !  Subsidiary to DQDOTA and DQDOTI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (MPUNFL-A)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  ! Called on multiple-precision underflow, i.e.  when the
  ! exponent of 'mp' number X would be less than -M.
  !
  !***
  ! **See also:**  DQDOTA, DQDOTI
  !***
  ! **Routines called:**  MPCHK

  !* REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  
  INTEGER X(*)
  !* FIRST EXECUTABLE STATEMENT  MPUNFL
  CALL MPCHK(1,4)
  ! THE UNDERFLOWING NUMBER IS SET TO ZERO
  ! AN ALTERNATIVE WOULD BE TO CALL MPMINR (X) AND RETURN,
  ! POSSIBLY UPDATING A COUNTER AND TERMINATING EXECUTION
  ! AFTER A PRESET NUMBER OF UNDERFLOWS.  ACTION COULD EASILY
  ! BE DETERMINED BY A FLAG IN LABELLED COMMON.
  X(1) = 0
END SUBROUTINE MPUNFL
