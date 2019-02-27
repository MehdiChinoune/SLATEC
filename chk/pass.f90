!DECK PASS
SUBROUTINE PASS(Lun,Icnt,Itest)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  PASS
  !***PURPOSE  Print a PASS/FAIL message for a particular quick check
  !            test.
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920210  PURPOSE added and code restructured.  (WRB)
  !***END PROLOGUE  PASS
  INTEGER Icnt, Itest, Lun
  !***FIRST EXECUTABLE STATEMENT  PASS
  IF ( Itest/=0 ) THEN
    WRITE (Lun,99001) Icnt
    99001 FORMAT (/' TEST NUMBER',I5,' PASSED')
  ELSE
    WRITE (Lun,99002) Icnt
    99002 FORMAT (/' *****TEST NUMBER',I5,' FAILED**********')
  ENDIF
  RETURN
END SUBROUTINE PASS
