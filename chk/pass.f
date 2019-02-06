*DECK PASS
      SUBROUTINE PASS (LUN, ICNT, ITEST)
C***BEGIN PROLOGUE  PASS
C***PURPOSE  Print a PASS/FAIL message for a particular quick check
C            test.
C***LIBRARY   SLATEC
C***AUTHOR  (UNKNOWN)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920210  PURPOSE added and code restructured.  (WRB)
C***END PROLOGUE  PASS
      INTEGER ICNT, ITEST, LUN
C***FIRST EXECUTABLE STATEMENT  PASS
      IF (ITEST .NE. 0) THEN
        WRITE (LUN,9000) ICNT
      ELSE
        WRITE (LUN,9100) ICNT
      ENDIF
      RETURN
 9000 FORMAT(/ ' TEST NUMBER', I5, ' PASSED')
 9100 FORMAT(/ ' *****TEST NUMBER' ,I5, ' FAILED**********')
      END
