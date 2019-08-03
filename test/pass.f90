!** PASS
SUBROUTINE PASS(Lun,Icnt,Itest)
  !> Print a PASS/FAIL message for a particular quick check
  !            test.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920210  PURPOSE added and code restructured.  (WRB)

  INTEGER :: Icnt, Itest, Lun
  !* FIRST EXECUTABLE STATEMENT  PASS
  IF( Itest/=0 ) THEN
    WRITE (Lun,99001) Icnt
    99001 FORMAT (/' TEST NUMBER',I5,' PASSED')
  ELSE
    WRITE (Lun,99002) Icnt
    99002 FORMAT (/' *****TEST NUMBER',I5,' FAILED**********')
  END IF
  !
  RETURN
END SUBROUTINE PASS