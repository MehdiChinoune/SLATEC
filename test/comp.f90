!** COMP
LOGICAL FUNCTION COMP(Ieract,Ierexp,Lout,Kprint)
  !>
  !***
  !  Compare actual and expected values of error flag.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Keywords:**  QUICK CHECK SERVICE ROUTINE
  !***
  ! **Author:**  Fritsch, F. N., (LLNL)
  !***
  ! **Description:**
  !
  !     COMPARE ACTUAL VALUE OF IERR WITH EXPECTED VALUE.
  !        PRINT ERROR MESSAGE IF THEY DON'T AGREE.
  !
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   820601  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   890706  Cosmetic changes to prologue.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  Revised prologue.  (FNF)
  !   900316  Minor modification to format 5010.  (FNF)
  !   910708  Minor modifications in use of KPRINT.  (WRB)

  INTEGER Ieract, Ierexp, Lout, Kprint
  !* FIRST EXECUTABLE STATEMENT  COMP
  IF ( Ieract==Ierexp ) THEN
    COMP = .TRUE.
    IF ( Kprint>=3 ) WRITE (Lout,99001)
    99001 FORMAT ('     OK.')
  ELSE
    COMP = .FALSE.
    IF ( Kprint>=3 ) WRITE (Lout,99002) Ieract
    99002 FORMAT (' *** COMPARE FAILED -- IERR =',I5)
  END IF
  !
  !------------- LAST LINE OF COMP FOLLOWS -----------------------------
END FUNCTION COMP
