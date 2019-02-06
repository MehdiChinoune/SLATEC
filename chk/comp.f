*DECK COMP
      LOGICAL FUNCTION COMP (IERACT, IEREXP, LOUT, KPRINT)
C***BEGIN PROLOGUE  COMP
C***SUBSIDIARY
C***PURPOSE  Compare actual and expected values of error flag.
C***LIBRARY   SLATEC
C***KEYWORDS  QUICK CHECK SERVICE ROUTINE
C***AUTHOR  Fritsch, F. N., (LLNL)
C***DESCRIPTION
C
C     COMPARE ACTUAL VALUE OF IERR WITH EXPECTED VALUE.
C        PRINT ERROR MESSAGE IF THEY DON'T AGREE.
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   820601  DATE WRITTEN
C   890618  REVISION DATE from Version 3.2
C   890706  Cosmetic changes to prologue.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  Revised prologue.  (FNF)
C   900316  Minor modification to format 5010.  (FNF)
C   910708  Minor modifications in use of KPRINT.  (WRB)
C***END PROLOGUE  COMP
      INTEGER  IERACT, IEREXP, LOUT, KPRINT
C***FIRST EXECUTABLE STATEMENT  COMP
      IF (IERACT .EQ. IEREXP)  THEN
         COMP = .TRUE.
         IF (KPRINT .GE. 3)  WRITE (LOUT, 5010)
 5010    FORMAT ('     OK.')
      ELSE
         COMP = .FALSE.
         IF (KPRINT .GE. 3)  WRITE (LOUT, 5020)  IERACT
 5020    FORMAT (' *** COMPARE FAILED -- IERR =',I5)
      ENDIF
C
      RETURN
C------------- LAST LINE OF COMP FOLLOWS -----------------------------
      END
