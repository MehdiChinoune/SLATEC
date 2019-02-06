*DECK EVERCK
      SUBROUTINE EVERCK (LOUT, KPRINT, FAIL)
C***BEGIN PROLOGUE  EVERCK
C***SUBSIDIARY
C***PURPOSE  Test error returns from PCHIP evaluators for PCHQK1.
C***LIBRARY   SLATEC (PCHIP)
C***TYPE      SINGLE PRECISION (EVERCK-S, DEVERK-D)
C***KEYWORDS  PCHIP EVALUATOR QUICK CHECK
C***AUTHOR  Fritsch, F. N., (LLNL)
C***DESCRIPTION
C
C --------- CODE TO TEST ERROR RETURNS FROM PCHIP EVALUATORS. ---------
C
C
C     FORTRAN LIBRARY ROUTINES USED:  (WRITE).
C     SLATEC LIBRARY ROUTINES USED:  CHFDV, CHFEV, PCHFD, PCHFE,
C                                    XERDMP, XGETF, XSETF.
C     OTHER ROUTINES USED:  COMP.
C
C***ROUTINES CALLED  CHFDV, CHFEV, COMP, PCHFD, PCHFE, XERDMP, XGETF,
C                    XSETF
C***REVISION HISTORY  (YYMMDD)
C   820601  DATE WRITTEN
C   820715  CONVERTED TO QUICK CHECK FOR SLATEC LIBRARY.
C   890207  ADDED CALLS TO ERROR HANDLER.
C   890316  Added call to XERDMP if KPRINT.GT.2 (FNF).
C   890629  Appended E0 to real constants to reduce S.P./D.P.
C           differences.
C   890706  Cosmetic changes to prologue.  (WRB)
C   891009  Removed unreferenced statement label.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900309  Added COMP to list of routines called.  (FNF)
C   900315  Revised prologue and improved some output formats.  (FNF)
C   900316  Deleted INCFD tests because some compilers object to them,
C           and made additional minor cosmetic changes.  (FNF)
C   900322  Made miscellaneous cosmetic changes.  (FNF)
C   910708  Minor modifications in use of KPRINT.  (WRB)
C   930504  Removed parens from constants in WRITE statements.  (FNF)
C***END PROLOGUE  EVERCK
C
C  Declare arguments.
C
      INTEGER  LOUT, KPRINT
      LOGICAL  FAIL
C
C  DECLARATIONS.
C
      INTEGER  I, IERR, KONTRL, N, NERR, NEXT(2)
      REAL  D(10), DUM, F(10), TEMP, X(10)
      LOGICAL  COMP, SKIP
C
C  INITIALIZE.
C
      PARAMETER (N = 10)
C***FIRST EXECUTABLE STATEMENT  EVERCK
      NERR = 0
C
      CALL XGETF (KONTRL)
      IF (KPRINT .LE. 2) THEN
         CALL XSETF (0)
      ELSE
         CALL XSETF (1)
      ENDIF
C
      IF (KPRINT .GE. 3)  WRITE (LOUT, 2000)
      IF (KPRINT .GE. 2)  WRITE (LOUT, 5000)
C
C  FIRST, TEST CHFEV AND CHFDV.
C
      IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -1
      CALL CHFEV (0.E0, 1.E0, 3.E0, 7.E0, 3.E0, 6.E0, 0, DUM, DUM,
     * NEXT, IERR)
      IF (.NOT. COMP (IERR, -1, LOUT, KPRINT) )  NERR = NERR + 1
C
      IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -2
      CALL CHFEV (1.E0, 1.E0, 3.E0, 7.E0, 3.E0, 6.E0, 1, DUM, DUM,
     * NEXT, IERR)
      IF (.NOT. COMP (IERR, -2, LOUT, KPRINT) )  NERR = NERR + 1
C
      IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -1
      CALL CHFDV (0.E0, 1.E0, 3.E0, 7.E0, 3.E0, 6.E0, 0, DUM, DUM, DUM,
     * NEXT, IERR)
      IF (.NOT. COMP (IERR, -1, LOUT, KPRINT) )  NERR = NERR + 1
C
      IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -2
      CALL CHFDV (1.E0, 1.E0, 3.E0, 7.E0, 3.E0, 6.E0, 1, DUM, DUM, DUM,
     * NEXT, IERR)
      IF (.NOT. COMP (IERR, -2, LOUT, KPRINT) )  NERR = NERR + 1
C
C  SET UP PCH DEFINITION.
C
      DO 10  I = 1, N
         X(I) = I
         F(I) = I + 2
         D(I) = 1.E0
   10 CONTINUE
C
C  SWAP POINTS 4 AND 7, SO X-ARRAY IS OUT OF ORDER.
C
      TEMP = X(4)
      X(4) = X(7)
      X(7) = TEMP
C
C  NOW, TEST PCHFE AND PCHFD.
C
      IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -1
      SKIP = .FALSE.
      CALL PCHFE (1, X, F, D, 1, SKIP, 0, DUM, DUM, IERR)
      IF (.NOT. COMP (IERR, -1, LOUT, KPRINT) )  NERR = NERR + 1
C
      IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -3
      SKIP = .FALSE.
      CALL PCHFE (N, X, F, D, 1, SKIP, 0, DUM, DUM, IERR)
      IF (.NOT. COMP (IERR, -3, LOUT, KPRINT) )  NERR = NERR + 1
C
      IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -4
      SKIP = .TRUE.
      CALL PCHFE (N, X, F, D, 1, SKIP, 0, DUM, DUM, IERR)
      IF (.NOT. COMP (IERR, -4, LOUT, KPRINT) )  NERR = NERR + 1
C
      IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -1
      SKIP = .FALSE.
      CALL PCHFD (1, X, F, D, 1, SKIP, 0, DUM, DUM, DUM, IERR)
      IF (.NOT. COMP (IERR, -1, LOUT, KPRINT) )  NERR = NERR + 1
C
      IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -3
      SKIP = .FALSE.
      CALL PCHFD (N, X, F, D, 1, SKIP, 0, DUM, DUM, DUM, IERR)
      IF (.NOT. COMP (IERR, -3, LOUT, KPRINT) )  NERR = NERR + 1
C
      IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -4
      SKIP = .TRUE.
      CALL PCHFD (N, X, F, D, 1, SKIP, 0, DUM, DUM, DUM, IERR)
      IF (.NOT. COMP (IERR, -4, LOUT, KPRINT) )  NERR = NERR + 1
C
C  SUMMARIZE RESULTS.
C
      IF (KPRINT .GT. 2)  CALL XERDMP
      IF (NERR .EQ. 0)  THEN
         FAIL = .FALSE.
         IF (KPRINT .GE. 2)  WRITE (LOUT, 5002)
      ELSE
         FAIL = .TRUE.
         IF (KPRINT .GE. 2)  WRITE (LOUT, 5003)  NERR
      ENDIF
C
C  TERMINATE.
C
      CALL XSETF (KONTRL)
      RETURN
C
C  FORMATS.
C
 2000 FORMAT ('1'//10X,'TEST ERROR RETURNS')
 5000 FORMAT (//10X,'EVERCK RESULTS'/10X,'--------------')
 5001 FORMAT (/' THIS CALL SHOULD RETURN IERR =',I3)
 5002 FORMAT (/' ALL ERROR RETURNS OK.')
 5003 FORMAT (//' ***** TROUBLE IN EVERCK *****'
     *        //5X,I5,' TESTS FAILED TO GIVE EXPECTED RESULTS.')
C------------- LAST LINE OF EVERCK FOLLOWS -----------------------------
      END
