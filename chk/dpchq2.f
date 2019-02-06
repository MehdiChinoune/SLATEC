*DECK DPCHQ2
      SUBROUTINE DPCHQ2 (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  DPCHQ2
C***PURPOSE  Test the PCHIP integrators DPCHIA and DPCHID.
C***LIBRARY   SLATEC (PCHIP)
C***TYPE      DOUBLE PRECISION (PCHQK2-S, DPCHQ2-D)
C***KEYWORDS  PCHIP INTEGRATOR QUICK CHECK
C***AUTHOR  Fritsch, F. N., (LLNL)
C***DESCRIPTION
C
C             DPCHIP QUICK CHECK NUMBER 2
C
C     TESTS THE INTEGRATORS:  DPCHIA, DPCHID.
C *Usage:
C
C        INTEGER  LUN, KPRINT, IPASS
C
C        CALL DPCHQ2 (LUN, KPRINT, IPASS)
C
C *Arguments:
C
C     LUN   :IN  is the unit number to which output is to be written.
C
C     KPRINT:IN  controls the amount of output, as specified in the
C                SLATEC Guidelines.
C
C     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
C                IPASS=0 indicates one or more tests failed.
C
C *Description:
C
C   This routine constructs data from a cubic, integrates it with DPCHIA
C   and compares the results with the correct answer.
C   Since DPCHIA calls DPCHID, this tests both integrators.
C
C***ROUTINES CALLED  D1MACH, DPCHIA
C***REVISION HISTORY  (YYMMDD)
C   820601  DATE WRITTEN
C   890306  Changed IPASS to the more accurate name IFAIL.  (FNF)
C   890316  1. Removed IMPLICIT statement.                  (FNF)
C           2. Eliminated unnecessary variable N1.          (FNF)
C           3. Miscellaneous cosmetic changes.              (FNF)
C   891004  Cosmetic changes to prologue.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900314  Improved some output formats.  (FNF)
C   900315  Revised prologue and improved some output formats.  (FNF)
C   900316  Additional minor cosmetic changes.  (FNF)
C   900321  Removed IFAIL from call sequence for SLATEC standards and
C           made miscellaneous cosmetic changes.  (FNF)
C   900323  Corrected list of routines called.  (FNF)
C   901130  Added 1P's to formats; changed to allow KPRINT.gt.3.  (FNF)
C   910708  Minor modifications in use of KPRINT.  (WRB)
C   930317  Improved output formats.  (FNF)
C***END PROLOGUE  DPCHQ2
C
C  Declare arguments.
C
      INTEGER  LUN, KPRINT, IPASS
C
C  DECLARE VARIABLES.
C
      INTEGER  I, IEREXP(17), IERR, IFAIL, N, NPAIRS
      DOUBLE PRECISION
     *      A(17), B(17), CALC, D(7), ERRMAX, ERROR, F(7), MACHEP,
     *      ONE, THREE, THRQTR, TOL, TRUE, TWO, X(7)
      LOGICAL  FAIL, SKIP
C
C  DECLARE EXTERNALS.
C
      DOUBLE PRECISION  DPCHIA, D1MACH
C
C  DEFINE TEST FUNCTIONS.
C
      DOUBLE PRECISION AX, FCN, DERIV, ANTDER
      FCN(AX) = THREE*AX*AX*(AX-TWO)
      DERIV(AX) = THREE*AX*(TWO*(AX-TWO) + AX)
      ANTDER(AX) = AX**3 * (THRQTR*AX - TWO)
C
C  INITIALIZE.
C
      DATA  THRQTR /0.75D0/,  ONE /1.D0/,  TWO /2.D0/,  THREE /3.D0/
      DATA  N /7/
      DATA  X /-4.D0, -2.D0, -0.9D0, 0.D0, 0.9D0, 2.D0, 4.D0/
      DATA  NPAIRS /17/
      DATA  A /-3.0D0, 3.0D0,-0.5D0,-0.5D0,-0.5D0,-4.0D0,-4.0D0, 3.0D0,
     *  -5.0D0,-5.0D0,-6.0D0, 6.0D0,-1.5D0,-1.5D0,-3.0D0, 3.0D0, 0.5D0/
      DATA  B / 3.0D0,-3.0D0, 1.0D0, 2.0D0, 5.0D0,-0.5D0, 4.0D0, 5.0D0,
     *  -3.0D0, 5.0D0,-5.0D0, 5.0D0,-0.5D0,-1.0D0,-2.5D0, 3.5D0, 0.5D0/
      DATA  IEREXP /0,0,0,0,2,0,0,2,1,3,3,3,0,0,0,0,0/
C
C  SET PASS/FAIL TOLERANCE.
C
C***FIRST EXECUTABLE STATEMENT  DPCHQ2
      MACHEP = D1MACH(4)
      TOL = 100.D0*MACHEP
C
C  SET UP PCH FUNCTION DEFINITION.
C
      DO 10  I = 1, N
         F(I) =   FCN(X(I))
         D(I) = DERIV(X(I))
   10 CONTINUE
C
      IF (KPRINT .GE. 3)  WRITE (LUN, 1000)
      IF (KPRINT .GE. 2)  WRITE (LUN, 1001)
      IF (KPRINT .GE. 3)  WRITE (LUN, 1002)  (X(I), F(I), D(I), I=1,N)
C
C  LOOP OVER (A,B)-PAIRS.
C
      IF (KPRINT .GE. 3)  WRITE (LUN, 2000)
C
      IFAIL = 0
C
      SKIP = .FALSE.
      DO 20  I = 1, NPAIRS
C               ---------------------------------------------
         CALC = DPCHIA (N, X, F, D, 1, SKIP, A(I), B(I), IERR)
C               ---------------------------------------------
         IF (IERR .GE. 0)  THEN
            FAIL = IERR .NE. IEREXP(I)
            TRUE = ANTDER(B(I)) - ANTDER(A(I))
            ERROR = CALC - TRUE
            IF (KPRINT .GE. 3)  THEN
               IF (FAIL)  THEN
                 WRITE (LUN, 2001) A(I), B(I), IERR, TRUE, CALC, ERROR,
     *                                          IEREXP(I)
               ELSE
                 WRITE (LUN, 2002) A(I), B(I), IERR, TRUE, CALC, ERROR
               ENDIF
            ENDIF
C
            ERROR = ABS(ERROR) / MAX(ONE, ABS(TRUE))
            IF (FAIL .OR. (ERROR.GT.TOL))  IFAIL = IFAIL + 1
            IF (I .EQ. 1)  THEN
               ERRMAX = ERROR
            ELSE
               ERRMAX = MAX(ERRMAX, ERROR)
            ENDIF
         ELSE
            IF (KPRINT .GE. 3)  WRITE (LUN, 2002)  A(I), B(I), IERR
            IFAIL = IFAIL + 1
         ENDIF
   20 CONTINUE
C
C  PRINT SUMMARY.
C
      IF (KPRINT .GE. 2)  THEN
         WRITE (LUN, 2003)  ERRMAX, TOL
         IF (IFAIL .NE. 0)  WRITE (LUN, 3001)  IFAIL
      ENDIF
C
C  TERMINATE.
C
      IF (IFAIL.EQ.0)  THEN
         IPASS = 1
         IF (KPRINT.GE.2) WRITE(LUN,99998)
      ELSE
         IPASS = 0
         IF (KPRINT.GE.1) WRITE(LUN,99999)
      ENDIF
C
      RETURN
C
C  FORMATS.
C
 1000 FORMAT ('1'//10X,'TEST DPCHIP INTEGRATORS')
 1001 FORMAT (//10X,'DPCHQ2 RESULTS'/10X,'--------------')
 1002 FORMAT (// 5X,'DATA:' //11X,'X',9X,'F',9X,'D' /(5X,3F10.3) )
 2000 FORMAT (// 5X,'TEST RESULTS:'
     *        //'    A     B    ERR     TRUE',16X,'CALC',15X,'ERROR')
 2001 FORMAT (2F6.1,I5,1P,2D20.10,D15.5,'  (',I1,') *****' )
 2002 FORMAT (2F6.1,I5,1P,2D20.10,D15.5)
 2003 FORMAT (/'  MAXIMUM RELATIVE ERROR IS:',1P,D15.5,
     *                       ',   TOLERANCE:',1P,D15.5)
 3001 FORMAT (/' *** TROUBLE ***',I5,' INTEGRATION TESTS FAILED.')
99998 FORMAT (/' ------------ DPCHIP PASSED  ALL INTEGRATION TESTS',
     *        ' ------------')
99999 FORMAT (/' ************ DPCHIP FAILED SOME INTEGRATION TESTS',
     *        ' ************')
C------------- LAST LINE OF DPCHQ2 FOLLOWS -----------------------------
      END
