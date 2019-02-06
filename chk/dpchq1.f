*DECK DPCHQ1
      SUBROUTINE DPCHQ1 (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  DPCHQ1
C***PURPOSE  Test the PCHIP evaluators DCHFDV, DCHFEV, DPCHFD, DPCHFE.
C***LIBRARY   SLATEC (PCHIP)
C***TYPE      DOUBLE PRECISION (PCHQK1-S, DPCHQ1-D)
C***KEYWORDS  PCHIP EVALUATOR QUICK CHECK
C***AUTHOR  Fritsch, F. N., (LLNL)
C***DESCRIPTION
C
C             DPCHIP QUICK CHECK NUMBER 1
C
C     TESTS THE EVALUATORS:  DCHFDV, DCHFEV, DPCHFD, DPCHFE.
C *Usage:
C
C        INTEGER  LUN, KPRINT, IPASS
C
C        CALL DPCHQ1 (LUN, KPRINT, IPASS)
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
C   This routine carries out three tests of the PCH evaluators:
C     DEVCHK tests the single-cubic evaluators.
C     DEVPCK tests the full PCH evaluators.
C     DEVERK exercises the error returns in all evaluators.
C
C***ROUTINES CALLED  DEVCHK, DEVERK, DEVPCK
C***REVISION HISTORY  (YYMMDD)
C   820601  DATE WRITTEN
C   890306  Changed IPASS to the more accurate name IFAIL.  (FNF)
C   890307  Removed conditional on call to DEVERK.
C   890706  Cosmetic changes to prologue.  (WRB)
C   891004  Correction in prologue.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900309  Added DEVERK to list of routines called.  (FNF)
C   900314  Improved some output formats.
C   900315  Revised prologue and improved some output formats.  (FNF)
C   900316  Additional minor cosmetic changes.  (FNF)
C   900321  Removed IFAIL from call sequence for SLATEC standards and
C           made miscellaneous cosmetic changes.  (FNF)
C   930317  Improved output formats.  (FNF)
C***END PROLOGUE  DPCHQ1
C
C  Declare arguments.
C
      INTEGER  LUN, KPRINT, IPASS
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  I1, I2, I3, I4, I5, I6, I7, I8, I9, IFAIL, NPTS
      DOUBLE PRECISION  WORK (4000)
      LOGICAL  FAIL
C
C***FIRST EXECUTABLE STATEMENT  DPCHQ1
      IF (KPRINT .GE. 2)  WRITE (LUN, 1000) KPRINT
C
C  TEST DCHFDV AND DCHFEV.
C
      IFAIL = 0
      NPTS = 1000
      I1 = 1  + NPTS
      I2 = I1 + NPTS
      I3 = I2 + NPTS
      CALL DEVCHK (LUN, KPRINT, NPTS, WORK(1), WORK(I1), WORK(I2),
     *                                          WORK(I3), FAIL)
      IF (FAIL)  IFAIL = IFAIL + 1
C
C  TEST DPCHFD AND DPCHFE.
C
      I1 = 1  +  10
      I2 = I1 +  10
      I3 = I2 + 100
      I4 = I3 + 100
      I5 = I4 + 100
      I6 = I5 +  51
      I7 = I6 +  51
      I8 = I7 +  51
      I9 = I8 +  51
      CALL DEVPCK (LUN, KPRINT, WORK(1), WORK(I1), WORK(I2), WORK(I3),
     *             WORK(I4), WORK(I5), WORK(I6), WORK(I7), WORK(I8),
     *             WORK(I9), FAIL)
      IF (FAIL)  IFAIL = IFAIL + 2
C
C  TEST ERROR RETURNS.
C
      CALL DEVERK (LUN, KPRINT, FAIL)
      IF (FAIL)  IFAIL = IFAIL + 4
C
C  PRINT SUMMARY AND TERMINATE.
C     At this point, IFAIL has the following value:
C        IFAIL = 0  IF ALL TESTS PASSED.
C        IFAIL BETWEEN 1 AND 7 IS THE SUM OF:
C           IFAIL=1  IF SINGLE CUBIC  TEST FAILED. (SEE DEVCHK OUTPUT.)
C           IFAIL=2  IF DPCHFD/DPCHFE TEST FAILED. (SEE DEVPCK OUTPUT.)
C           IFAIL=4  IF ERROR RETURN  TEST FAILED. (SEE DEVERK OUTPUT.)
C
      IF ((KPRINT.GE.2).AND.(IFAIL.NE.0))  WRITE (LUN, 3001)  IFAIL
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
 1000 FORMAT ('1'/' ------------ DPCHIP QUICK CHECK OUTPUT',
     *        ' ------------' //20X,'( KPRINT =',I2,' )')
 3001 FORMAT (/' *** TROUBLE ***',I5,' EVALUATION TESTS FAILED.')
99998 FORMAT (/' ------------ DPCHIP PASSED  ALL EVALUATION TESTS',
     *        ' ------------')
99999 FORMAT (/' ************ DPCHIP FAILED SOME EVALUATION TESTS',
     *        ' ************')
C------------- LAST LINE OF DPCHQ1 FOLLOWS -----------------------------
      END
