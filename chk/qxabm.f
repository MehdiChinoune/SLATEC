*DECK QXABM
      SUBROUTINE QXABM (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  QXABM
C***SUBSIDIARY
C***PURPOSE  Test the DEPAC routine DEABM.
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (QXABM-S, QXDABM-D)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  Chow, Jeff, (LANL)
C***DESCRIPTION
C
C *Usage:
C
C        INTEGER  LUN, KPRINT, IPASS
C
C        CALL QXABM (LUN, KPRINT, IPASS)
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
C   DEABM is tested by solving the equations of motion of a body
C   moving in a plane about a spherical earth, namely
C           (D/DT)(D/DT)X = -G*X/R**3
C           (D/DT)(D/DT)Y = -G*Y/R**3
C   where G = 1, R = SQRT(X**2 + Y**2) and
C           X(0) = 1
C           (D/DT)X(0) = 0
C           Y(0) = 0
C           (D/DT)Y(0) = 1.
C
C***ROUTINES CALLED  DEABM, FDEQC, R1MACH
C***REVISION HISTORY  (YYMMDD)
C   810801  DATE WRITTEN
C   890618  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900415  Code extensively revised.  (WRB)
C***END PROLOGUE  QXABM
C
C     Declare arguments.
C
      INTEGER  LUN, KPRINT, IPASS
C
C     Declare local variables.
C
      INTEGER  IDID, INFO(15), IPAR, IWORK(51), N, LIW, LRW
      REAL ABSERR, R, R1MACH, RELERR, RELTOL, RPAR, RWORK(214), T, TOUT,
     +     U(4)
      EXTERNAL FDEQC
C***FIRST EXECUTABLE STATEMENT  QXABM
      IF (KPRINT .GE. 2)  WRITE (LUN, 9000)
C
C     Initialize problem.
C
      N = 4
      LRW = 214
      LIW = 51
      T = 0.0E0
      TOUT = 8.0E0*ATAN(1.0E0)
      U(1) = 1.0E0
      U(2) = 0.0E0
      U(3) = 0.0E0
      U(4) = 1.0E0
      IPASS = 1
      RELTOL = SQRT(R1MACH(4))
      RELERR = 0.1E0*RELTOL
      ABSERR = RELERR**1.5E0
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = 1
      INFO(4) = 0
      IF (KPRINT .GT. 2) WRITE (LUN, 9010) RELERR, ABSERR, T, (1.0E0)
C
  100 CALL DEABM (FDEQC, N, T, U, TOUT, INFO, RELERR, ABSERR,
     +            IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR)
      R = SQRT(U(1)*U(1)+U(2)*U(2))
      IF (ABS(R-1.0E0) .GT. RELTOL) IPASS = 0
      IF (KPRINT .GT. 2) WRITE (LUN, 9020) T, R
      INFO(1) = 1
      IF (IDID .EQ. 1) GO TO 100
C
C     Finish up.
C
      IF (IDID .LT. 1) IPASS = 0
      IF (KPRINT.GT.1 .AND. IDID.LT.1)  WRITE (LUN, 9030) IDID
      IF (KPRINT.GT.1 .AND. IPASS.EQ.1) WRITE (LUN, 9040)
      IF (KPRINT.GE.1 .AND. IPASS.EQ.0) WRITE (LUN, 9050)
      RETURN
C
C     FORMATs.
C
 9000 FORMAT ('1'/' ------------  DEABM QUICK CHECK OUTPUT',
     +        ' ------------')
 9010 FORMAT (/ ' RELERR = ', E16.8, '   ABSERR =', E16.8 /
     +        12X, 'T', 19X, 'R' / 2E20.8)
 9020 FORMAT (2E20.8)
 9030 FORMAT (1X, 'ERROR RETURN FROM DEABM.  IDID = ', I3)
 9040 FORMAT (/ ' ------------  DEABM PASSED TESTS  ------------')
 9050 FORMAT (/ ' ************  DEABM FAILED TESTS  ************')
      END
