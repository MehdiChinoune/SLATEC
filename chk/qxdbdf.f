*DECK QXDBDF
      SUBROUTINE QXDBDF (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  QXDBDF
C***PURPOSE  Test the DEPAC routine DDEBDF.
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (QXBDF-S, QXDBDF-D)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  Chow, Jeff, (LANL)
C***DESCRIPTION
C
C *Usage:
C
C        INTEGER  LUN, KPRINT, IPASS
C
C        CALL QXDBDF (LUN, KPRINT, IPASS)
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
C   DDEBDF is tested by solving the equations of motion of a body
C   moving in a plane about a spherical earth, namely
C           (D/DT)(D/DT)X = -G*X/R**3
C           (D/DT)(D/DT)Y = -G*Y/R**3
C   where G = 1, R = SQRT(X**2 + Y**2) and
C           X(0) = 1
C           (D/DT)X(0) = 0
C           Y(0) = 0
C           (D/DT)Y(0) = 1.
C
C***ROUTINES CALLED  D1MACH, DDEBDF, DFDEQC, DJAC
C***REVISION HISTORY  (YYMMDD)
C   810801  DATE WRITTEN
C   890618  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900415  Code extensively revised.  (WRB)
C***END PROLOGUE  QXDBDF
C
C     Declare arguments.
C
      INTEGER  LUN, KPRINT, IPASS
C
C     Declare local variables.
C
      INTEGER  IDID, INFO(15), IPAR, IWORK(60), N, LIW, LRW, NSTEP
      DOUBLE PRECISION ABSERR, D1MACH, R, RELTOL, RELERR, RPAR,
     +                 RWORK(306), T, TOUT, U(4)
      EXTERNAL DFDEQC, DJAC
C***FIRST EXECUTABLE STATEMENT  QXDBDF
      IF (KPRINT .GE. 2)  WRITE (LUN, 9000)
C
C     Initialize problem.
C
      N = 4
      LRW = 306
      LIW = 60
      T = 0.0D0
      TOUT = 8.0D0*ATAN(1.0D0)
      U(1) = 1.0D0
      U(2) = 0.0D0
      U(3) = 0.0D0
      U(4) = 1.0D0
      IPASS = 1
      NSTEP = 0
      RELTOL = MAX(SQRT(D1MACH(4)),1.D-9)
      RELERR = MAX(0.0001D0*RELTOL,1.D-12)
      ABSERR = RELERR**1.5D0
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = 1
      INFO(4) = 0
      INFO(5) = 1
      INFO(6) = 0
      IF (KPRINT .GT. 2) WRITE (LUN, 9010) RELERR, ABSERR, T, (1.0D0)
C
  100 CALL DDEBDF (DFDEQC, N, T, U, TOUT, INFO, RELERR, ABSERR,
     +             IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, DJAC)
      R = SQRT(U(1)*U(1)+U(2)*U(2))
      IF (ABS(R-1.0D0) .GT. RELTOL) IPASS = 0
      IF (KPRINT .GT. 2) WRITE (LUN, 9020) T, R
      INFO(1) = 1
      IF (IDID .EQ. 1) GO TO 100
C
C     For the double precision version, we allow the integrator to take
C     up to 2000 steps before we declare failure.
C
      IF (IDID .EQ. -1) THEN
         NSTEP = NSTEP + 500
         IF (NSTEP .LT. 2000) GOTO 100
      ENDIF
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
 9000 FORMAT ('1'/' ------------  DDEBDF QUICK CHECK OUTPUT',
     +        ' ------------')
 9010 FORMAT (/ ' RELERR = ', D16.8, '   ABSERR =', D16.8 /
     +        12X, 'T', 19X, 'R' / 2D20.8)
 9020 FORMAT (2D20.8)
 9030 FORMAT (1X, 'ERROR RETURN FROM DDEBDF.  IDID = ', I3)
 9040 FORMAT (/ ' ------------  DDEBDF PASSED TESTS  ------------')
 9050 FORMAT (/ ' ************  DDEBDF FAILED TESTS  ************')
      END
