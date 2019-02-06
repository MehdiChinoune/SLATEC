*DECK DSOSQX
      SUBROUTINE DSOSQX (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  DSOSQX
C***PURPOSE  Quick check for DSOS.
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (SOSNQX-S, DSOSQX-D)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C   This subroutine performs a quick check on the subroutine DSOS.
C
C***ROUTINES CALLED  D1MACH, DNRM2, DSOS, DSOSFN, PASS
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890618  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Code cleaned up and TYPE section added.  (RWC, WRB)
C***END PROLOGUE  DSOSQX
C     .. Scalar Arguments ..
      INTEGER IPASS, KPRINT, LUN
C     .. Local Scalars ..
      DOUBLE PRECISION AER, FNORM, FNORMS, RER, TOLF
      INTEGER ICNT, IFLAG, IFLAGS, LIW, LWA, N
C     .. Local Arrays ..
      DOUBLE PRECISION FVEC(2), WA(17), X(2)
      INTEGER ITEST(2), IW(6)
C     .. External Functions ..
      DOUBLE PRECISION D1MACH, DNRM2, DSOSFN
      EXTERNAL D1MACH, DNRM2, DSOSFN
C     .. External Subroutines ..
      EXTERNAL DSOS, PASS
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C***FIRST EXECUTABLE STATEMENT  DSOSQX
      IFLAGS = 3
      FNORMS = 0.0D0
      N = 2
      LWA = 17
      LIW = 6
      TOLF = SQRT(D1MACH(4))
      RER = SQRT(D1MACH(4))
      AER = 0.0D0
      IF (KPRINT .GE. 2) WRITE (LUN,9000)
C
C     Test the code with proper input values.
C
      IFLAG = 0
      X(1) = -1.2D0
      X(2) = 1.0D0
      CALL DSOS (DSOSFN,N,X,RER,AER,TOLF,IFLAG,WA,LWA,IW,LIW)
      ICNT = 1
      FVEC(1) = DSOSFN(X,1)
      FVEC(2) = DSOSFN(X,2)
      FNORM = DNRM2(N,FVEC,1)
      ITEST(ICNT) = 0
      IF (IFLAG.LE.IFLAGS .AND. FNORM-FNORMS.LE.RER) ITEST(ICNT) = 1
C
      IF (KPRINT .NE. 0) THEN
         IF (KPRINT.GE.3 .OR. (KPRINT.GE.2.AND.ITEST(ICNT).NE.1))
     +       WRITE (LUN,9010) IFLAGS,FNORMS,IFLAG,FNORM
         IF (KPRINT.GE.2 .OR. (KPRINT.EQ.1.AND.ITEST(ICNT).NE.1))
     +       CALL PASS (LUN,ICNT,ITEST(ICNT))
      ENDIF
C
C     Test improper input parameters.
C
      LWA = 15
      IFLAG = 0
      X(1) = -1.2D0
      X(2) = 1.0D0
      CALL DSOS (DSOSFN,N,X,RER,AER,TOLF,IFLAG,WA,LWA,IW,LIW)
      ICNT = 2
      ITEST(ICNT) = 0
      IF (IFLAG .EQ. 9) ITEST(ICNT) = 1
      IF (KPRINT.GE.2 .OR. (KPRINT.EQ.1.AND.ITEST(ICNT).NE.1))
     +    CALL PASS (LUN,ICNT,ITEST(ICNT))
C
C     Set IPASS.
C
      IPASS = ITEST(1)*ITEST(2)
      IF (KPRINT.GE.1 .AND. IPASS.NE.1) WRITE (LUN,9020)
      IF (KPRINT.GE.2 .AND. IPASS.EQ.1) WRITE (LUN,9030)
      RETURN
 9000 FORMAT ('1' / '  DSOS QUICK CHECK' /)
 9010 FORMAT (' EXPECTED VALUE OF IFLAG AND RESIDUAL NORM', I5, D20.5 /
     +        ' RETURNED VALUE OF IFLAG AND RESIDUAL NORM', I5, D20.5 /)
 9020 FORMAT (/' **********WARNING -- DSOS FAILED SOME TESTS**********')
 9030 FORMAT (/' ----------DSOS PASSED ALL TESTS----------')
      END
