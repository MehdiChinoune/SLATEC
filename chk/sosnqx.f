*DECK SOSNQX
      SUBROUTINE SOSNQX (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  SOSNQX
C***PURPOSE  Quick check for SOS.
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (SOSNQX-S, DSOSQX-D)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C   This subroutine performs a quick check on the subroutine SOS.
C
C***ROUTINES CALLED  PASS, R1MACH, SNRM2, SOS, SOSFNC
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890618  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Code cleaned up and TYPE section added.  (RWC, WRB)
C***END PROLOGUE  SOSNQX
C     .. Scalar Arguments ..
      INTEGER IPASS, KPRINT, LUN
C     .. Local Scalars ..
      REAL AER, FNORM, FNORMS, RER, TOLF
      INTEGER ICNT, IFLAG, IFLAGS, LIW, LWA, N
C     .. Local Arrays ..
      REAL FVEC(2), WA(17), X(2)
      INTEGER ITEST(2), IW(6)
C     .. External Functions ..
      REAL R1MACH, SNRM2, SOSFNC
      EXTERNAL R1MACH, SNRM2, SOSFNC
C     .. External Subroutines ..
      EXTERNAL PASS, SOS
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C***FIRST EXECUTABLE STATEMENT  SOSNQX
      IFLAGS = 3
      FNORMS = 0.0E0
      N = 2
      LWA = 17
      LIW = 6
      TOLF = SQRT(R1MACH(4))
      RER = SQRT(R1MACH(4))
      AER = 0.0E0
      IF (KPRINT .GE. 2) WRITE (LUN,9000)
C
C     Test the code with proper input values.
C
      IFLAG = 0
      X(1) = -1.2E0
      X(2) = 1.0E0
      CALL SOS (SOSFNC,N,X,RER,AER,TOLF,IFLAG,WA,LWA,IW,LIW)
      ICNT = 1
      FVEC(1) = SOSFNC(X,1)
      FVEC(2) = SOSFNC(X,2)
      FNORM = SNRM2(N,FVEC,1)
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
      X(1) = -1.2E0
      X(2) = 1.0E0
      CALL SOS (SOSFNC,N,X,RER,AER,TOLF,IFLAG,WA,LWA,IW,LIW)
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
 9000 FORMAT ('1' / '  SOS QUICK CHECK' /)
 9010 FORMAT (' EXPECTED VALUE OF IFLAG AND RESIDUAL NORM', I5, E20.5 /
     +        ' RETURNED VALUE OF IFLAG AND RESIDUAL NORM', I5, E20.5 /)
 9020 FORMAT (/' **********WARNING -- SOS FAILED SOME TESTS**********')
 9030 FORMAT (/' ----------SOS PASSED ALL TESTS----------')
      END
