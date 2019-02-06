*DECK SNSQQK
      SUBROUTINE SNSQQK (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  SNSQQK
C***PURPOSE  Quick check for SNSQE and SNSQ.
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (SNSQQK-S, DNSQQK-D)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C   This subroutine performs a quick check on the subroutine SNSQE
C   (and SNSQ).
C
C***ROUTINES CALLED  ENORM, PASS, R1MACH, SNSQE, SQFCN2, SQJAC2
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   891009  Removed unreferenced variable.  (WRB)
C   891009  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Code cleaned up and TYPE section added.  (RWC, WRB)
C***END PROLOGUE  SNSQQK
C     .. Scalar Arguments ..
      INTEGER IPASS, KPRINT, LUN
C     .. Local Scalars ..
      REAL FNORM, FNORMS, TOL
      INTEGER ICNT, INFO, INFOS, IOPT, LWA, N, NPRINT
C     .. Local Arrays ..
      REAL FVEC(2), WA(19), X(2)
      INTEGER ITEST(3)
C     .. External Functions ..
      REAL ENORM, R1MACH
      EXTERNAL ENORM, R1MACH
C     .. External Subroutines ..
      EXTERNAL PASS, SNSQE, SQFCN2, SQJAC2
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C***FIRST EXECUTABLE STATEMENT  SNSQQK
      INFOS = 1
      FNORMS = 0.0E0
      N = 2
      LWA = 19
      NPRINT = -1
      TOL = SQRT(R1MACH(4))
      IF (KPRINT .GE. 2) WRITE (LUN,9000)
C
C     Option 1, the user provides the Jacobian.
C
      IOPT = 1
      X(1) = -1.2E0
      X(2) = 1.0E0
      CALL SNSQE (SQFCN2,SQJAC2,IOPT,N,X,FVEC,TOL,NPRINT,INFO,WA,LWA)
      ICNT = 1
      FNORM = ENORM(N,FVEC)
      ITEST(ICNT) = 0
      IF ((INFO.EQ.INFOS) .AND. (FNORM-FNORMS.LE.TOL)) ITEST(ICNT) = 1
C
      IF (KPRINT .NE. 0) THEN
         IF ((KPRINT.GE.2 .AND. ITEST(ICNT).NE.1) .OR. KPRINT.GE.3)
     +       WRITE (LUN,9010) INFOS,FNORMS,INFO,FNORM
         IF ((KPRINT.GE.2) .OR. (KPRINT.EQ.1 .AND. ITEST(ICNT).NE.1))
     +       CALL PASS (LUN, ICNT, ITEST(ICNT))
      ENDIF
C
C     Option 2, the code approximates the Jacobian.
C
      IOPT = 2
      X(1) = -1.2E0
      X(2) = 1.0E0
      CALL SNSQE (SQFCN2,SQJAC2,IOPT,N,X,FVEC,TOL,NPRINT,INFO,WA,LWA)
      ICNT = 2
      FNORM = ENORM(N,FVEC)
      ITEST(ICNT) = 0
      IF ((INFO.EQ.INFOS) .AND. (FNORM-FNORMS.LE.TOL)) ITEST(ICNT) = 1
C
      IF (KPRINT .NE. 0) THEN
         IF (KPRINT.GE.3 .OR. (KPRINT.GE.2.AND.ITEST(ICNT).NE.1))
     +       WRITE (LUN,9010) INFOS, FNORMS, INFO, FNORM
         IF (KPRINT.GE.2 .OR. (KPRINT.EQ.1.AND.ITEST(ICNT).NE.1))
     +       CALL PASS (LUN, ICNT, ITEST(ICNT))
      ENDIF
C
C     Test improper input parameters.
C
      LWA = 15
      IOPT = 1
      X(1) = -1.2E0
      X(2) = 1.0E0
      CALL SNSQE (SQFCN2,SQJAC2,IOPT,N,X,FVEC,TOL,NPRINT,INFO,WA,LWA)
      ICNT = 3
      ITEST(ICNT) = 0
      IF (INFO .EQ. 0) ITEST(ICNT) = 1
      IF (KPRINT.GE.2 .OR. (KPRINT.EQ.1.AND.ITEST(ICNT).NE.1))
     +    CALL PASS (LUN, ICNT, ITEST(ICNT))
C
C     Set IPASS.
C
      IPASS = ITEST(1)*ITEST(2)*ITEST(3)
      IF (KPRINT.GE.1 .AND. IPASS.NE.1) WRITE (LUN,9020)
      IF (KPRINT.GE.2 .AND. IPASS.EQ.1) WRITE (LUN,9030)
      RETURN
 9000 FORMAT ('1' / '  SNSQE QUICK CHECK'/)
 9010 FORMAT (' EXPECTED VALUE OF INFO AND RESIDUAL NORM', I5, E20.5 /
     +        ' RETURNED VALUE OF INFO AND RESIDUAL NORM', I5, E20.5 /)
 9020 FORMAT (/' **********WARNING -- SNSQE/SNSQ FAILED SOME TESTS****',
     +        '******')
 9030 FORMAT (/' ----------SNSQE/SNSQ PASSED ALL TESTS----------')
      END
