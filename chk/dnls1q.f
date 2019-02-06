*DECK DNLS1Q
      SUBROUTINE DNLS1Q (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  DNLS1Q
C***PURPOSE  Quick check for DNLS1E, DNLS1 and DCOV.
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (SNLS1Q-S, DNLS1Q-D)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C   This subroutine performs a quick check on the subroutines DNLS1E
C   (and DNLS1) and DCOV.
C
C***ROUTINES CALLED  DENORM, DFCN1, DFCN2, DFCN3, DFDJC3, PASS, D1MACH,
C                    DCOV, DNLS1E
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   930214  Declarations sections added, code revised to test error
C           returns for all values of KPRINT and code polished.  (WRB)
C***END PROLOGUE  DNLS1Q
C     .. Scalar Arguments ..
      INTEGER IPASS, KPRINT, LUN
C     .. Local Scalars ..
      DOUBLE PRECISION FNORM, FNORMS, ONE, SIGMA, TEMP1, TEMP2, TEMP3,
     +                 TOL, TOL2, ZERO
      INTEGER I, IFLAG, INFO, INFOS, IOPT, KONTRL, LDFJAC, LWA, M, N,
     +        NERR, NPRINT
      LOGICAL FATAL
C     .. Local Arrays ..
      DOUBLE PRECISION FJAC(10,2), FJROW(2), FJTJ(3), FVEC(10), WA(40),
     +                 X(2)
      INTEGER IW(2)
C     .. External Functions ..
      DOUBLE PRECISION D1MACH, DENORM
      INTEGER NUMXER
      EXTERNAL D1MACH, DENORM, NUMXER
C     .. External Subroutines ..
      EXTERNAL DFCN1, DFCN2, DFCN3, DFDJC3, PASS, DCOV, DNLS1E, XGETF,
     +         XSETF
C     .. Intrinsic Functions ..
      INTRINSIC ABS, SQRT
C***FIRST EXECUTABLE STATEMENT  DNLS1Q
      IF (KPRINT .GE. 2) WRITE (LUN,9000)
C
      IPASS = 1
      INFOS = 1
      FNORMS = 1.1151779D+01
      M = 10
      N = 2
      LWA = 40
      LDFJAC = 10
      NPRINT = -1
      IFLAG = 1
      ZERO = 0.0D0
      ONE = 1.0D0
      TOL = MAX(SQRT(40.0D0*D1MACH(4)),1.0D-12)
      TOL2 = SQRT(TOL)
C
C     OPTION=2, the full Jacobian is stored and the user provides the
C     Jacobian.
C
      IOPT = 2
      X(1) = 3.0D-1
      X(2) = 4.0D-1
      CALL DNLS1E(DFCN2,IOPT,M,N,X,FVEC,TOL,NPRINT,INFO,IW,WA,LWA)
      FNORM = DENORM(M,FVEC)
      IF (INFO.EQ.INFOS .AND. ABS(FNORM-FNORMS)/FNORMS.LE.TOL2) THEN
        FATAL = .FALSE.
        IF (KPRINT .GE. 3) CALL PASS(LUN,1,1)
      ELSE
        IPASS = 0
        FATAL = .TRUE.
        IF (KPRINT .GE. 2) CALL PASS(LUN,1,0)
      ENDIF
      IF ((FATAL.AND.KPRINT.GE.2) .OR. KPRINT.GE.3)
     +   WRITE (LUN,9010) INFOS,FNORMS,INFO,FNORM
C
C     Form JAC-transpose*JAC.
C
      SIGMA = FNORM*FNORM/(M-N)
      IFLAG = 2
      CALL DFCN2(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
      DO 20 I = 1,3
        FJTJ(I) = ZERO
   20 CONTINUE
      DO 30 I = 1,M
        FJTJ(1) = FJTJ(1) + FJAC(I,1)**2
        FJTJ(2) = FJTJ(2) + FJAC(I,1)*FJAC(I,2)
        FJTJ(3) = FJTJ(3) + FJAC(I,2)**2
   30 CONTINUE
C
C     Calculate the covariance matrix.
C
      CALL DCOV(DFCN2,IOPT,M,N,X,FVEC,FJAC,LDFJAC,INFO,WA(1),WA(N+1),
     +          WA(2*N+1),WA(3*N+1))
C
C     Form JAC-transpose*JAC * covariance matrix (should = SIGMA*I).
C
      TEMP1 = (FJTJ(1)*FJAC(1,1)+FJTJ(2)*FJAC(1,2))/SIGMA
      TEMP2 = (FJTJ(1)*FJAC(1,2)+FJTJ(2)*FJAC(2,2))/SIGMA
      TEMP3 = (FJTJ(2)*FJAC(1,2)+FJTJ(3)*FJAC(2,2))/SIGMA
      IF (INFO.EQ.INFOS .AND. ABS(TEMP1-ONE).LT.TOL2 .AND.
     +    ABS(TEMP2).LT.TOL2 .AND. ABS(TEMP3-ONE).LT.TOL2) THEN
        FATAL = .FALSE.
        IF (KPRINT .GE. 3) CALL PASS(LUN,2,1)
      ELSE
        IPASS = 0
        FATAL = .TRUE.
        IF (KPRINT .GE. 2) CALL PASS(LUN,2,0)
      ENDIF
      IF ((FATAL.AND.KPRINT.GE.2) .OR. KPRINT.GE.3)
     +  WRITE (LUN,9020) INFOS,INFO,TEMP1,TEMP2,TEMP3
C
C     OPTION=1, the full Jacobian is stored and the code approximates
C     the Jacobian.
C
      IOPT = 1
      X(1) = 3.0D-1
      X(2) = 4.0D-1
      CALL DNLS1E(DFCN1,IOPT,M,N,X,FVEC,TOL,NPRINT,INFO,IW,WA,LWA)
      FNORM = DENORM(M,FVEC)
      IF (INFO.EQ.INFOS .AND. ABS(FNORM-FNORMS)/FNORMS.LE.TOL2) THEN
        FATAL = .FALSE.
        IF (KPRINT .GE. 3) CALL PASS(LUN,3,1)
      ELSE
        IPASS = 0
        FATAL = .TRUE.
        IF (KPRINT .GE. 2) CALL PASS(LUN,3,0)
      ENDIF
      IF ((FATAL.AND.KPRINT.GE.2) .OR. KPRINT.GE.3)
     +  WRITE (LUN,9010) INFOS,FNORMS,INFO,FNORM
C
C     Form JAC-transpose*JAC.
C
      SIGMA = FNORM*FNORM/(M-N)
      IFLAG = 1
      CALL DFDJC3(DFCN1,M,N,X,FVEC,FJAC,LDFJAC,IFLAG,ZERO,WA)
      DO 60 I = 1,3
        FJTJ(I) = ZERO
   60 CONTINUE
      DO 70 I = 1,M
        FJTJ(1) = FJTJ(1) + FJAC(I,1)**2
        FJTJ(2) = FJTJ(2) + FJAC(I,1)*FJAC(I,2)
        FJTJ(3) = FJTJ(3) + FJAC(I,2)**2
   70 CONTINUE
C
C     Calculate the covariance matrix.
C
      CALL DCOV(DFCN1,IOPT,M,N,X,FVEC,FJAC,LDFJAC,INFO,WA(1),WA(N+1),
     +          WA(2*N+1),WA(3*N+1))
C
C     Form JAC-transpose*JAC * covariance matrix (should = SIGMA*I).
C
      TEMP1 = (FJTJ(1)*FJAC(1,1)+FJTJ(2)*FJAC(1,2))/SIGMA
      TEMP2 = (FJTJ(1)*FJAC(1,2)+FJTJ(2)*FJAC(2,2))/SIGMA
      TEMP3 = (FJTJ(2)*FJAC(1,2)+FJTJ(3)*FJAC(2,2))/SIGMA
      IF (INFO.EQ.INFOS .AND. ABS(TEMP1-ONE).LT.TOL2 .AND.
     +    ABS(TEMP2).LT.TOL2 .AND. ABS(TEMP3-ONE).LT.TOL2) THEN
        FATAL = .FALSE.
        IF (KPRINT .GE. 3) CALL PASS(LUN,4,1)
      ELSE
        IPASS = 0
        FATAL = .TRUE.
        IF (KPRINT .GE. 2) CALL PASS(LUN,4,0)
      ENDIF
      IF ((FATAL.AND.KPRINT.GE.2) .OR. KPRINT.GE.3)
     +  WRITE (LUN,9020) INFOS,INFO,TEMP1,TEMP2,TEMP3
C
C     OPTION=3, the full Jacobian is not stored.  Only the product of
C     the Jacobian transpose and Jacobian is stored.  The user provides
C     the Jacobian one row at a time.
C
      IOPT = 3
      X(1) = 3.0D-1
      X(2) = 4.0D-1
      CALL DNLS1E(DFCN3,IOPT,M,N,X,FVEC,TOL,NPRINT,INFO,IW,WA,LWA)
      FNORM = DENORM(M,FVEC)
      IF (INFO.EQ.INFOS .AND. ABS(FNORM-FNORMS)/FNORMS.LE.TOL2) THEN
        FATAL = .FALSE.
        IF (KPRINT .GE. 3) CALL PASS(LUN,5,1)
      ELSE
        IPASS = 0
        FATAL = .TRUE.
        IF (KPRINT .GE. 2) CALL PASS(LUN,5,0)
      ENDIF
      IF ((FATAL.AND.KPRINT.GE.2) .OR. KPRINT.GE.3)
     +  WRITE (LUN,9010) INFOS,FNORMS,INFO,FNORM
C
C     Form JAC-transpose*JAC.
C
      SIGMA = FNORM*FNORM/(M-N)
      DO 100 I = 1,3
        FJTJ(I) = ZERO
  100 CONTINUE
      IFLAG = 3
      DO 110 I = 1,M
        CALL DFCN3(IFLAG,M,N,X,FVEC,FJROW,I)
        FJTJ(1) = FJTJ(1) + FJROW(1)**2
        FJTJ(2) = FJTJ(2) + FJROW(1)*FJROW(2)
        FJTJ(3) = FJTJ(3) + FJROW(2)**2
  110 CONTINUE
C
C     Calculate the covariance matrix.
C
      CALL DCOV(DFCN3,IOPT,M,N,X,FVEC,FJAC,LDFJAC,INFO,WA(1),WA(N+1),
     +          WA(2*N+1),WA(3*N+1))
C
C     Form JAC-transpose*JAC * covariance matrix (should = SIGMA*I).
C
      TEMP1 = (FJTJ(1)*FJAC(1,1)+FJTJ(2)*FJAC(1,2))/SIGMA
      TEMP2 = (FJTJ(1)*FJAC(1,2)+FJTJ(2)*FJAC(2,2))/SIGMA
      TEMP3 = (FJTJ(2)*FJAC(1,2)+FJTJ(3)*FJAC(2,2))/SIGMA
      IF (INFO.EQ.INFOS .AND. ABS(TEMP1-ONE).LT.TOL2 .AND.
     +    ABS(TEMP2).LT.TOL2 .AND. ABS(TEMP3-ONE).LT.TOL2) THEN
        FATAL = .FALSE.
        IF (KPRINT .GE. 3) CALL PASS(LUN,6,1)
      ELSE
        IPASS = 0
        FATAL = .TRUE.
        IF (KPRINT .GE. 2) CALL PASS(LUN,6,0)
      ENDIF
      IF ((FATAL.AND.KPRINT.GE.2) .OR. KPRINT.GE.3)
     +  WRITE (LUN,9020) INFOS,INFO,TEMP1,TEMP2,TEMP3
C
C     Test improper input parameters.
C
      CALL XGETF (KONTRL)
      IF (KPRINT .LE. 2) THEN
         CALL XSETF (0)
      ELSE
         CALL XSETF (1)
      ENDIF
      FATAL = .FALSE.
      CALL XERCLR
C
      IF (KPRINT .GE. 3) WRITE (LUN, 9060)
C
      LWA = 35
      IOPT = 2
      X(1) = 3.0D-1
      X(2) = 4.0D-1
      CALL DNLS1E(DFCN2,IOPT,M,N,X,FVEC,TOL,NPRINT,INFO,IW,WA,LWA)
      IF (INFO.NE.0 .OR. NUMXER(NERR).NE.2) FATAL = .TRUE.
C
      M = 0
      CALL DCOV(DFCN2,IOPT,M,N,X,FVEC,FJAC,LDFJAC,INFO,WA(1),WA(N+1),
     +          WA(2*N+1),WA(3*N+1))
      IF (INFO.NE.0 .OR. NUMXER(NERR).NE.2) FATAL = .TRUE.
C
C     Restore KONTRL and check to see if the tests of error detection
C     passed.
C
      CALL XSETF (KONTRL)
      IF (FATAL) THEN
         IPASS = 0
         IF (KPRINT .GE. 2) THEN
            WRITE (LUN, 9070)
         ENDIF
      ELSE
         IF (KPRINT .GE. 3) THEN
            WRITE (LUN, 9080)
         ENDIF
      ENDIF
C
C     Print PASS/FAIL message.
C
      IF (IPASS.EQ.1 .AND. KPRINT.GE.2) WRITE (LUN,9100)
      IF (IPASS.EQ.0 .AND. KPRINT.GE.1) WRITE (LUN,9110)
C
  130 RETURN
C
 9000 FORMAT ('1' / ' Test DNLS1E, DNLS1 and DCOV')
 9010 FORMAT (' EXPECTED VALUE OF INFO AND RESIDUAL NORM', I5, D20.9 /
     +        ' RETURNED VALUE OF INFO AND RESIDUAL NORM', I5, D20.9 /)
 9020 FORMAT (' EXPECTED AND RETURNED VALUE OF INFO', I5, 10X, I5 /
     +        ' RETURNED PRODUCT OF (J-TRANS*J)*COVARIANCE MATRIX/SIGMA'
     +        / ' (SHOULD = THE IDENTITY -- 1.0, 0.0, 1.0)' / 3D20.9 /)
 9060 FORMAT (/ ' TRIGGER 2 ERROR MESSAGES',/)
 9070 FORMAT (' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
 9080 FORMAT (' ALL INCORRECT ARGUMENT TESTS PASSED')
 9100 FORMAT (/' *************DNLS1E PASSED ALL TESTS*****************')
 9110 FORMAT (/' ************DNLS1E FAILED SOME TESTS*****************')
      END
