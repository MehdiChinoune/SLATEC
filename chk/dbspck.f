*DECK DBSPCK
      SUBROUTINE DBSPCK (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  DBSPCK
C***PURPOSE  Quick check for the B-Spline package.
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (BSPCK-S, DBSPCK-D)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C   DBSPCK is a quick check routine for the B-Spline package which
C   tests consistency between results from higher level routines.
C   Those routines not explicitly called are exercised at some lower
C   level.  The routines exercised are DBFQAD, DBINT4, DBINTK, DBNFAC,
C   DBNSLV, DBSGQ8, DBSPDR, DBSPEV, DBSPPP, DBSPVD, DBSPVN, DBSQAD,
C   DBVALU, DINTRV, DPFQAD, DPPGQ8, DPPQAD and DPPVAL.
C
C***ROUTINES CALLED  D1MACH, DBFQAD, DBINT4, DBINTK, DBSPDR, DBSPEV,
C                    DBSPPP, DBSPVD, DBSPVN, DBSQAD, DBVALU, DFB,
C                    DINTRV, DPFQAD, DPPQAD, DPPVAL
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891004  Removed unreachable code.  (WRB)
C   891009  Removed unreferenced variables.  (WRB)
C   891009  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   930214  Declarations sections added, code revised to test error
C           returns for all values of KPRINT and code polished.  (WRB)
C***END PROLOGUE  DBSPCK
C     .. Scalar Arguments ..
      INTEGER IPASS, KPRINT, LUN
C     .. Local Scalars ..
      DOUBLE PRECISION ATOL, BQUAD, BV, DEN, DN, ER, FBCL, FBCR, PI,
     +     PQUAD, QUAD, SPV, TOL, X1, X2, XL, XX
      INTEGER I, IBCL, IBCR, ID, IERR, IKNT, ILEFT, ILO, INBV, INEV,
     +        INPPV, IWORK, J, JHIGH, JJ, K, KK, KNT, KNTOPT, KONTRL,
     +        LDC, LDCC, LXI, MFLAG, N, NDATA, NERR, NMK, NN
      LOGICAL FATAL
C     .. Local Arrays ..
      DOUBLE PRECISION ADIF(52), BC(13), C(4, 10), CC(4, 4), Q(3),
     +     QQ(77), QSAVE(2), SV(4), T(17), W(65), X(11), XI(11), Y(11)
C     .. External Functions ..
      DOUBLE PRECISION D1MACH, DBVALU, DFB, DPPVAL
      INTEGER NUMXER
C     .. External Subroutines ..
      EXTERNAL DBFQAD, DBINT4, DBINTK, DBSPDR, DBSPEV, DBSPPP, DBSPVD,
     +         DBSPVN, DBSQAD, DFB, DINTRV, DPFQAD, DPPQAD, XGETF, XSETF
C     .. Intrinsic Functions ..
      INTRINSIC ABS, MAX, SIN
C***FIRST EXECUTABLE STATEMENT  DBSPCK
      IF (KPRINT .GE. 2) WRITE (LUN, 9000)
C
      IPASS = 1
      PI = 3.14159265358979324D0
      TOL = 1000.0D0*MAX(D1MACH(4),1.0D-18)
C
C     Generate data.
C
      NDATA = 11
      DEN = NDATA - 1
      DO 20 I = 1,NDATA
        X(I) = (I-1)/DEN
        Y(I) = SIN(PI*X(I))
   20 CONTINUE
      X(3) = 2.0D0/DEN
      Y(3) = SIN(PI*X(3))
C
C     Compute splines for two knot arrays.
C
      DO 110 IKNT = 1,2
        KNT = 3 - IKNT
        IBCL = 1
        IBCR = 2
        FBCL = PI
        FBCR = 0.0D0
        CALL DBINT4 (X,Y,NDATA,IBCL,IBCR,FBCL,FBCR,KNT,T,BC,N,K,W)
C
C       Error test on DBINT4.
C
        INBV = 1
        DO 30 I = 1,NDATA
          XX = X(I)
          BV = DBVALU(T,BC,N,K,0,XX,INBV,W)
          ER = ABS(Y(I)-BV)
          IF (ER .GT. TOL) THEN
            IPASS = 0
            IF (KPRINT .GE. 2) WRITE (LUN, 9010)
          ENDIF
   30   CONTINUE
        INBV = 1
        BV = DBVALU(T,BC,N,K,1,X(1),INBV,W)
        ER = ABS(PI-BV)
        IF (ER .GT. TOL) THEN
          IPASS = 0
          IF (KPRINT .GE. 2) WRITE (LUN, 9020)
        ENDIF
        BV = DBVALU(T,BC,N,K,2,X(NDATA),INBV,W)
        ER = ABS(BV)
        IF (ER .GT. TOL) THEN
          IPASS = 0
          IF (KPRINT .GE. 2) WRITE (LUN, 9030)
        ENDIF
C
C       Test for equality of area from 4 routines.
C
        X1 = X(1)
        X2 = X(NDATA)
        CALL DBSQAD (T,BC,N,K,X1,X2,BQUAD,W)
        LDC = 4
        CALL DBSPPP (T,BC,N,K,LDC,C,XI,LXI,W)
        CALL DPPQAD (LDC,C,XI,LXI,K,X1,X2,Q(1))
        CALL DBFQAD (DFB,T,BC,N,K,0,X1,X2,TOL,Q(2),IERR,W)
        CALL DPFQAD (DFB,LDC,C,XI,LXI,K,0,X1,X2,TOL,Q(3),IERR)
C
C       Error test for quadratures.
C
        DO 90 I = 1,3
          ER = ABS(BQUAD-Q(I))
          IF (ER .GT. TOL) THEN
            IPASS = 0
            IF (KPRINT .GE. 2) WRITE (LUN, 9040)
          ENDIF
   90   CONTINUE
        QSAVE(KNT) = BQUAD
  110 CONTINUE
      ER = ABS(QSAVE(1)-QSAVE(2))
      IF (ER .GT. TOL) THEN
        IPASS = 0
        IF (KPRINT .GE. 2) WRITE (LUN, 9060)
      ENDIF
C
C     Check DBSPDR and DBSPEV against DBVALU, DPPVAL and DBSPVD.
C
      CALL DBSPDR (T,BC,N,K,K,ADIF)
      INEV = 1
      INBV = 1
      INPPV = 1
      ILO = 1
      DO 170 I = 1,6
        XX = X(I+I-1)
        CALL DBSPEV (T,ADIF,N,K,K,XX,INEV,SV,W)
        ATOL = TOL
        DO 130 J = 1,K
          SPV = DBVALU (T,BC,N,K,J-1,XX,INBV,W)
          ER = ABS(SPV-SV(J))
          X2 = ABS(SV(J))
          IF (X2 .GT. 1.0D0) ER = ER/X2
          IF (ER .GT. ATOL) THEN
            IPASS = 0
            IF (KPRINT .GE. 2) WRITE (LUN, 9070)
          ENDIF
          ATOL = 10.0D0*ATOL
  130   CONTINUE
        ATOL = TOL
        DO 140 J = 1,K
          SPV = DPPVAL (LDC,C,XI,LXI,K,J-1,XX,INPPV)
          ER = ABS(SPV-SV(J))
          X2 = ABS(SV(J))
          IF (X2 .GT. 1.0D0) ER = ER/X2
          IF (ER .GT. ATOL) THEN
            IPASS = 0
            IF (KPRINT .GE. 2) WRITE (LUN, 9080)
          ENDIF
          ATOL = 10.0D0*ATOL
  140   CONTINUE
        ATOL = TOL
        LDCC = 4
        X1 = XX
        IF (I+I-1 .EQ. NDATA) X1 = T(N)
        NN = N + K
        CALL DINTRV (T,NN,X1,ILO,ILEFT,MFLAG)
        DO 160 J = 1,K
          CALL DBSPVD (T,K,J,XX,ILEFT,LDCC,CC,W)
          ER = 0.0D0
          DO 150 JJ = 1,K
            ER = ER + BC(ILEFT-K+JJ)*CC(JJ,J)
  150     CONTINUE
          ER = ABS(ER-SV(J))
          X2 = ABS(SV(J))
          IF (X2 .GT. 1.0D0) ER = ER/X2
          IF (ER .GT. ATOL) THEN
            IPASS = 0
            IF (KPRINT .GE. 2) WRITE (LUN, 9090)
          ENDIF
          ATOL = 10.0D0*ATOL
  160   CONTINUE
  170 CONTINUE
      DO 220 K = 2,4
        N = NDATA
        NMK = N - K
        DO 190 I = 1,K
          T(I) = X(1)
          T(N+I) = X(N)
  190   CONTINUE
        XL = X(N) - X(1)
        DN = N - K + 1
        DO 200 I = 1,NMK
          T(K+I) = X(1) + I*XL/DN
  200   CONTINUE
        CALL DBINTK (X,Y,T,N,K,BC,QQ,W)
C
C       Error test on DBINTK.
C
        INBV = 1
        DO 210 I = 1,N
          XX = X(I)
          BV = DBVALU(T,BC,N,K,0,XX,INBV,W)
          ER = ABS(Y(I)-BV)
          IF (ER .GT. TOL) THEN
            IPASS = 0
            IF (KPRINT .GE. 2) WRITE (LUN, 9100)
          ENDIF
  210   CONTINUE
  220 CONTINUE
C
C     Trigger error conditions.
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
      IF (KPRINT .GE. 3) WRITE (LUN, 9050)
C
      W(1) = 11.0D0
      W(2) = 4.0D0
      W(3) = 2.0D0
      W(4) = 0.5D0
      W(5) = 4.0D0
      ILO = 1
      INEV = 1
      INBV = 1
      CALL DINTRV (T,N+1,W(4),ILO,ILEFT,MFLAG)
      DO 320 I = 1,5
        W(I) = -W(I)
        N = W(1)
        K = W(2)
        ID = W(3)
        XX = W(4)
        LDC = W(5)
        IF (I .LE. 4) THEN
          BV = DBVALU (T,BC,N,K,ID,XX,INBV,QQ)
          IF (NUMXER(NERR) .NE. 2) THEN
            IPASS = 0
            FATAL = .TRUE.
          ENDIF
          CALL XERCLR
C
          CALL DBSPEV (T,ADIF,N,K,ID,XX,INEV,SV,QQ)
          IF (NUMXER(NERR) .NE. 2) THEN
            IPASS = 0
            FATAL = .TRUE.
          ENDIF
          CALL XERCLR
C
          JHIGH = N - 10
          CALL DBSPVN (T,JHIGH,K,ID,XX,ILEFT,SV,QQ,IWORK)
          IF (NUMXER(NERR) .NE. 2) THEN
            IPASS = 0
            FATAL = .TRUE.
          ENDIF
          CALL XERCLR
C
          CALL DBFQAD (DFB,T,BC,N,K,ID,XX,X2,TOL,QUAD,IERR,QQ)
          IF (NUMXER(NERR) .NE. 2) THEN
            IPASS = 0
            FATAL = .TRUE.
          ENDIF
          CALL XERCLR
        ENDIF
C
        IF (I.NE.3 .AND. I.NE.4) THEN
          CALL DBSPPP (T,BC,N,K,LDC,C,XI,LXI,QQ)
          IF (NUMXER(NERR) .NE. 2) THEN
            IPASS = 0
            FATAL = .TRUE.
          ENDIF
          CALL XERCLR
        ENDIF
C
        IF (I .LE. 3) THEN
          CALL DBSPDR (T,BC,N,K,ID,ADIF)
          IF (NUMXER(NERR) .NE. 2) THEN
            IPASS = 0
            FATAL = .TRUE.
          ENDIF
          CALL XERCLR
        ENDIF
C
        IF (I.NE.3 .AND. I.NE.5) THEN
          CALL DBSQAD (T,BC,N,K,XX,X2,BQUAD,QQ)
          IF (NUMXER(NERR) .NE. 2) THEN
            IPASS = 0
            FATAL = .TRUE.
          ENDIF
          CALL XERCLR
        ENDIF
C
        IF (I .GT. 1) THEN
          CALL DBSPVD (T,K,ID,XX,ILEFT,LDC,C,QQ)
          IF (NUMXER(NERR) .NE. 2) THEN
            IPASS = 0
            FATAL = .TRUE.
          ENDIF
          CALL XERCLR
        ENDIF
C
        IF (I .LE. 2) THEN
          CALL DBINTK (X,Y,T,N,K,BC,QQ,ADIF)
          IF (NUMXER(NERR) .NE. 2) THEN
            IPASS = 0
            FATAL = .TRUE.
          ENDIF
          CALL XERCLR
        ENDIF
C
        IF (I .NE. 4) THEN
          KNTOPT = LDC - 2
          IBCL = K - 2
          CALL DBINT4 (X,Y,N,IBCL,ID,FBCL,FBCR,KNTOPT,T,BC,NN,KK,QQ)
          IF (NUMXER(NERR) .NE. 2) THEN
            IPASS = 0
            FATAL = .TRUE.
          ENDIF
          CALL XERCLR
        ENDIF
        W(I) = -W(I)
  320 CONTINUE
      KNTOPT = 1
      X(1) = 1.0D0
      CALL DBINT4 (X,Y,N,IBCL,IBCR,FBCL,FBCR,KNTOPT,T,BC,N,K,QQ)
      IF (NUMXER(NERR) .NE. 2) THEN
        IPASS = 0
        FATAL = .TRUE.
      ENDIF
      CALL XERCLR
C
      CALL DBINTK (X,Y,T,N,K,BC,QQ,ADIF)
      IF (NUMXER(NERR) .NE. 2) THEN
        IPASS = 0
        FATAL = .TRUE.
      ENDIF
      CALL XERCLR
C
      X(1) = 0.0D0
      ATOL = 1.0D0
      KNTOPT = 3
      DO 330 I = 1,3
        QQ(I) = -0.30D0 + 0.10D0*(I-1)
        QQ(I+3) = 1.1D0 + 0.10D0*(I-1)
  330 CONTINUE
      QQ(1) = 1.0D0
      CALL DBINT4 (X,Y,NDATA,1,1,FBCL,FBCR,3,T,BC,N,K,QQ)
      IF (NUMXER(NERR) .NE. 2) THEN
        IPASS = 0
        FATAL = .TRUE.
      ENDIF
      CALL XERCLR
C
      CALL DBFQAD (DFB,T,BC,N,K,ID,X1,X2,ATOL,QUAD,IERR,QQ)
      IF (NUMXER(NERR) .NE. 2) THEN
        IPASS = 0
        FATAL = .TRUE.
      ENDIF
      CALL XERCLR
C
      INPPV = 1
      DO 350 I = 1,5
        W(I) = -W(I)
        LXI = W(1)
        K = W(2)
        ID = W(3)
        XX = W(4)
        LDC = W(5)
        SPV = DPPVAL (LDC,C,XI,LXI,K,ID,XX,INPPV)
        IF ((I.NE.4.AND.NUMXER(NERR).NE.2) .OR.
     +      (I.EQ.4.AND.NUMXER(NERR).NE.0)) THEN
          IPASS = 0
          FATAL = .TRUE.
        ENDIF
        CALL XERCLR
C
        CALL DPFQAD (DFB,LDC,C,XI,LXI,K,ID,XX,X2,TOL,QUAD,IERR)
        IF ((I.NE.4.AND.NUMXER(NERR).NE.2) .OR.
     +      (I.EQ.4.AND.NUMXER(NERR).NE.0)) THEN
          IPASS = 0
          FATAL = .TRUE.
        ENDIF
        CALL XERCLR
C
        IF (I .NE. 3) THEN
          CALL DPPQAD (LDC,C,XI,LXI,K,XX,X2,PQUAD)
          IF ((I.NE.4.AND.NUMXER(NERR).NE.2) .OR.
     +        (I.EQ.4.AND.NUMXER(NERR).NE.0)) THEN
            IPASS = 0
            FATAL = .TRUE.
          ENDIF
          CALL XERCLR
        ENDIF
C
        W(I) = -W(I)
  350 CONTINUE
      LDC = W(5)
      CALL DPFQAD (DFB,LDC,C,XI,LXI,K,ID,X1,X2,ATOL,QUAD,IERR)
      IF (NUMXER(NERR) .NE. 2) THEN
        IPASS = 0
        FATAL = .TRUE.
      ENDIF
      CALL XERCLR
C
C     Restore KONTRL and check to see if the tests of error detection
C     passed.
C
      CALL XSETF (KONTRL)
      IF (FATAL) THEN
         IF (KPRINT .GE. 2) THEN
            WRITE (LUN, 9110)
         ENDIF
      ELSE
         IF (KPRINT .GE. 3) THEN
            WRITE (LUN, 9120)
         ENDIF
      ENDIF
C
C     Print PASS/FAIL message.
C
      IF (IPASS.EQ.1 .AND. KPRINT.GE.2) WRITE (LUN, 9200)
      IF (IPASS.EQ.0 .AND. KPRINT.GE.1) WRITE (LUN, 9210)
      RETURN
C
 9000 FORMAT ('1 QUICK CHECK FOR SPLINE ROUTINES',//)
 9010 FORMAT (' ERROR TEST FOR INTERPOLATION BY DBINT4 NOT SATISFIED')
 9020 FORMAT (' ERROR TEST FOR INTERPOLATION BY DBINT4 NOT SATISFIED ',
     +        'BY FIRST DERIVATIVE')
 9030 FORMAT (' ERROR TEST FOR INTERPOLATION BY DBINT4 NOT SATISFIED ',
     +        'BY SECOND DERIVATIVE')
 9040 FORMAT (' ERROR IN QUADRATURE CHECKS')
 9050 FORMAT (/' TRIGGER 52 ERROR CONDITIONS',/)
 9060 FORMAT (' ERROR IN QUADRATURE CHECK USING TWO SETS OF KNOTS')
 9070 FORMAT (' COMPARISONS FROM DBSPEV AND DBVALU DO NOT AGREE')
 9080 FORMAT (' COMPARISONS FROM DBSPEV AND DPPVAL DO NOT AGREE')
 9090 FORMAT (' COMPARISONS FROM DBSPEV AND DBSPVD DO NOT AGREE')
 9100 FORMAT (' ERROR TEST FOR INTERPOLATION BY DBINTK NOT SATISFIED')
 9110 FORMAT (/ ' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
 9120 FORMAT (/ ' ALL INCORRECT ARGUMENT TESTS PASSED')
 9200 FORMAT (/' **********B-SPLINE PACKAGE PASSED ALL TESTS**********')
 9210 FORMAT (/' *********B-SPLINE PACKAGE FAILED SOME TESTS**********')
      END
