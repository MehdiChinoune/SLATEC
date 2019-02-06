*DECK FFTQX
      SUBROUTINE FFTQX (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  FFTQX
C***PURPOSE  Quick check for the NCAR FFT routines.
C***LIBRARY   SLATEC
C***KEYWORDS  QUICK CHECK
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***DESCRIPTION
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C                       VERSION 4  APRIL 1985
C
C                         A TEST DRIVER FOR
C          A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE FAST FOURIER
C           TRANSFORM OF PERIODIC AND OTHER SYMMETRIC SEQUENCES
C
C                              BY
C
C                       PAUL N SWARZTRAUBER
C
C    NATIONAL CENTER FOR ATMOSPHERIC RESEARCH  BOULDER, COLORADO 80307
C
C        WHICH IS SPONSORED BY THE NATIONAL SCIENCE FOUNDATION
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C
C             THIS PROGRAM TESTS THE PACKAGE OF FAST FOURIER
C     TRANSFORMS FOR BOTH COMPLEX AND REAL PERIODIC SEQUENCES AND
C     CERTAIN OTHER SYMMETRIC SEQUENCES THAT ARE LISTED BELOW.
C
C     1.   RFFTI     INITIALIZE  RFFTF AND RFFTB
C     2.   RFFTF     FORWARD TRANSFORM OF A REAL PERIODIC SEQUENCE
C     3.   RFFTB     BACKWARD TRANSFORM OF A REAL COEFFICIENT ARRAY
C
C     4.   EZFFTI    INITIALIZE EZFFTF AND EZFFTB
C     5.   EZFFTF    A SIMPLIFIED REAL PERIODIC FORWARD TRANSFORM
C     6.   EZFFTB    A SIMPLIFIED REAL PERIODIC BACKWARD TRANSFORM
C
C     7.   SINTI     INITIALIZE SINT
C     8.   SINT      SINE TRANSFORM OF A REAL ODD SEQUENCE
C
C     9.   COSTI     INITIALIZE COST
C     10.  COST      COSINE TRANSFORM OF A REAL EVEN SEQUENCE
C
C     11.  SINQI     INITIALIZE SINQF AND SINQB
C     12.  SINQF     FORWARD SINE TRANSFORM WITH ODD WAVE NUMBERS
C     13.  SINQB     UNNORMALIZED INVERSE OF SINQF
C
C     14.  COSQI     INITIALIZE COSQF AND COSQB
C     15.  COSQF     FORWARD COSINE TRANSFORM WITH ODD WAVE NUMBERS
C     16.  COSQB     UNNORMALIZED INVERSE OF COSQF
C
C     17.  CFFTI     INITIALIZE CFFTF AND CFFTB
C     18.  CFFTF     FORWARD TRANSFORM OF A COMPLEX PERIODIC SEQUENCE
C     19.  CFFTB     UNNORMALIZED INVERSE OF CFFTF
C
C***ROUTINES CALLED  CFFTB, CFFTF, CFFTI, COSQB, COSQF, COSQI, COST,
C                    COSTI, EZFFTB, EZFFTF, EZFFTI, PIMACH, R1MACH,
C                    RFFTB, RFFTF, RFFTI, SINQB, SINQF, SINQI, SINT,
C                    SINTI
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   890718  Changed computation of PI to use PIMACH.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   901205  Changed usage of R1MACH(3) to R1MACH(4).  (RWC)
C   910708  Minor modifications in use of KPRINT.  (WRB)
C   920211  Code cleaned up, an error in printing an error message fixed
C           and comments on PASS/FAIL of individual tests added.  (WRB)
C   920618  Code upgraded to "Version 4".  (BKS, WRB)
C   930315  Modified RFFT* tests to compute "slow-transform" in DOUBLE
C           PRECISION.  (WRB)
C***END PROLOGUE  FFTQX
C     .. Scalar Arguments ..
      INTEGER IPASS, KPRINT, LUN
C     .. Local Scalars ..
      DOUBLE PRECISION ARG, ARG1, ARG2, DT, PI, SUM, SUM1, SUM2
      REAL AZERO, AZEROH, CF, COSQBT, COSQFB, COSQFT, COSTFB, COSTT,
     +     DCFB, DCFFTB, DCFFTF, DEZB1, DEZF1, DEZFB, ERRMAX, RFTB,
     +     RFTF, RFTFB, SIGN, SINQBT, SINQFB, SINQFT, SINTFB, SINTT,
     +     SQRT2, TPI
      INTEGER I, J, K, MODN, N, NM1, NNS, NP1, NS2, NS2M, NZ
C     .. Local Arrays ..
      COMPLEX CX(200), CY(200)
      REAL A(100), AH(100), B(100), BH(100), W(2000), X(200), XH(200),
     +     Y(200)
      INTEGER ND(10)
C     .. External Functions ..
      REAL R1MACH
      EXTERNAL R1MACH
C     .. External Subroutines ..
      EXTERNAL CFFTB, CFFTF, CFFTI, COSQB, COSQF, COSQI, COST, COSTI,
     +         EZFFTB, EZFFTF, EZFFTI, RFFTB, RFFTF, RFFTI, SINQB,
     +         SINQF, SINQI, SINT, SINTI
C     .. Intrinsic Functions ..
      INTRINSIC ABS, CABS, CMPLX, COS, MAX, MOD, SIN, SQRT
C     .. Data statements ..
      DATA ND(1), ND(2), ND(3), ND(4), ND(5), ND(6), ND(7)/120, 54, 49,
     +     32, 4, 3, 2/
C***FIRST EXECUTABLE STATEMENT  FFTQX
      SQRT2 = SQRT(2.0)
      ERRMAX = 2.0*SQRT(R1MACH(4))
      NNS = 7
      PI = 4.0D0*ATAN(1.0D0)
      IF (KPRINT .GE. 2) WRITE (LUN, 9000)
      IPASS = 1
      DO 660 NZ=1,NNS
        N = ND(NZ)
        IF (KPRINT .GE. 2) WRITE (LUN, 9010) N
        MODN = MOD(N, 2)
        NP1 = N + 1
        NM1 = N - 1
        DO 100 J=1,NP1
          X(J) = SIN(J*SQRT2)
          Y(J) = X(J)
          XH(J) = X(J)
  100   CONTINUE
C
C       Test Subroutines RFFTI, RFFTF and RFFTB
C
        CALL RFFTI(N, W)
        DT = (PI+PI)/N
        NS2 = (N+1)/2
        IF (NS2 .LT. 2) GO TO 130
        DO 120 K=2,NS2
          SUM1 = 0.0D0
          SUM2 = 0.0D0
          ARG = (K-1)*DT
          DO 110 I=1,N
            ARG1 = (I-1)*ARG
            SUM1 = SUM1 + X(I)*COS(ARG1)
            SUM2 = SUM2 + X(I)*SIN(ARG1)
  110     CONTINUE
          Y(2*K-2) = SUM1
          Y(2*K-1) = -SUM2
  120   CONTINUE
  130   SUM1 = 0.0D0
        SUM2 = 0.0D0
        DO 140 I=1,NM1,2
          SUM1 = SUM1 + X(I)
          SUM2 = SUM2 + X(I+1)
  140   CONTINUE
        IF (MODN .EQ. 1) SUM1 = SUM1 + X(N)
        Y(1) = SUM1 + SUM2
        IF (MODN .EQ. 0) Y(N) = SUM1 - SUM2
        CALL RFFTF(N, X, W)
        RFTF = 0.0
        DO 150 I=1,N
          RFTF = MAX(RFTF, ABS(X(I)-Y(I)))
          X(I) = XH(I)
  150   CONTINUE
        RFTF = RFTF/N
        IF (RFTF .LE. ERRMAX) THEN
          IF (KPRINT .GE. 3) WRITE (LUN, 9020)
        ELSE
          IPASS = 0
          IF (KPRINT .GE. 2) WRITE (LUN, 9030)
        END IF
        SIGN = 1.0
        DO 180 I=1,N
          SUM = 0.5D0*X(1)
          ARG = (I-1)*DT
          IF (NS2 .LT. 2) GO TO 170
          DO 160 K=2,NS2
            ARG1 = (K-1)*ARG
            SUM = SUM + X(2*K-2)*COS(ARG1) - X(2*K-1)*SIN(ARG1)
  160     CONTINUE
  170     IF (MODN .EQ. 0) SUM = SUM + 0.5D0*SIGN*X(N)
          Y(I) = SUM + SUM
          SIGN = -SIGN
  180   CONTINUE
        CALL RFFTB(N, X, W)
        RFTB = 0.0
        DO 190 I=1,N
          RFTB = MAX(RFTB, ABS(X(I)-Y(I)))
          X(I) = XH(I)
          Y(I) = XH(I)
  190   CONTINUE
        RFTB = RFTB/N
        IF (RFTB .LE. ERRMAX) THEN
          IF (KPRINT .GE. 3) WRITE (LUN, 9040)
        ELSE
          IPASS = 0
          IF (KPRINT .GE. 2) WRITE (LUN, 9050)
        END IF
C
        CALL RFFTB(N, Y, W)
        CALL RFFTF(N, Y, W)
        CF = 1.0/N
        RFTFB = 0.0
        DO 200 I=1,N
          RFTFB = MAX(RFTFB, ABS(CF*Y(I)-X(I)))
  200   CONTINUE
        IF (RFTFB .LE. ERRMAX) THEN
          IF (KPRINT .GE. 3) WRITE (LUN, 9060)
        ELSE
          IPASS = 0
          IF (KPRINT .GE. 2) WRITE (LUN, 9070)
        END IF
C
C       Test Subroutines SINTI and SINT
C
        DT = PI/N
        DO 210 I=1,NM1
          X(I) = XH(I)
  210   CONTINUE
        DO 230 I=1,NM1
          Y(I) = 0.0
          ARG1 = I*DT
          DO 220 K=1,NM1
            Y(I) = Y(I) + X(K)*SIN((K)*ARG1)
  220     CONTINUE
          Y(I) = Y(I) + Y(I)
  230   CONTINUE
        CALL SINTI(NM1, W)
        CALL SINT(NM1, X, W)
        CF = 0.5/N
        SINTT = 0.0
        DO 240 I=1,NM1
          SINTT = MAX(SINTT, ABS(X(I)-Y(I)))
          X(I) = XH(I)
          Y(I) = X(I)
  240   CONTINUE
        SINTT = CF*SINTT
        IF (SINTT .LE. ERRMAX) THEN
          IF (KPRINT .GE. 3) WRITE (LUN, 9080)
        ELSE
          IPASS = 0
          IF (KPRINT .GE. 2) WRITE (LUN, 9090)
        END IF
        CALL SINT(NM1, X, W)
        CALL SINT(NM1, X, W)
        SINTFB = 0.0
        DO 250 I=1,NM1
          SINTFB = MAX(SINTFB, ABS(CF*X(I)-Y(I)))
  250   CONTINUE
        IF (SINTFB .LE. ERRMAX) THEN
          IF (KPRINT .GE. 3) WRITE (LUN, 9100)
        ELSE
          IPASS = 0
          IF (KPRINT .GE. 2) WRITE (LUN, 9110)
        END IF
C
C       Test Subroutines COSTI and COST
C
        DO 260 I=1,NP1
          X(I) = XH(I)
  260   CONTINUE
        SIGN = 1.0
        DO 280 I=1,NP1
          Y(I) = 0.5*(X(1)+SIGN*X(N+1))
          ARG = (I-1)*DT
          DO 270 K=2,N
            Y(I) = Y(I) + X(K)*COS((K-1)*ARG)
  270     CONTINUE
          Y(I) = Y(I) + Y(I)
          SIGN = -SIGN
  280   CONTINUE
        CALL COSTI(NP1, W)
        CALL COST(NP1, X, W)
        COSTT = 0.0
        DO 290 I=1,NP1
          COSTT = MAX(COSTT, ABS(X(I)-Y(I)))
          X(I) = XH(I)
          Y(I) = XH(I)
  290   CONTINUE
        COSTT = CF*COSTT
        IF (COSTT .LE. ERRMAX) THEN
          IF (KPRINT .GE. 3) WRITE (LUN, 9120)
        ELSE
          IPASS = 0
          IF (KPRINT .GE. 2) WRITE (LUN, 9130)
        END IF
C
        CALL COST(NP1, X, W)
        CALL COST(NP1, X, W)
        COSTFB = 0.0
        DO 300 I=1,NP1
          COSTFB = MAX(COSTFB, ABS(CF*X(I)-Y(I)))
  300   CONTINUE
        IF (COSTFB .LE. ERRMAX) THEN
          IF (KPRINT .GE. 3) WRITE (LUN, 9140)
        ELSE
          IPASS = 0
          IF (KPRINT .GE. 2) WRITE (LUN, 9150)
        END IF
C
C       Test Subroutines SINQI, SINQF and SINQB
C
        CF = 0.25/N
        DO 310 I=1,N
          Y(I) = XH(I)
  310   CONTINUE
        DT = PI/(N+N)
        DO 330 I=1,N
          X(I) = 0.0
          ARG = I*DT
          DO 320 K=1,N
            X(I) = X(I) + Y(K)*SIN((K+K-1)*ARG)
  320     CONTINUE
          X(I) = 4.0*X(I)
  330   CONTINUE
        CALL SINQI(N, W)
        CALL SINQB(N, Y, W)
        SINQBT = 0.0
        DO 340 I=1,N
          SINQBT = MAX(SINQBT, ABS(Y(I)-X(I)))
          X(I) = XH(I)
  340   CONTINUE
        SINQBT = CF*SINQBT
        IF (SINQBT .LE. ERRMAX) THEN
          IF (KPRINT .GE. 3) WRITE (LUN, 9160)
        ELSE
          IPASS = 0
          IF (KPRINT .GE. 2) WRITE (LUN, 9170)
        END IF
C
        SIGN = 1.0
        DO 360 I=1,N
          ARG = (I+I-1)*DT
          Y(I) = 0.5*SIGN*X(N)
          DO 350 K=1,NM1
            Y(I) = Y(I) + X(K)*SIN((K)*ARG)
  350     CONTINUE
          Y(I) = Y(I) + Y(I)
          SIGN = -SIGN
  360   CONTINUE
        CALL SINQF(N, X, W)
        SINQFT = 0.0
        DO 370 I=1,N
          SINQFT = MAX(SINQFT, ABS(X(I)-Y(I)))
          Y(I) = XH(I)
          X(I) = XH(I)
  370   CONTINUE
        IF (SINQFT .LE. ERRMAX) THEN
          IF (KPRINT .GE. 3) WRITE (LUN, 9180)
        ELSE
          IPASS = 0
          IF (KPRINT .GE. 2) WRITE (LUN, 9190)
        END IF
C
        CALL SINQF(N, Y, W)
        CALL SINQB(N, Y, W)
        SINQFB = 0.0
        DO 380 I=1,N
          SINQFB = MAX(SINQFB, ABS(CF*Y(I)-X(I)))
  380   CONTINUE
        IF (SINQFB .LE. ERRMAX) THEN
          IF (KPRINT .GE. 3) WRITE (LUN, 9200)
        ELSE
          IPASS = 0
          IF (KPRINT .GE. 2) WRITE (LUN, 9210)
        END IF
C
C       Test Subroutines COSQI, COSQF and COSQB
C
        DO 390 I=1,N
          Y(I) = XH(I)
  390   CONTINUE
        DO 410 I=1,N
          X(I) = 0.0
          ARG = (I-1)*DT
          DO 400 K=1,N
            X(I) = X(I) + Y(K)*COS((K+K-1)*ARG)
  400     CONTINUE
          X(I) = 4.0*X(I)
  410   CONTINUE
        CALL COSQI(N, W)
        CALL COSQB(N, Y, W)
        COSQBT = 0.0
        DO 420 I=1,N
          COSQBT = MAX(COSQBT, ABS(X(I)-Y(I)))
          X(I) = XH(I)
  420   CONTINUE
        COSQBT = CF*COSQBT
        IF (COSQBT .LE. ERRMAX) THEN
          IF (KPRINT .GE. 3) WRITE (LUN, 9220)
        ELSE
          IPASS = 0
          IF (KPRINT .GE. 2) WRITE (LUN, 9230)
        END IF
C
        DO 440 I=1,N
          Y(I) = 0.5*X(1)
          ARG = (I+I-1)*DT
          DO 430 K=2,N
            Y(I) = Y(I) + X(K)*COS((K-1)*ARG)
  430     CONTINUE
          Y(I) = Y(I) + Y(I)
  440   CONTINUE
        CALL COSQF(N, X, W)
        COSQFT = 0.0
        DO 450 I=1,N
          COSQFT = MAX(COSQFT, ABS(Y(I)-X(I)))
          X(I) = XH(I)
          Y(I) = XH(I)
  450   CONTINUE
        COSQFT = CF*COSQFT
        IF (COSQFT .LE. ERRMAX) THEN
          IF (KPRINT .GE. 3) WRITE (LUN, 9240)
        ELSE
          IPASS = 0
          IF (KPRINT .GE. 2) WRITE (LUN, 9250)
        END IF
C
        CALL COSQB(N, X, W)
        CALL COSQF(N, X, W)
        COSQFB = 0.0
        DO 460 I=1,N
          COSQFB = MAX(COSQFB, ABS(CF*X(I)-Y(I)))
  460   CONTINUE
        IF (COSQFB .LE. ERRMAX) THEN
          IF (KPRINT .GE. 3) WRITE (LUN, 9260)
        ELSE
          IPASS = 0
          IF (KPRINT .GE. 2) WRITE (LUN, 9270)
        END IF
C
C       Test Subroutines EZFFTI, EZFFTF and EZFFTB
C
        CALL EZFFTI(N, W)
        DO 470 I=1,N
          X(I) = XH(I)
  470   CONTINUE
        TPI = 2.0*PI
        DT = TPI/N
        NS2 = (N+1)/2
        CF = 2.0/N
        NS2M = NS2 - 1
        IF (NS2M .LE. 0) GO TO 500
        DO 490 K=1,NS2M
          SUM1 = 0.0D0
          SUM2 = 0.0D0
          ARG = K*DT
          DO 480 I=1,N
            ARG1 = (I-1)*ARG
            SUM1 = SUM1 + X(I)*COS(ARG1)
            SUM2 = SUM2 + X(I)*SIN(ARG1)
  480     CONTINUE
          A(K) = CF*SUM1
          B(K) = CF*SUM2
  490   CONTINUE
  500   NM1 = N - 1
        SUM1 = 0.0D0
        SUM2 = 0.0D0
        DO 510 I=1,NM1,2
          SUM1 = SUM1 + X(I)
          SUM2 = SUM2 + X(I+1)
  510   CONTINUE
        IF (MODN .EQ. 1) SUM1 = SUM1 + X(N)
        AZERO = 0.5*CF*(SUM1+SUM2)
        IF (MODN .EQ. 0) A(NS2) = 0.5*CF*(SUM1-SUM2)
        CALL EZFFTF(N, X, AZEROH, AH, BH, W)
        DEZF1 = ABS(AZEROH-AZERO)
        IF (MODN .EQ. 0) DEZF1 = MAX(DEZF1, ABS(A(NS2)-AH(NS2)))
        IF (NS2M .LE. 0) GO TO 530
        DO 520 I=1,NS2M
          DEZF1 = MAX(DEZF1, ABS(AH(I)-A(I)), ABS(BH(I)-B(I)))
  520   CONTINUE
        IF (DEZF1 .LE. ERRMAX) THEN
          IF (KPRINT .GE. 3) WRITE (LUN, 9280)
        ELSE
          IPASS = 0
          IF (KPRINT .GE. 2) WRITE (LUN, 9290)
        END IF
C
  530   NS2 = N/2
        IF (MODN .EQ. 0) B(NS2) = 0.0
        DO 550 I=1,N
          SUM = AZERO
          ARG1 = (I-1)*DT
          DO 540 K=1,NS2
            ARG2 = K*ARG1
            SUM = SUM + A(K)*COS(ARG2) + B(K)*SIN(ARG2)
  540     CONTINUE
          X(I) = SUM
  550   CONTINUE
        CALL EZFFTB(N, Y, AZERO, A, B, W)
        DEZB1 = 0.0
        DO 560 I=1,N
          DEZB1 = MAX(DEZB1, ABS(X(I)-Y(I)))
          X(I) = XH(I)
  560   CONTINUE
        IF (DEZB1 .LE. ERRMAX) THEN
          IF (KPRINT .GE. 3) WRITE (LUN, 9300)
        ELSE
          IPASS = 0
          IF (KPRINT .GE. 2) WRITE (LUN, 9310)
        END IF
C
        CALL EZFFTF(N, X, AZERO, A, B, W)
        CALL EZFFTB(N, Y, AZERO, A, B, W)
        DEZFB = 0.0
        DO 570 I=1,N
          DEZFB = MAX(DEZFB, ABS(X(I)-Y(I)))
  570   CONTINUE
        IF (DEZFB .LE. ERRMAX) THEN
          IF (KPRINT .GE. 3) WRITE (LUN, 9320)
        ELSE
          IPASS = 0
          IF (KPRINT .GE. 2) WRITE (LUN, 9330)
        END IF
C
C       Test Subroutines CFFTI, CFFTF and CFFTB
C
        DO 580 I=1,N
          CX(I) = CMPLX(COS(SQRT2*I), SIN(SQRT2*(I*I)))
  580   CONTINUE
        DT = (PI+PI)/N
        DO 600 I=1,N
          ARG1 = -(I-1)*DT
          CY(I) = (0.0,0.0)
          DO 590 K=1,N
            ARG2 = (K-1)*ARG1
            CY(I) = CY(I) + CMPLX(COS(ARG2),SIN(ARG2))*CX(K)
  590     CONTINUE
  600   CONTINUE
        CALL CFFTI(N, W)
        CALL CFFTF(N, CX, W)
        DCFFTF = 0.0
        DO 610 I=1,N
          DCFFTF = MAX(DCFFTF, CABS(CX(I)-CY(I)))
          CX(I) = CX(I)/N
  610   CONTINUE
        DCFFTF = DCFFTF/N
        IF (DCFFTF .LE. ERRMAX) THEN
          IF (KPRINT .GE. 3) WRITE (LUN, 9340)
        ELSE
          IPASS = 0
          IF (KPRINT .GE. 2) WRITE (LUN, 9350)
        END IF
C
        DO 630 I=1,N
          ARG1 = (I-1)*DT
          CY(I) = (0.0,0.0)
          DO 620 K=1,N
            ARG2 = (K-1)*ARG1
            CY(I) = CY(I) + CMPLX(COS(ARG2),SIN(ARG2))*CX(K)
  620     CONTINUE
  630   CONTINUE
        CALL CFFTB(N, CX, W)
        DCFFTB = 0.0
        DO 640 I=1,N
          DCFFTB = MAX(DCFFTB, CABS(CX(I)-CY(I)))
          CX(I) = CY(I)
  640   CONTINUE
        IF (DCFFTB .LE. ERRMAX) THEN
          IF (KPRINT .GE. 3) WRITE (LUN, 9360)
        ELSE
          IPASS = 0
          IF (KPRINT .GE. 2) WRITE (LUN, 9370)
        END IF
C
        CF = 1.0/N
        CALL CFFTF(N, CX, W)
        CALL CFFTB(N, CX, W)
        DCFB = 0.0
        DO 650 I=1,N
          DCFB = MAX(DCFB, CABS(CF*CX(I)-CY(I)))
  650   CONTINUE
        IF (DCFB .LE. ERRMAX) THEN
          IF (KPRINT .GE. 3) WRITE (LUN, 9380)
        ELSE
          IPASS = 0
          IF (KPRINT .GE. 2) WRITE (LUN, 9390)
        END IF
        IF (KPRINT .GE. 3) THEN
          WRITE (LUN, 9400) N, RFTF, RFTB, RFTFB, SINTT, SINTFB,
     +      COSTT, COSTFB, SINQFT, SINQBT, SINQFB, COSQFT, COSQBT,
     +      COSQFB, DEZF1, DEZB1, DEZFB, DCFFTF, DCFFTB, DCFB
        END IF
  660 CONTINUE
      IF (KPRINT.GE.2 .AND. IPASS.EQ.1) WRITE (LUN, 9410)
      IF (KPRINT.GE.1 .AND. IPASS.EQ.0) WRITE (LUN, 9420)
      RETURN
C
 9000 FORMAT ('1' / ' FFT QUICK CHECK')
 9010 FORMAT (/ ' Test FFT routines with a sequence of length ', I3)
 9020 FORMAT (' Test of RFFTF PASSED')
 9030 FORMAT (' Test of RFFTF FAILED')
 9040 FORMAT (' Test of RFFTB PASSED')
 9050 FORMAT (' Test of RFFTB FAILED')
 9060 FORMAT (' Test of RFFTF and RFFTB PASSED')
 9070 FORMAT (' Test of RFFTF and RFFTB FAILED')
 9080 FORMAT (' First test of SINT PASSED')
 9090 FORMAT (' First test of SINT FAILED')
 9100 FORMAT (' Second test of SINT PASSED')
 9110 FORMAT (' Second test of SINT FAILED')
 9120 FORMAT (' First test of COST PASSED')
 9130 FORMAT (' First test of COST FAILED')
 9140 FORMAT (' Second test of COST PASSED')
 9150 FORMAT (' Second test of COST FAILED')
 9160 FORMAT (' Test of SINQB PASSED')
 9170 FORMAT (' Test of SINQB FAILED')
 9180 FORMAT (' Test of SINQF PASSED')
 9190 FORMAT (' Test of SINQF FAILED')
 9200 FORMAT (' Test of SINQF and SINQB PASSED')
 9210 FORMAT (' Test of SINQF and SINQB FAILED')
 9220 FORMAT (' Test of COSQB PASSED')
 9230 FORMAT (' Test of COSQB FAILED')
 9240 FORMAT (' Test of COSQF PASSED')
 9250 FORMAT (' Test of COSQF FAILED')
 9260 FORMAT (' Test of COSQF and COSQB PASSED')
 9270 FORMAT (' Test of COSQF and COSQB FAILED')
 9280 FORMAT (' Test of EZFFTF PASSED')
 9290 FORMAT (' Test of EZFFTF FAILED')
 9300 FORMAT (' Test of EZFFTB PASSED')
 9310 FORMAT (' Test of EZFFTB FAILED')
 9320 FORMAT (' Test of EZFFTF and EZFFTB PASSED')
 9330 FORMAT (' Test of EZFFTF and EZFFTB FAILED')
 9340 FORMAT (' Test of CFFTF PASSED')
 9350 FORMAT (' Test of CFFTF FAILED')
 9360 FORMAT (' Test of CFFTB PASSED')
 9370 FORMAT (' Test of CFFTB FAILED')
 9380 FORMAT (' Test of CFFTF and CFFTB PASSED')
 9390 FORMAT (' Test of CFFTF and CFFTB FAILED')
 9400 FORMAT ('0N', I5, '  RFFTF  ', E9.3, '  RFFTB  ', E9.3,
     +        '  RFFTFB ',E9.3 /
     +        7X, '  SINT   ', E9.3, '  SINTFB ', E9.3 /
     +        7X, '  COST   ', E9.3 , '  COSTFB '  , E9.3 /
     +        7X, '  SINQF  ', E9.3, '  SINQB  ', E9.3, '  SINQFB ',
     +        E9.3 /
     +        7X, '  COSQF  ', E9.3, '  COSQB  ', E9.3, '  COSQFB ',
     +        E9.3 /
     +        7X, '  DEZF1  ', E9.3, '  DEZB1  ', E9.3, '  DEZFB  ',
     +        E9.3 /
     +        7X, '  CFFTF  ', E9.3, '  CFFTB  ', E9.3, '  CFFTFB ',
     +        E9.3)
 9410 FORMAT (/ ' ***********FFT ROUTINES PASSED ALL TESTS************')
 9420 FORMAT (/ ' ***********FFT ROUTINES FAILED SOME TESTS***********')
      END
