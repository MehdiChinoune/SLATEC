*DECK SDASQC
      SUBROUTINE SDASQC (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  SDASQC
C***PURPOSE  Quick check for SLATEC routine SDASSL.
C***LIBRARY   SLATEC (DASSL)
C***CATEGORY  I1A2
C***TYPE      SINGLE PRECISION (SDASQC-S, DDASQC-D)
C***KEYWORDS  QUICK CHECK, SDASSL
C***AUTHOR  PETZOLD, LINDA R., (LLNL)
C             COMPUTING AND MATHEMATICS RESEARCH DIVISION
C             LAWRENCE LIVERMORE NATIONAL LABORATORY
C             L - 316, P.O. BOX 808,
C             LIVERMORE, CA.    94550
C***DESCRIPTION
C       Demonstration program for SDASSL.
C
C       SDASSL is used to solve two simple problems,
C       one with a full Jacobian, the other with a banded Jacobian,
C       with the 2 appropriate values of info(5) in each case.
C       If the errors are too large, or other difficulty occurs,
C       a warning message is printed.  All output is on unit LUN.
C
C       To run the demonstration problems with full print, call
C       SDASQC with KPRINT = 3.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  EDIT2, SDASSL, SDJAC1, SDJAC2, SDRES1, SDRES2
C***REVISION HISTORY  (YYMMDD)
C   851113  DATE WRITTEN
C   880608  Revised to meet new LDOC standards.
C   880615  Revised to meet new prologue standards.
C   891016  Converted to 4.0 format (FNF).
C   901001  Minor improvements to prologue.  (FNF)
C   901003  Added some error prints and improved output formats.  (FNF)
C   901009  Corrected GAMS classification code.  (FNF)
C   901009  Changed AMAX1 to MAX.  (FNF)
C   901030  Made all declarations explicit; added 1P's to formats. (FNF)
C***END PROLOGUE  SDASQC
C
      INTEGER  LUN, KPRINT, IPASS
C
      EXTERNAL  EDIT2, SDASSL, SDJAC1, SDJAC2, SDRES1, SDRES2
C
      INTEGER  I, IDID, INFO(15), IOUT, IPAR(1), IRES, IWORK(45),
     *   J190, J290, LIW, LRW, ML, MU, NEQ, NERR, NFE, NJE, NOUT,
     *   NQU, NST
      REAL  ATOL, DELTA(25), DTOUT, ER, ER1, ER2, ERM, ERO, HU, RPAR(1),
     *   RTOL, RWORK(550), T, TOUT, TOUT1, Y(25), YPRIME(25), YT1, YT2
C
      DATA TOUT1/1.0E0/, DTOUT/1.0E0/
C
C***FIRST EXECUTABLE STATEMENT  SDASQC
      IPASS = 1
      NERR = 0
      RTOL = 0.0E0
      ATOL = 1.0E-3
      LRW = 550
      LIW = 45
C
C FIRST PROBLEM
C
      NEQ = 2
      NOUT = 10
      IF (KPRINT .GE. 2)  WRITE (LUN,110) NEQ,RTOL,ATOL
 110  FORMAT('1'/1X,' DEMONSTRATION PROGRAM FOR SDASSL',///
     1  1X,' PROBLEM 1..   LINEAR DIFFERENTIAL/ALGEBRAIC SYSTEM..',/
     2  1X,' X1DOT + 10.0*X1 = 0,  X1 + X2 = 1',/
     3  1X,' X1(0) = 1.0, X2(0) = 0.0',/
     4  1X,' NEQ =',I2/
     5  1X,' RTOL =',1P,E10.1,'   ATOL =',E10.1)
C
      DO 190 J190 = 1,2
      DO 115 I = 1,15
 115    INFO(I) = 0
      IF(J190 .EQ. 2) INFO(5) = 1
C
      IF (KPRINT .GT. 2)  WRITE (LUN,120) INFO(5)
 120  FORMAT(////1X,' INFO(5) =',I3//
     1  6X,'T',15X,'X1',14X,'X2',12X,'NQ',6X,'H',12X/)
C
      T = 0.0E0
      Y(1) = 1.0E0
      Y(2) = 0.0E0
      YPRIME(1) = -10.0E0
      YPRIME(2) =  10.0E0
      TOUT = TOUT1
      ERO = 0.0E0
      DO 170 IOUT = 1,NOUT
        CALL SDASSL(SDRES1,NEQ,T,Y,YPRIME,TOUT,INFO,RTOL,ATOL,IDID,
     1     RWORK,LRW,IWORK,LIW,RPAR,IPAR,SDJAC1)
        HU = RWORK(7)
        NQU = IWORK(8)
        IF (KPRINT .GT. 2) THEN
          WRITE (LUN,140) T,Y(1),Y(2),NQU,HU
        ENDIF
 140    FORMAT(1X,1P,E15.5,E16.5,E16.5,I6,E14.3)
C
        IF (IDID .LT. 0) GO TO 175
        YT1 = EXP(-10.0E0*T)
        YT2 = 1.0E0 - YT1
        ER1 = ABS(YT1 - Y(1))
        ER2 = ABS(YT2 - Y(2))
        ER = MAX(ER1,ER2)/ATOL
        ERO = MAX(ERO,ER)
        IF (ER .GT. 1000.0E0)  THEN
          IF (KPRINT .GE. 2)  WRITE (LUN,150) T
 150      FORMAT(//' WARNING.. ERROR EXCEEDS 1000 * TOLERANCE',
     +             '  WHEN  T =',1P,E13.5//)
C
          NERR = NERR + 1
        ENDIF
 170    TOUT = TOUT + DTOUT
 175  CONTINUE
      IF (IDID .LT. 0)  THEN
        IF (KPRINT .GE. 2)  WRITE (LUN, 176)  IDID, T
 176    FORMAT (//'TROUBLE..  SDASSL RETURNED  IDID =',I4,
     +            '  WHEN  T =',1P,E13.5)
C
        NERR = NERR + 1
      ENDIF
      NST = IWORK(11)
      NFE = IWORK(12)
      NJE = IWORK(13)
      IF (KPRINT .GT. 2)  WRITE (LUN,180) NST,NFE,NJE,ERO
 180  FORMAT(//1X,' FINAL STATISTICS FOR THIS RUN..',/
     1  1X,' NUMBER OF STEPS =',I5/
     2  1X,' NUMBER OF F-S   =',I5/
     4  1X,' NUMBER OF J-S   =',I5/
     5  1X,' ERROR OVERRUN =',1P,E10.2)
C
 190  CONTINUE
C
C SECOND PROBLEM
C
      NEQ = 25
      ML = 5
      MU = 0
      IWORK(1) = ML
      IWORK(2) = MU
      NOUT = 5
      IF (KPRINT .GE. 2)  WRITE (LUN,210) NEQ,ML,MU,RTOL,ATOL
 210  FORMAT('1'/1X,' DEMONSTRATION PROGRAM FOR SDASSL',///
     1  1X,' PROBLEM 2.. YDOT = A * Y , WHERE ',
     2  ' A IS A BANDED LOWER TRIANGULAR MATRIX',/
     2  1X,'  DERIVED FROM 2-D ADVECTION PDE',/
     3  1X,' NEQ =',I3,'   ML =',I2,'   MU =',I2/
     4  1X,' RTOL =',1P,E10.1,'   ATOL =',E10.1)
C
      DO 290 J290 = 1,2
      DO 215 I = 1,15
 215    INFO(I) = 0
      INFO(6) = 1
      IF(J290 .EQ. 2) INFO(5) = 1
C
      IF (KPRINT .GT. 2)  WRITE (LUN,220) INFO(5)
 220  FORMAT(////1X,' INFO(5) =',I3//
     1  6X,'T',14X,'MAX.ERR.',5X,'NQ',6X,'H'/)
C
      T = 0.0E0
      DO 230 I = 2,NEQ
 230    Y(I) = 0.0E0
      Y(1) = 1.0E0
      DO 235 I = 1,NEQ
 235    DELTA(I) = 0.0E0
C        Following is to initialize YPRIME.
      CALL SDRES2(T,Y,DELTA,YPRIME,IRES,RPAR,IPAR)
      TOUT = 0.01E0
      ERO = 0.0E0
      DO 270 IOUT = 1,NOUT
        CALL SDASSL(SDRES2,NEQ,T,Y,YPRIME,TOUT,INFO,RTOL,ATOL,IDID,
     1     RWORK,LRW,IWORK,LIW,RPAR,IPAR,SDJAC2)
        CALL EDIT2(Y,T,ERM)
        HU = RWORK(7)
        NQU = IWORK(8)
        IF (KPRINT .GT. 2)  WRITE (LUN,240) T,ERM,NQU,HU
 240    FORMAT(1X,1P,E15.5,E14.3,I6,E14.3)
C
        IF (IDID .LT. 0) GO TO 275
        ER = ERM/ATOL
        ERO = MAX(ERO,ER)
        IF (ER .GT. 1000.0E0) THEN
          IF (KPRINT .GE. 2)  WRITE (LUN,150) T
          NERR = NERR + 1
        ENDIF
 270    TOUT = TOUT*10.0E0
 275  CONTINUE
      IF (IDID .LT. 0)  THEN
        IF (KPRINT .GE. 2)  WRITE (LUN, 176)  IDID, T
        NERR = NERR + 1
      ENDIF
      NST = IWORK(11)
      NFE = IWORK(12)
      NJE = IWORK(13)
      IF (KPRINT .GT. 2)  WRITE (LUN,180) NST,NFE,NJE,ERO
C
 290  CONTINUE
      IF (KPRINT .GE. 2)  WRITE (LUN,300) NERR
C
 300  FORMAT(////1X,' NUMBER OF ERRORS ENCOUNTERED =',I3)
C
      IF (NERR .GT. 0) THEN
        IPASS = 0
      ENDIF
      IF ((IPASS .EQ. 1) .AND. (KPRINT .GT. 1)) WRITE (LUN,700)
      IF ((IPASS .EQ. 0) .AND. (KPRINT .NE. 0)) WRITE (LUN,800)
 700  FORMAT (/,' ----------SDASSL PASSED ALL TESTS----------')
 800  FORMAT (/,' **********SDASSL FAILED SOME TESTS*********')
      RETURN
      END
