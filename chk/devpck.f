*DECK DEVPCK
      SUBROUTINE DEVPCK (LOUT, KPRINT, X, Y, F, FX, FY, XE, YE, FE, DE,
     +   FE2, FAIL)
C***BEGIN PROLOGUE  DEVPCK
C***SUBSIDIARY
C***PURPOSE  Test usage of increment argument in DPCHFD and DPCHFE for
C            DPCHQ1.
C***LIBRARY   SLATEC (PCHIP)
C***TYPE      DOUBLE PRECISION (EVPCCK-S, DEVPCK-D)
C***KEYWORDS  PCHIP EVALUATOR QUICK CHECK
C***AUTHOR  Fritsch, F. N., (LLNL)
C***DESCRIPTION
C
C ---- CODE TO TEST USAGE OF INCREMENT ARGUMENT IN DPCHFD AND DPCHFE ---
C
C     EVALUATES A BICUBIC FUNCTION AND ITS FIRST PARTIAL DERIVATIVES
C     ON A 4X6 MESH CONTAINED IN A 10X10 ARRAY.
C
C     INTERPOLATION OF THESE DATA ALONG MESH LINES IN EITHER DIMENSION
C     SHOULD AGREE WITH CORRECT FUNCTION WITHIN ROUNDOFF ERROR.
C
C     ARRAYS ARE ARGUMENTS ONLY TO ALLOW SHARING STORAGE WITH OTHER
C     TEST ROUTINES.
C
C     NOTE:  RUN WITH KPRINT=4 FOR FULL GORY DETAILS (10 PAGES WORTH).
C
C
C     FORTRAN INTRINSICS USED:  ABS.
C     FORTRAN LIBRARY ROUTINES USED:  (WRITE).
C     SLATEC LIBRARY ROUTINES USED:  DPCHFD, DPCHFE, D1MACH.
C
C***ROUTINES CALLED  D1MACH, DPCHFD, DPCHFE
C***REVISION HISTORY  (YYMMDD)
C   820601  DATE WRITTEN
C   820714  CONVERTED TO QUICK CHECK FOR SLATEC LIBRARY.
C   820715  1. CORRECTED SOME FORMATS.
C           2. ADDED CALL TO D1MACH TO SET MACHEP.
C   890406  1. Modified to make sure final elements of X and XE
C             agree, to avoid possible failure due to roundoff
C             error.
C           2. Added printout of TOL in case of failure.
C           3. Removed unnecessary IMPLICIT declaration.
C           4. Corrected a few S.P. constants to D.P.
C           5. Minor cosmetic changes.
C   890706  Cosmetic changes to prologue.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891004  Cosmetic changes to prologue.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  Revised prologue and improved some output formats.  (FNF)
C   900316  Additional minor cosmetic changes.  (FNF)
C   900321  Made miscellaneous cosmetic changes.  (FNF)
C   901130  Made many changes to output:  (FNF)
C           1. Reduced amount of output for KPRINT=3.  (Now need to
C              use KPRINT=4 for full output.)
C           2. Added 1P's to formats and revised some to reduce maximum
C              line length.
C   910708  Minor modifications in use of KPRINT.  (WRB)
C   930317  Improved output formats.  (FNF)
C***END PROLOGUE  DEVPCK
C
C  Declare arguments.
C
      INTEGER  LOUT, KPRINT
      LOGICAL  FAIL
      DOUBLE PRECISION
     *      X(10), Y(10), F(10,10), FX(10,10), FY(10,10),
     *      XE(51), YE(51), FE(51), DE(51), FE2(51)
C
C  DECLARATIONS.
C
      INTEGER  I, IER2, IERR, INC, J, K, NE, NERR, NMAX, NX, NY
      LOGICAL  FAILD, FAILE, FAILOC, SKIP
      DOUBLE PRECISION
     *      DERMAX, DERR, DTRUE, DX, FDIFF, FDIFMX, FERMAX, FERR,
     *      FTRUE, MACHEP, TOL, PDERMX, PDIFMX, PFERMX, ZERO
      DOUBLE PRECISION  D1MACH
C
C  DEFINE TEST FUNCTION AND DERIVATIVES.
C
      DOUBLE PRECISION  AX, AY, FCN, DFDX, DFDY
      FCN (AX,AY) =  AX*(AY*AY)*(AX*AX + 1.D0)
      DFDX(AX,AY) = (AY*AY)*(3.D0*AX*AX + 1.D0)
      DFDY(AX,AY) =   2.D0*AX*AY*(AX*AX + 1.D0)
C
      DATA  NMAX /10/,  NX /4/,  NY /6/
      DATA  NE /51/
      DATA  ZERO /0.D0/
C
C  INITIALIZE.
C
C***FIRST EXECUTABLE STATEMENT  DEVPCK
      MACHEP = D1MACH(4)
C       Following tolerance is looser than S.P. version to avoid
C       spurious failures on some systems.
      TOL = 25.D0*MACHEP
C
      FAIL = .FALSE.
C
C  SET UP 4-BY-6 MESH IN A 10-BY-10 ARRAY:
C     X =  0.25(0.25)1.   ;
C     Y = -0.75(0.5 )1.75 .
C
      DO 1  I = 1, NX-1
         X(I) = 0.25D0*I
    1 CONTINUE
      X(NX) = 1.D0
      DO 5  J = 1, NY
         Y(J) = 0.5D0*J - 1.25D0
         DO 4  I = 1, NX
             F(I,J) = FCN (X(I), Y(J))
            FX(I,J) = DFDX(X(I), Y(J))
            FY(I,J) = DFDY(X(I), Y(J))
    4    CONTINUE
    5 CONTINUE
C
C  SET UP EVALUATION POINTS:
C     XE =  0.(0.02)1. ;
C     YE = -2.(0.08)2. .
C
      DX = 1.D0/(NE-1)
      DO 8  K = 1, NE-1
         XE(K) = DX*(K-1)
         YE(K) = 4.D0*XE(K) - 2.D0
    8 CONTINUE
      XE(NE) = 1.D0
      YE(NE) = 2.D0
C
      IF (KPRINT .GE. 3)  WRITE (LOUT, 1000)
      IF (KPRINT .GE. 2)  WRITE (LOUT, 1001)
C
C  EVALUATE ON HORIZONTAL MESH LINES (Y FIXED, X RUNNING) ..............
C
      NERR = 0
      INC = 1
      SKIP = .FALSE.
      DO 20  J = 1, NY
C        --------------------------------------------------------------
         CALL DPCHFD (NX, X, F(1,J), FX(1,J), INC, SKIP, NE, XE, FE, DE,
     *               IERR)
C        --------------------------------------------------------------
         IF (KPRINT .GE. 3)
     *       WRITE (LOUT, 2000)  INC, 'J', J, 'Y', Y(J), IERR
         IF (IERR .LT. 0)  GO TO 15
         IF (KPRINT .GT. 3)  WRITE (LOUT, 2001)  'X'
C
C        DPCHFE SHOULD AGREE EXACTLY WITH DPCHFD.
C
C        -----------------------------------------------------------
         CALL DPCHFE (NX, X, F(1,J), FX(1,J), INC, SKIP, NE, XE, FE2,
     *               IER2)
C        -----------------------------------------------------------
C
         DO 10  K = 1, NE
            FTRUE =  FCN(XE(K), Y(J))
            FERR = FE(K) - FTRUE
            DTRUE = DFDX(XE(K), Y(J))
            DERR = DE(K) - DTRUE
            IF (KPRINT .GT. 3)
     *         WRITE (LOUT, 2002)  XE(K), FTRUE, FE(K), FERR,
     *                                    DTRUE, DE(K), DERR
            IF (K .EQ. 1)  THEN
C              INITIALIZE.
               FERMAX = ABS(FERR)
               PFERMX = XE(1)
               DERMAX = ABS(DERR)
               PDERMX = XE(1)
               FDIFMX = ABS(FE2(1) - FE(1))
               PDIFMX = XE(1)
            ELSE
C              SELECT.
               FERR = ABS(FERR)
               IF (FERR .GT. FERMAX)  THEN
                  FERMAX = FERR
                  PFERMX = XE(K)
               ENDIF
               DERR = ABS(DERR)
               IF (DERR .GT. DERMAX)  THEN
                  DERMAX = DERR
                  PDERMX = XE(K)
               ENDIF
               FDIFF = ABS(FE2(K) - FE(K))
               IF (FDIFF .GT. FDIFMX)  THEN
                  FDIFMX = FDIFF
                  PDIFMX = XE(K)
               ENDIF
            ENDIF
   10    CONTINUE
C
         FAILD = (FERMAX.GT.TOL) .OR. (DERMAX.GT.TOL)
         FAILE = FDIFMX .NE. ZERO
         FAILOC = FAILD .OR. FAILE .OR. (IERR.NE.13) .OR. (IER2.NE.IERR)
C
         IF (FAILOC .AND. (KPRINT.GE.2))
     *      WRITE (LOUT, 2003)  'J', J, 'Y', Y(J)
C
         IF ((KPRINT.GE.3) .OR. (FAILD.AND.(KPRINT.EQ.2)) )
     *      WRITE (LOUT, 2004)  FERMAX, PFERMX, DERMAX, PDERMX
         IF (FAILD .AND. (KPRINT.GE.2))  WRITE (LOUT, 2014)  TOL
C
         IF ((KPRINT.GE.3) .OR. (FAILE.AND.(KPRINT.EQ.2)) )
     *      WRITE (LOUT, 2005)  FDIFMX, PDIFMX
C
         IF ((IERR.NE.13) .AND. (KPRINT.GE.2))
     *      WRITE (LOUT, 2006)  'D', IERR, 13
C
         IF ((IER2.NE.IERR) .AND. (KPRINT.GE.2))
     *      WRITE (LOUT, 2006)  'E', IER2, IERR
         GO TO 19
C
   15    CONTINUE
         FAILOC = .TRUE.
         IF (KPRINT .GE. 2)  WRITE (LOUT, 3000) IERR
C
   19    CONTINUE
         IF (FAILOC)  NERR = NERR + 1
         FAIL = FAIL .OR. FAILOC
   20 CONTINUE
C
      IF (KPRINT .GE. 2)  THEN
         IF (NERR .GT. 0)  THEN
            WRITE (LOUT, 3001)  NERR, 'J'
         ELSE
            WRITE (LOUT, 4000)  'J'
         ENDIF
      ENDIF
C
C  EVALUATE ON VERTICAL MESH LINES (X FIXED, Y RUNNING) ................
C
      NERR = 0
      INC = NMAX
      SKIP = .FALSE.
      DO 40  I = 1, NX
C        --------------------------------------------------------------
         CALL DPCHFD (NY, Y, F(I,1), FY(I,1), INC, SKIP, NE, YE, FE, DE,
     *               IERR)
C        --------------------------------------------------------------
         IF (KPRINT .GE. 3)
     *       WRITE (LOUT, 2000)  INC, 'I', I, 'X', X(I), IERR
         IF (IERR .LT. 0)  GO TO 35
         IF (KPRINT .GT. 3)  WRITE (LOUT, 2001)  'Y'
C
C        DPCHFE SHOULD AGREE EXACTLY WITH DPCHFD.
C
C        -----------------------------------------------------------
         CALL DPCHFE (NY, Y, F(I,1), FY(I,1), INC, SKIP, NE, YE, FE2,
     *               IER2)
C        -----------------------------------------------------------
C
         DO 30  K = 1, NE
            FTRUE =  FCN(X(I), YE(K))
            FERR = FE(K) - FTRUE
            DTRUE = DFDY(X(I), YE(K))
            DERR = DE(K) - DTRUE
            IF (KPRINT .GT. 3)
     *         WRITE (LOUT, 2002)  YE(K), FTRUE, FE(K), FERR,
     *                                    DTRUE, DE(K), DERR
            IF (K .EQ. 1)  THEN
C              INITIALIZE.
               FERMAX = ABS(FERR)
               PFERMX = YE(1)
               DERMAX = ABS(DERR)
               PDERMX = YE(1)
               FDIFMX = ABS(FE2(1) - FE(1))
               PDIFMX = YE(1)
            ELSE
C              SELECT.
               FERR = ABS(FERR)
               IF (FERR .GT. FERMAX)  THEN
                  FERMAX = FERR
                  PFERMX = YE(K)
               ENDIF
               DERR = ABS(DERR)
               IF (DERR .GT. DERMAX)  THEN
                  DERMAX = DERR
                  PDERMX = YE(K)
               ENDIF
               FDIFF = ABS(FE2(K) - FE(K))
               IF (FDIFF .GT. FDIFMX)  THEN
                  FDIFMX = FDIFF
                  PDIFMX = YE(K)
               ENDIF
            ENDIF
   30    CONTINUE
C
         FAILD = (FERMAX.GT.TOL) .OR. (DERMAX.GT.TOL)
         FAILE = FDIFMX .NE. ZERO
         FAILOC = FAILD .OR. FAILE .OR. (IERR.NE.20) .OR. (IER2.NE.IERR)
C
         IF (FAILOC .AND. (KPRINT.GE.2))
     *      WRITE (LOUT, 2003)  'I', I, 'X', X(I)
C
         IF ((KPRINT.GE.3) .OR. (FAILD.AND.(KPRINT.EQ.2)) )
     *      WRITE (LOUT, 2004)  FERMAX, PFERMX, DERMAX, PDERMX
         IF (FAILD .AND. (KPRINT.GE.2))  WRITE (LOUT, 2014)  TOL
C
         IF ((KPRINT.GE.3) .OR. (FAILE.AND.(KPRINT.EQ.2)) )
     *      WRITE (LOUT, 2005)  FDIFMX, PDIFMX
C
         IF ((IERR.NE.20) .AND. (KPRINT.GE.2))
     *      WRITE (LOUT, 2006)  'D', IERR, 20
C
         IF ((IER2.NE.IERR) .AND. (KPRINT.GE.2))
     *      WRITE (LOUT, 2006)  'E', IER2, IERR
         GO TO 39
C
   35    CONTINUE
         FAILOC = .TRUE.
         IF (KPRINT .GE. 2)  WRITE (LOUT, 3000) IERR
C
   39    CONTINUE
         IF (FAILOC)  NERR = NERR + 1
         FAIL = FAIL .OR. FAILOC
   40 CONTINUE
C
      IF (KPRINT .GE. 2)  THEN
         IF (NERR .GT. 0)  THEN
            WRITE (LOUT, 3001)  NERR, 'I'
         ELSE
            WRITE (LOUT, 4000)  'I'
         ENDIF
      ENDIF
C
C  TERMINATE.
C
      RETURN
C
C  FORMATS.
C
 1000 FORMAT ('1'//10X,'TEST DPCHFE AND DPCHFD')
 1001 FORMAT (//10X,'DEVPCK RESULTS'/10X,'--------------')
 2000 FORMAT (//20X,'DPCHFD INCREMENT TEST -- INCFD = ',I2
     *        /15X,'ON ',A1,'-LINE ',I2,',  ',A1,' =',F8.4,
     *           '  --  IERR =',I3)
 2001 FORMAT ( /3X,A1,'E',10X,'F',8X,'FE',9X,'DIFF',
     *                    13X,'D',8X,'DE',9X,'DIFF')
 2002 FORMAT (F7.2,2(2X,2F10.5,1P,D15.5,0P))
 2003 FORMAT (/' ***** DPCHFD AND/OR DPCHFE FAILED ON ',A1,'-LINE ',I1,
     *                             ',  ',A1,' =',F8.4)
 2004 FORMAT (/19X,'  MAXIMUM ERROR IN FUNCTION =',1P,
     *                                   1P,D13.5,0P,' (AT',F6.2,'),'
     *        /33X,    'IN DERIVATIVE =',1P,D13.5,0P,' (AT',F6.2,').' )
 2005 FORMAT ( '  MAXIMUM DIFFERENCE BETWEEN DPCHFE AND DPCHFD =',
     *                                   1P,D13.5,0P,' (AT',F6.2,').' )
 2006 FORMAT (/'  DPCHF',A1,' RETURNED IERR = ',I2,' INSTEAD OF ',I2)
 2014 FORMAT ('  *** BOTH SHOULD BE .LE. TOL =',1P,D12.5,' ***')
 3000 FORMAT (//' ***** ERROR ***** DPCHFD RETURNED IERR =',I5//)
 3001 FORMAT (//' ***** ERROR ***** DPCHFD AND/OR DPCHFE FAILED ON',I2,
     *                                1X, A1,'-LINES.'//)
 4000 FORMAT (/' DPCHFD AND DPCHFE OK ON ',A1,'-LINES.')
C------------- LAST LINE OF DEVPCK FOLLOWS -----------------------------
      END
