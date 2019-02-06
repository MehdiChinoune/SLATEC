*DECK EVCHCK
      SUBROUTINE EVCHCK (LOUT, KPRINT, NPTS, XEV, FEV, DEV, FEV2, FAIL)
C***BEGIN PROLOGUE  EVCHCK
C***SUBSIDIARY
C***PURPOSE  Test evaluation accuracy of CHFDV and CHFEV for PCHQK1.
C***LIBRARY   SLATEC (PCHIP)
C***TYPE      SINGLE PRECISION (EVCHCK-S, DEVCHK-D)
C***KEYWORDS  PCHIP EVALUATOR QUICK CHECK
C***AUTHOR  Fritsch, F. N., (LLNL)
C***DESCRIPTION
C
C -------- CODE TO TEST EVALUATION ACCURACY OF CHFDV AND CHFEV --------
C
C     USING FUNCTION AND DERIVATIVE VALUES FROM A CUBIC (COMPUTED IN
C     DOUBLE PRECISION) AT NINT DIFFERENT (X1,X2) PAIRS:
C     1. CHECKS THAT CHFDV AND CHFEV BOTH REPRODUCE ENDPOINT VALUES.
C     2. EVALUATES AT NPTS POINTS, 10 OF WHICH ARE OUTSIDE THE INTERVAL
C        AND:
C        A. CHECKS ACCURACY OF CHFDV FUNCTION AND DERIVATIVE VALUES
C           AGAINST EXACT VALUES.
C        B. CHECKS THAT RETURNED VALUES OF NEXT SUM TO 10.
C        C. CHECKS THAT FUNCTION VALUES FROM CHFEV AGREE WITH THOSE
C           FROM CHFDV.
C
C
C     FORTRAN INTRINSICS USED:  ABS, MAX, MIN.
C     FORTRAN LIBRARY ROUTINES USED:  SQRT, (READ), (WRITE).
C     SLATEC LIBRARY ROUTINES USED:  CHFDV, CHFEV, R1MACH, RAND.
C     OTHER ROUTINES USED:  FDTRUE.
C
C***ROUTINES CALLED  CHFDV, CHFEV, FDTRUE, R1MACH, RAND
C***REVISION HISTORY  (YYMMDD)
C   820601  DATE WRITTEN
C   820624  CONVERTED TO QUICK CHECK FOR SLATEC LIBRARY.
C   820630  1. MODIFIED DEFINITIONS OF RELATIVE ERROR AND TEST
C             TOLERANCES.
C           2. VARIOUS IMPROVEMENTS TO OUTPUT FORMATS.
C   820716  1. SET MACHEP VIA A CALL TO R1MACH.
C           2. CHANGED FROM FORTLIB'S RANF TO SLATEC'S RAND.
C   890629  1. Appended E0 to real constants to reduce S.P./D.P.
C             differences.
C           2. Other minor cosmetic changes.
C   890831  Modified array declarations.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891004  Cosmetic changes to prologue.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  Revised prologue and improved some output formats.  (FNF)
C           Also moved formats to end to be consistent with other PCHIP
C           quick checks.
C   900316  Additional minor cosmetic changes.  (FNF)
C   900321  Made miscellaneous cosmetic changes.  (FNF)
C   901130  Added 1P's to formats and revised some to reduce maximum
C           line length.  (FNF)
C   910708  Minor modifications in use of KPRINT.  (WRB)
C   910801  Added EXTERNAL statement for RAND due to problem on IBM
C           RS 6000.  (WRB)
C***END PROLOGUE  EVCHCK
C
C  Declare arguments.
C
      INTEGER  LOUT, KPRINT, NPTS
      REAL  XEV(*), FEV(*), DEV(*), FEV2(*)
      LOGICAL  FAIL
C
C  DECLARATIONS.
C
      INTEGER  I, IERR, IINT, NEXT(2), NEXT2(2), NINT
      REAL  AED, AED2, AEDMAX, AEDMIN, AEF, AEF2, AEFMAX, AEFMIN,
     *      CHECK(2), CHECKF(2), CHECKD(2), D1, D2, DERMAX, DTRUE, DX,
     *      EPS1, EPS2, F1, F2, FACT, FERMAX, FLOORD, FLOORF, FOUR,
     *      FTRUE, LEFT(3), MACHEP,
     *      ONE, RED, RED2, REDMAX, REDMIN, REF, REF2, REFMAX,
     *      REFMIN, RIGHT(3), SMALL, TEN, TOL1, TOL2,
     *      X1, X2, XADMAX, XADMIN, XAFMAX, XAFMIN, XRDMAX,
     *      XRDMIN, XRFMAX, XRFMIN, ZERO
      LOGICAL  FAILOC, FAILNX
C
      REAL  R1MACH
C       The following should stay REAL (no D.P. equivalent).
      REAL  RAND
      EXTERNAL  RAND
C
C  DEFINE RELATIVE ERROR WITH FLOOR.
C
      REAL  RERR, ERR, VALUE, FLOOR
      RERR(ERR,VALUE,FLOOR) = ERR / MAX(ABS(VALUE), FLOOR)
C
C  INITIALIZE.
C
      DATA  ZERO /0.E0/, ONE /1.E0/, FOUR /4.E0/, TEN /10.E0/
      DATA  SMALL  /1.0E-10/
      DATA  NINT /3/
      DATA   LEFT /-1.5E0, 2.0E-10, 1.0E0 /
      DATA  RIGHT / 2.5E0, 3.0E-10, 1.0E+8/
C
C***FIRST EXECUTABLE STATEMENT  EVCHCK
      MACHEP = R1MACH(4)
      EPS1 = FOUR*MACHEP
      EPS2 = TEN*MACHEP
C
      FAIL = .FALSE.
C
      IF (KPRINT .GE. 2)  WRITE (LOUT, 3000)
C
C  CYCLE OVER INTERVALS.
C
      DO 90  IINT = 1, NINT
      X1 =  LEFT(IINT)
      X2 = RIGHT(IINT)
C
      FACT = MAX(SQRT(X2-X1), ONE)
      TOL1 = EPS1 * FACT
      TOL2 = EPS2 * FACT
C
C  COMPUTE AND PRINT ENDPOINT VALUES.
C
      CALL FDTRUE (X1, F1, D1)
      CALL FDTRUE (X2, F2, D2)
C
      IF (KPRINT .GE. 3)  THEN
         IF (IINT .EQ. 1)  WRITE (LOUT, 2000)
         WRITE (LOUT, '(/)')
         WRITE (LOUT, 2001)  'X1', X1, 'X2', X2
         WRITE (LOUT, 2001)  'F1', F1, 'F2', F2
         WRITE (LOUT, 2001)  'D1', D1, 'D2', D2
      ENDIF
C
      IF (KPRINT .GE. 2)  WRITE (LOUT, 3001)  X1, X2
C
C  COMPUTE FLOORS FOR RELATIVE ERRORS.
C
      FLOORF = MAX( MIN(ABS(F1),ABS(F2)), SMALL)
      FLOORD = MAX( MIN(ABS(D1),ABS(D2)), SMALL)
C
C  CHECK REPRODUCTION OF ENDPOINT VALUES.
C
      XEV(1) = X1
      XEV(2) = X2
C     -----------------------------------------------------------
      CALL CHFDV (X1, X2, F1, F2, D1, D2, 2, XEV, CHECKF, CHECKD,
     *            NEXT, IERR)
C     -----------------------------------------------------------
      AEF  = CHECKF(1)-F1
      REF  = RERR(AEF , F1, FLOORF)
      AEF2 = CHECKF(2)-F2
      REF2 = RERR(AEF2, F2, FLOORF)
      AED  = CHECKD(1)-D1
      RED  = RERR(AED , D1, FLOORD)
      AED2 = CHECKD(2)-D2
      RED2 = RERR(AED2, D2, FLOORD)
C
      FAILOC = MAX(ABS(REF),ABS(REF2),ABS(RED),ABS(RED2)) .GT. TOL1
      FAIL = FAIL .OR. FAILOC
C
      IF (KPRINT .GE. 3)  THEN
         WRITE (LOUT, 2002)  NEXT, AEF, AEF2, AED, AED2
         WRITE (LOUT, 2003)  REF, REF2, RED, RED2
      ENDIF
C
      IF (FAILOC .AND. (KPRINT.GE.2))  WRITE (LOUT, 3002)
C
C  CHFEV SHOULD AGREE EXACTLY WITH CHFDV.
C                     -------
C     --------------------------------------------------------------
      CALL CHFEV (X1, X2, F1, F2, D1, D2, 2, XEV, CHECK, NEXT, IERR)
C     --------------------------------------------------------------
      FAILOC = (CHECK(1).NE.CHECKF(1)) .OR. (CHECK(2).NE.CHECKF(2))
      FAIL = FAIL .OR. FAILOC
C
      IF (FAILOC .AND. (KPRINT.GE.2))  WRITE (LOUT, 3003)
C
C  EVALUATE AT NPTS 'UNIFORMLY RANDOM' POINTS IN (X1,X2).
C     THIS VERSION EXTENDS EVALUATION DOMAIN BY ADDING 4 SUBINTERVALS
C     TO LEFT AND 6 TO RIGHT OF [X1,X2].
C
      DX = (X2-X1)/(NPTS-10)
      DO 20  I = 1, NPTS
         XEV(I) = (X1 + (I-5)*DX) + DX*RAND(ZERO)
   20 CONTINUE
C     --------------------------------------------------------
      CALL CHFDV (X1, X2, F1, F2, D1, D2, NPTS, XEV, FEV, DEV,
     *            NEXT, IERR)
C     --------------------------------------------------------
      IF (IERR .NE. 0)  THEN
         FAILOC = .TRUE.
         IF (KPRINT .GE. 2)  WRITE (LOUT, 4003)  IERR
      ELSE
C
C     CUMULATE LARGEST AND SMALLEST ERRORS FOR SUMMARY.
C
      DO 30  I = 1, NPTS
         CALL FDTRUE (XEV(I), FTRUE, DTRUE)
         AEF = FEV(I) - FTRUE
         REF = RERR(AEF, FTRUE, FLOORF)
         AED = DEV(I) - DTRUE
         RED = RERR(AED, DTRUE, FLOORD)
C
         IF (I .EQ. 1)  THEN
C            INITIALIZE.
            AEFMIN = AEF
            AEFMAX = AEF
            AEDMIN = AED
            AEDMAX = AED
            REFMIN = REF
            REFMAX = REF
            REDMIN = RED
            REDMAX = RED
            XAFMIN = XEV(1)
            XAFMAX = XEV(1)
            XADMIN = XEV(1)
            XADMAX = XEV(1)
            XRFMIN = XEV(1)
            XRFMAX = XEV(1)
            XRDMIN = XEV(1)
            XRDMAX = XEV(1)
         ELSE
C            SELECT.
            IF (AEF .LT. AEFMIN)  THEN
               AEFMIN = AEF
               XAFMIN = XEV(I)
            ELSE IF (AEF .GT. AEFMAX)  THEN
               AEFMAX = AEF
               XAFMAX = XEV(I)
            ENDIF
            IF (AED .LT. AEDMIN)  THEN
               AEDMIN = AED
               XADMIN = XEV(I)
            ELSE IF (AED .GT. AEDMAX)  THEN
               AEDMAX = AED
               XADMAX = XEV(I)
            ENDIF
            IF (REF .LT. REFMIN)  THEN
               REFMIN = REF
               XRFMIN = XEV(I)
            ELSE IF (REF .GT. REFMAX)  THEN
               REFMAX = REF
               XRFMAX = XEV(I)
            ENDIF
            IF (RED .LT. REDMIN)  THEN
               REDMIN = RED
               XRDMIN = XEV(I)
            ELSE IF (RED .GT. REDMAX)  THEN
               REDMAX = RED
               XRDMAX = XEV(I)
            ENDIF
         ENDIF
   30    CONTINUE
C
         FERMAX = MAX (ABS(REFMAX), ABS(REFMIN))
         DERMAX = MAX (ABS(REDMAX), ABS(REDMIN))
C
         FAILNX = (NEXT(1) + NEXT(2)) .NE. 10
         FAILOC = FAILNX .OR. (MAX(FERMAX, DERMAX) .GT. TOL2)
      ENDIF
      FAIL = FAIL .OR. FAILOC
C
C  PRINT SUMMARY.
C
      IF (KPRINT .GE. 3)  THEN
         WRITE (LOUT, 2004)  NPTS-10, NEXT
C
         WRITE (LOUT, 2005)  'MIN', AEFMIN, REFMIN, AEDMIN, REDMIN
         WRITE (LOUT, 2006) XAFMIN, XRFMIN, XADMIN, XRDMIN
         WRITE (LOUT, 2005)  'MAX', AEFMAX, REFMAX, AEDMAX, REDMAX
         WRITE (LOUT, 2006) XAFMAX, XRFMAX, XADMAX, XRDMAX
      ENDIF
C
      IF (KPRINT .GE. 2)  THEN
         IF (FAILOC) THEN
            IF (FERMAX .GT. TOL2)  WRITE (LOUT, 3006) 'F', FERMAX, TOL2
            IF (DERMAX .GT. TOL2)  WRITE (LOUT, 3006) 'D', DERMAX, TOL2
            IF (FAILNX)  WRITE (LOUT, 4006)  NEXT
         ELSE
            WRITE (LOUT, 5006)
         ENDIF
      ENDIF
C
C  CHECK THAT CHFEV AGREES WITH CHFDV.
C
C     -----------------------------------------------------------------
      CALL CHFEV (X1, X2, F1, F2, D1, D2, NPTS, XEV, FEV2, NEXT2, IERR)
C     -----------------------------------------------------------------
      IF (IERR .NE. 0)  THEN
         FAILOC = .TRUE.
         IF (KPRINT .GE. 2)  WRITE (LOUT, 3007)  IERR
      ELSE
         AEFMAX = ABS(FEV2(1) - FEV(1))
         XAFMAX = XEV(1)
         DO 40  I = 2, NPTS
            AEF = ABS(FEV2(I) - FEV(I))
            IF (AEF .GT. AEFMAX)  THEN
               AEFMAX = AEF
               XAFMAX = XEV(I)
            ENDIF
   40    CONTINUE
         FAILNX = (NEXT2(1).NE.NEXT(1)) .OR. (NEXT2(2).NE.NEXT(2))
         FAILOC = FAILNX .OR. (AEFMAX.NE.ZERO)
         IF (KPRINT .GE. 2)  THEN
            IF (FAILOC)  THEN
               WRITE (LOUT, 3008)
               IF (AEFMAX.NE.ZERO)  WRITE (LOUT, 3009)  AEFMAX, XAFMAX
               IF (FAILNX)  WRITE (LOUT, 4009)  NEXT2, NEXT
            ELSE
               WRITE (LOUT, 5009)
            ENDIF
         ENDIF
      ENDIF
C
      FAIL = FAIL .OR. FAILOC
C
C  GO BACK FOR ANOTHER INTERVAL.
C
   90 CONTINUE
C
      RETURN
C
C  FORMATS.
C
 2000 FORMAT (/10X,'CHFDV ACCURACY TEST')
 2001 FORMAT (10X,A2,' =',1P,E18.10,5X,A2,' =',E18.10)
 2002 FORMAT (/' ERRORS AT ENDPOINTS:',40X,'(NEXT =',2I3,')'
     *        // 1P,4X,'F1:',E13.5,4X,'F2:',E13.5,
     *              4X,'D1:',E13.5,4X,'D2:',E13.5)
 2003 FORMAT (1P,4(7X,E13.5))
 2004 FORMAT (/' ERRORS AT ',I5,' INTERIOR POINTS + 10 OUTSIDE:',
     *                15X,'(NEXT =',2I3,')'
     *        //30X,'FUNCTION',17X,'DERIVATIVE'
     *         /15X,2(11X,'ABS',9X,'REL') )
 2005 FORMAT (/5X,A3,'IMUM ERROR:  ',1P,2E12.4,2X,2E12.4)
 2006 FORMAT ( 5X,'LOCATED AT X =  ',1P,2E12.4,2X,2E12.4)
 3000 FORMAT (//10X,'EVCHCK RESULTS'/10X,'--------------')
 3001 FORMAT (/10X,'INTERVAL = (',1P,E12.5,',',E12.5,' ):' )
 3002 FORMAT (/' ***** CHFDV FAILED TO REPRODUCE ENDPOINT VALUES.')
 3003 FORMAT (/' ***** CHFEV DOES NOT AGREE WITH CHFDV AT ENDPOINTS.')
 3006 FORMAT (/' ***** MAXIMUM RELATIVE ERROR IN ',A1,' =',1P,E12.5,','
     *        /        17X,'EXCEEDS TOLERANCE =',E12.5)
 3007 FORMAT (/' ***** ERROR ***** CHFEV RETURNED IERR =',I5)
 3008 FORMAT (/' ***** CHFEV DID NOT AGREE WITH CHFDV:')
 3009 FORMAT (7X,'MAXIMUM DIFFERENCE ',1P,E12.5,
     *                '; OCCURRED AT X =',E12.5)
 4003 FORMAT (/' ***** ERROR ***** CHFDV RETURNED IERR =',I5)
 4006 FORMAT (/' ***** REPORTED NEXT =',2I5,'   RATHER THAN    4    6')
 4009 FORMAT (7X,'REPORTED NEXT =',2I3,'   RATHER THAN ',2I3)
 5006 FORMAT (/' CHFDV RESULTS OK.')
 5009 FORMAT (/' CHFEV AGREES WITH CHFDV.')
C------------- LAST LINE OF EVCHCK FOLLOWS -----------------------------
      END
