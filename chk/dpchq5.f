*DECK DPCHQ5
      SUBROUTINE DPCHQ5 (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  DPCHQ5
C***PURPOSE  Test the PCH to B-spline conversion routine DPCHBS.
C***LIBRARY   SLATEC (PCHIP)
C***TYPE      DOUBLE PRECISION (PCHQK5-S, DPCHQ5-D)
C***KEYWORDS  PCHIP CONVERSION ROUTINE QUICK CHECK
C***AUTHOR  Fritsch, F. N., (LLNL)
C***DESCRIPTION
C
C             DPCHIP QUICK CHECK NUMBER 5
C
C     TESTS THE CONVERSION ROUTINE:  DPCHBS.
C *Usage:
C
C        INTEGER  LUN, KPRINT, IPASS
C
C        CALL DPCHQ5 (LUN, KPRINT, IPASS)
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
C   This routine tests a constructed data set with four different
C   KNOTYP settings.  It computes the function and derivatives of the
C   resulting B-representation via DBVALU and compares with PCH data.
C
C *Caution:
C   This routine assumes DBVALU has already been successfully tested.
C
C***ROUTINES CALLED  DBVALU, DPCHBS, D1MACH
C***REVISION HISTORY  (YYMMDD)
C   900411  DATE WRITTEN
C   900412  Corrected minor errors in initial implementation.
C   900430  Produced double precision version.
C   900501  Corrected declarations.
C   930317  Improved output formats.  (FNF)
C***END PROLOGUE  DPCHQ5
C
C*Internal Notes:
C  TOL  is the tolerance to use for quantities that should only
C       theoretically be equal.
C  TOLZ is the tolerance to use for quantities that should be exactly
C       equal.
C
C**End
C
C  Declare arguments.
C
      INTEGER  LUN, KPRINT, IPASS
C
C  Declare externals.
C
      DOUBLE PRECISION  DBVALU, D1MACH
      EXTERNAL  DBVALU, DPCHBS, D1MACH
C
C  Declare variables.
C
      INTEGER  I, IERR, IFAIL, INBV, J, KNOTYP, K, N, NDIM, NKNOTS
      PARAMETER  (N = 9)
      DOUBLE PRECISION  BCOEF(2*N), D(N), DCALC, DERR, DERMAX, F(N),
     *      FCALC, FERR, FERMAX, T(2*N+4), TERR, TERMAX, TOL, TOLZ,
     *      TSAVE(2*N+4), WORK(16*N), X(N), ZERO
      PARAMETER  (ZERO = 0.0D0)
      LOGICAL  FAIL
C
C  Define relative error function.
C
      DOUBLE PRECISION  ANS, ERR, RELERR
      RELERR (ERR, ANS) = ABS(ERR) / MAX(1.0D-5,ABS(ANS))
C
C  Define test data.
C
      DATA  X /-2.2D0,   -1.2D0,   -1.0D0,   -0.5D0,   -0.01D0,
     *          0.5D0,    1.0D0,    2.0D0,    2.2D0/
      DATA  F / 0.0079D0, 0.2369D0, 0.3679D0, 0.7788D0, 0.9999D0,
     *          0.7788D0, 0.3679D0, 0.1083D0, 0.0079D0/
      DATA  D / 0.0000D0, 0.3800D0, 0.7173D0, 0.5820D0, 0.0177D0,
     *         -0.5696D0,-0.5135D0,-0.0778D0,-0.0025D0/
C
C  Initialize.
C
C***FIRST EXECUTABLE STATEMENT  DPCHQ5
      IFAIL = 0
      TOL = 100*D1MACH(4)
      TOLZ = ZERO
C
      IF (KPRINT.GE.3)  WRITE (LUN, 1000)
      IF (KPRINT.GE.2)  WRITE (LUN, 1001)
C
C  Loop over a series of values of KNOTYP.
C
      IF (KPRINT.GE.3)  WRITE (LUN, 1010)
      DO 300  KNOTYP = 2, -1, -1
C        ------------
         CALL DPCHBS (N, X, F, D, 1, KNOTYP, NKNOTS, T, BCOEF, NDIM, K,
     *               IERR)
C        ------------
         IF (KPRINT.GE.3)
     *       WRITE (LUN, 2000) KNOTYP, NKNOTS, NDIM, K, IERR
         IF ( IERR.NE.0 )  THEN
            IFAIL = IFAIL + 1
            IF (KPRINT.GE.3)  WRITE (LUN, 2001)
         ELSE
C             Compare evaluated results with inputs to DPCHBS.
            INBV = 1
            FERMAX = ZERO
            DERMAX = ZERO
            IF (KPRINT.GE.3)  THEN
               WRITE (LUN, 2002)
               WRITE (LUN, 2003)  T(1), T(2)
               J = 1
            ENDIF
            DO 100  I = 1, N
               FCALC = DBVALU (T, BCOEF, NDIM, K, 0, X(I), INBV, WORK)
               FERR = F(I) - FCALC
               FERMAX = MAX(FERMAX, RELERR(FERR,F(I)) )
               DCALC = DBVALU (T, BCOEF, NDIM, K, 1, X(I), INBV, WORK)
               DERR = D(I) - DCALC
               DERMAX = MAX(DERMAX, RELERR(DERR,D(I)) )
               IF (KPRINT.GE.3)  THEN
                  J = J + 2
                  WRITE (LUN, 2004)  X(I), T(J), T(J+1), F(I), FERR,
     *                                                   D(I), DERR
               ENDIF
  100       CONTINUE
            IF (KPRINT.GE.3)  THEN
               J = J + 2
               WRITE (LUN, 2003)  T(J), T(J+1)
            ENDIF
            FAIL = (FERMAX.GT.TOL).OR.(DERMAX.GT.TOL)
            IF (FAIL)  IFAIL = IFAIL + 1
            IF ((KPRINT.GE.3).OR.(KPRINT.GE.2).AND.FAIL)
     *         WRITE (LUN, 2005)  FERMAX, DERMAX, TOL
         ENDIF
C
C          Special check for KNOTYP=-1.
         IF (KNOTYP.EQ.0)  THEN
C             Save knot vector for next test.
            DO 200  I = 1, NKNOTS
               TSAVE(I) = T(I)
  200       CONTINUE
         ELSE IF (KNOTYP.EQ.-1)  THEN
C             Check that knot vector is unchanged.
            TERMAX = ZERO
            DO 250  I = 1, NKNOTS
               TERR = ABS(T(I) - TSAVE(I))
               TERMAX = MAX(TERMAX, TERR)
  250       CONTINUE
            IF (TERMAX.GT.TOLZ)  THEN
               IFAIL = IFAIL + 1
               IF (KPRINT.GE.2)  WRITE (LUN, 2007)  TERMAX, TOLZ
            ENDIF
         ENDIF
  300 CONTINUE
C
C  PRINT SUMMARY AND TERMINATE.
C
      IF ((KPRINT.GE.2).AND.(IFAIL.NE.0))  WRITE (LUN, 3001)  IFAIL
C
      IF (IFAIL.EQ.0)  THEN
         IPASS = 1
         IF (KPRINT.GE.2) WRITE(LUN,99998)
      ELSE
         IPASS = 0
         IF (KPRINT.GE.1) WRITE(LUN,99999)
      ENDIF
C
      RETURN
C
C  FORMATS.
C
 1000 FORMAT ('1'//10X,'TEST PCH TO B-SPLINE CONVERTER')
 1001 FORMAT (//10X,'DPCHQ5 RESULTS'/10X,'--------------')
 1010 FORMAT (/4X,'(Results should be the same for all KNOTYP values.)')
 2000 FORMAT (/4X,'KNOTYP =',I2,':  NKNOTS =',I3,',  NDIM =',I3,
     *                        ',  K =',I2,',  IERR =',I3)
 2001 FORMAT (' *** Failed -- bad IERR value.')
 2002 FORMAT (/15X,'X',9X,'KNOTS',10X,'F',7X,'FERR',8X,'D',7X,'DERR')
 2003 FORMAT (18X,2F8.2)
 2004 FORMAT (10X,3F8.2,F10.4,1P,D10.2,0P,F10.4,1P,D10.2)
 2005 FORMAT (/5X,'Maximum relative errors:'
     *       /15X,'F-error =',1P,D13.5,5X,'D-error =',D13.5
     *        /5X,'Both should be less than  TOL =',D13.5)
 2007 FORMAT (/' *** T-ARRAY MAXIMUM CHANGE =',1P,D13.5,
     *           ';  SHOULD NOT EXCEED TOLZ =',D13.5)
 3001 FORMAT (/' *** TROUBLE ***',I5,' CONVERSION TESTS FAILED.')
99998 FORMAT (/' ------------ DPCHIP PASSED  ALL CONVERSION TESTS',
     *        ' ------------')
99999 FORMAT (/' ************ DPCHIP FAILED SOME CONVERSION TESTS',
     *        ' ************')
C------------- LAST LINE OF DPCHQ5 FOLLOWS -----------------------------
      END
