*DECK DFCQX
      SUBROUTINE DFCQX (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  DFCQX
C***PURPOSE  Quick check for DFC.
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (FCQX-S, DFCQX-D)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C   Quick check subprogram for the subroutine DFC.
C
C   Fit discrete data by an S-shaped curve.  Evaluate the fitted curve,
C   its first two derivatives, and probable error curve.
C
C   Use subprogram DFC to obtain the constrained cubic B-spline
C   representation of the curve.
C
C   The values of the coefficients of the B-spline as computed by DFC
C   and the values of the fitted curve as computed by DBVALU in the
C   de Boor package are tested for accuracy with the expected values.
C   See the example program in the report sand78-1291, pp. 22-27.
C
C   The dimensions in the following arrays are as small as possible for
C   the problem being solved.
C
C***ROUTINES CALLED  D1MACH, DBVALU, DCOPY, DCV, DFC, DMOUT, DVOUT,
C                    IVOUT
C***REVISION HISTORY  (YYMMDD)
C   780801  DATE WRITTEN
C   890718  Changed references from DBVLUE to DBVALU.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891004  Changed computation of XVAL.  (WRB)
C   891004  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   901010  Restructured using IF-THEN-ELSE-ENDIF, modified tolerances
C           to use D1MACH(4) rather than D1MACH(3) and cleaned up
C           FORMATs.  (RWC)
C   930214  Declarations sections added, code revised to test error
C           returns for all values of KPRINT and code polished.  (WRB)
C***END PROLOGUE  DFCQX
C     .. Scalar Arguments ..
      INTEGER IPASS, KPRINT, LUN
C     .. Local Scalars ..
      DOUBLE PRECISION DIFF, ONE, T, TOL, XVAL, ZERO
      INTEGER KONTRL, I, IDIGIT, II, J, L, LAST, MODE, N, NBKPT,
     +        NCONST, NDATA, NDEG, NERR, NORD, NVAL
      LOGICAL FATAL
C     .. Local Arrays ..
      DOUBLE PRECISION BKPT(13), CHECK(51), COEFCK(9), COEFF(9),
     +                 SDDATA(9), V(51,5), W(529), WORK(12), XCONST(11),
     +                 XDATA(9), YCONST(11), YDATA(9)
      INTEGER IW(30), NDERIV(11)
C     .. External Functions ..
      DOUBLE PRECISION D1MACH, DBVALU, DCV
      INTEGER NUMXER
      EXTERNAL DBVALU, DCV, NUMXER, D1MACH
C     .. External Subroutines ..
      EXTERNAL DCOPY, DFC, DMOUT, DVOUT, IVOUT, XGETF, XSETF
C     .. Intrinsic Functions ..
      INTRINSIC ABS, DBLE, SQRT
C     .. Data statements ..
C
      DATA XDATA(1),XDATA(2),XDATA(3),XDATA(4),XDATA(5),XDATA(6),
     +     XDATA(7),XDATA(8),XDATA(9)
     +     /0.15D0,0.27D0,0.33D0,0.40D0,0.43D0,0.47D0,
     +      0.53D0,0.58D0,0.63D0/
      DATA YDATA(1),YDATA(2),YDATA(3),YDATA(4),YDATA(5),YDATA(6),
     +     YDATA(7),YDATA(8),YDATA(9)
     +     /0.025D0,0.05D0,0.13D0,0.27D0,0.37D0,0.47D0,
     +      0.64D0,0.77D0,0.87D0/
      DATA SDDATA(1)/0.015D0/, NDATA/9/, NORD/4/, NBKPT/13/, LAST/10/
      DATA BKPT(1),BKPT(2),BKPT(3),BKPT(4),BKPT(5),BKPT(6),BKPT(7),
     +     BKPT(8),BKPT(9),BKPT(10),BKPT(11),BKPT(12),BKPT(13)
     +     /-0.6D0,-0.4D0,-0.2D0,0.0D0,0.2D0,0.4D0,0.6D0,
     +      0.8D0,0.9D0,1.0D0,1.1D0,1.2D0,1.3D0/
C
C     Store the data to be used to check the accuracy of the computed
C     results.  See SAND78-1291, p.26.
C
      DATA COEFCK(1),COEFCK(2),COEFCK(3),COEFCK(4),COEFCK(5),COEFCK(6),
     +     COEFCK(7),COEFCK(8),COEFCK(9)
     +     /1.186380846D-13,-2.826166426D-14,-4.333929094D-15,
     +      1.722113311D-01, 9.421965984D-01, 9.684708719D-01,
     +      9.894902905D-01, 1.005254855D+00, 9.894902905D-01/
      DATA CHECK(1),CHECK(2),CHECK(3),CHECK(4),CHECK(5),CHECK(6),
     +     CHECK(7),CHECK(8),CHECK(9)
     +     /2.095830752D-16, 2.870188850D-05, 2.296151081D-04,
     +      7.749509897D-04, 1.836920865D-03, 3.587736064D-03,
     +      6.199607918D-03, 9.844747759D-03, 1.469536692D-02/
      DATA CHECK(10),CHECK(11),CHECK(12),CHECK(13),CHECK(14),CHECK(15),
     +     CHECK(16),CHECK(17),CHECK(18)
     +     /2.092367672D-02, 2.870188851D-02, 3.824443882D-02,
     +      4.993466504D-02, 6.419812979D-02, 8.146039566D-02,
     +      1.021470253D-01, 1.266835812D-01, 1.554956261D-01/
      DATA CHECK(19),CHECK(20),CHECK(21),CHECK(22),CHECK(23),CHECK(24),
     +     CHECK(25),CHECK(26),CHECK(27)
     +     /1.890087225D-01, 2.276484331D-01, 2.718403204D-01,
     +      3.217163150D-01, 3.762338189D-01, 4.340566020D-01,
     +      4.938484342D-01, 5.542730855D-01,6.139943258D-01/
      DATA CHECK(28),CHECK(29),CHECK(30),CHECK(31),CHECK(32),CHECK(33),
     +     CHECK(34),CHECK(35),CHECK(36)
     +     /6.716759250D-01, 7.259816530D-01, 7.755752797D-01,
     +      8.191205752D-01, 8.556270903D-01, 8.854875002D-01,
     +      9.094402609D-01, 9.282238286D-01, 9.425766596D-01/
      DATA CHECK(37),CHECK(38),CHECK(39),CHECK(40),CHECK(41),CHECK(42),
     +     CHECK(43),CHECK(44),CHECK(45)
     +     /9.532372098D-01, 9.609439355D-01, 9.664352927D-01,
     +      9.704497377D-01, 9.737257265D-01, 9.768786393D-01,
     +      9.800315521D-01, 9.831844649D-01, 9.863373777D-01/
      DATA CHECK(46),CHECK(47),CHECK(48),CHECK(49),CHECK(50),
     +     CHECK(51)
     +     /9.894902905D-01, 9.926011645D-01, 9.954598055D-01,
     +      9.978139804D-01, 9.994114563D-01, 1.000000000D+00/
C***FIRST EXECUTABLE STATEMENT  DFCQX
      IF (KPRINT .GE. 2) WRITE (LUN,9000)
      IPASS = 1
C
C     Broadcast SDDATA(1) value to all of SDDATA(*).
C
      CALL DCOPY(NDATA,SDDATA,0,SDDATA,1)
      ZERO = 0
      ONE = 1
      NDEG = NORD - 1
C
C     Write the various constraints for the fitted curve.
C
      NCONST = 0
      T = BKPT(NORD)
C
C     Constrain function to be zero at left-most breakpoint.
C
      NCONST = NCONST + 1
      XCONST(NCONST) = T
      YCONST(NCONST) = ZERO
      NDERIV(NCONST) = 2 + 4*0
C
C     Constrain first derivative to be nonnegative at left-most
C     breakpoint.
C
      NCONST = NCONST + 1
      XCONST(NCONST) = T
      YCONST(NCONST) = ZERO
      NDERIV(NCONST) = 1 + 4*1
C
C     Constrain second derivatives to be nonnegative at left set of
C     breakpoints.
C
      DO 10 I = 1,3
        L = NDEG + I
        T = BKPT(L)
        NCONST = NCONST + 1
        XCONST(NCONST) = T
        YCONST(NCONST) = ZERO
        NDERIV(NCONST) = 1 + 4*2
   10 CONTINUE
C
C     Constrain function value at right-most breakpoint to be one.
C
      NCONST = NCONST + 1
      T = BKPT(LAST)
      XCONST(NCONST) = T
      YCONST(NCONST) = ONE
      NDERIV(NCONST) = 2 + 4*0
C
C     Constrain slope to agree at left- and right-most breakpoints.
C
      NCONST = NCONST + 1
      XCONST(NCONST) = BKPT(NORD)
      YCONST(NCONST) = BKPT(LAST)
      NDERIV(NCONST) = 3 + 4*1
C
C     Constrain second derivatives to be nonpositive at right set of
C     breakpoints.
C
      DO 20 I = 1,4
        NCONST = NCONST + 1
        L = LAST - 4 + I
        XCONST(NCONST) = BKPT(L)
        YCONST(NCONST) = ZERO
        NDERIV(NCONST) = 0 + 4*2
   20 CONTINUE
C
      IDIGIT = -4
C
      IF (KPRINT .GE. 3) THEN
        CALL DVOUT (NBKPT,BKPT,'('' ARRAY OF KNOTS.'')',IDIGIT)
        CALL DVOUT (NDATA,XDATA,
     +              '('' INDEPENDENT VARIABLE VALUES'')',IDIGIT)
        CALL DVOUT (NDATA,YDATA,'('' DEPENDENT VARIABLE VALUES'')',
     +              IDIGIT)
        CALL DVOUT (NDATA,SDDATA,
     +              '('' DEPENDENT VARIABLE UNCERTAINTY'')',IDIGIT)
        CALL DVOUT (NCONST,XCONST,
     +              '('' INDEPENDENT VARIABLE CONSTRAINT VALUES'')',
     +              IDIGIT)
        CALL DVOUT (NCONST,YCONST,'('' CONSTRAINT VALUES'')',IDIGIT)
        CALL IVOUT (NCONST,NDERIV,'('' CONSTRAINT INDICATOR'')',IDIGIT)
      ENDIF
C
C     Declare amount of working storage allocated to DFC.
C
      IW(1) = 529
      IW(2) = 30
C
C     Set mode to indicate a new problem and request the variance
C     function.
C
      MODE = 2
C
C     Obtain the coefficients of the B-spline.
C
      CALL DFC(NDATA,XDATA,YDATA,SDDATA,NORD,NBKPT,BKPT,NCONST,XCONST,
     +        YCONST,NDERIV,MODE,COEFF,W,IW)
C
C     Check coefficients.
C
      TOL = MAX(7.0D0*SQRT(D1MACH(4)),1.0D-8)
      DIFF = 0.0D0
      DO 30 I = 1,NDATA
        DIFF = MAX(DIFF,ABS(COEFF(I)-COEFCK(I)))
   30 CONTINUE
      IF (DIFF .LE. TOL) THEN
        FATAL = .FALSE.
        IF (KPRINT .GE. 3) WRITE (LUN,9010)
      ELSE
        IPASS = 0
        FATAL = .TRUE.
        IF (KPRINT .GE. 2) WRITE (LUN,9020)
      ENDIF
C
      IF ((FATAL.AND.KPRINT.GE.2) .OR. KPRINT.GE.3) THEN
        CALL DVOUT (NDATA,COEFCK,
     +              '(/'' PREDICTED COEFFICIENTS OF THE B-SPLINE '//
     +              'FROM SAMPLE'')',IDIGIT)
        CALL DVOUT (NDATA,COEFF,
     +              '(/'' COEFFICIENTS OF THE B-SPLINE COMPUTED '//
     +              'BY DFC'')',IDIGIT)
      ENDIF
C
C     Compute value, first two derivatives and probable uncertainty.
C
      N = NBKPT - NORD
      NVAL = 51
      DO 70 I = 1,NVAL
C
C       The function DBVALU is in the de Boor B-spline package.
C
        XVAL = DBLE(I-1)/(NVAL-1)
        II = 1
        DO 60 J = 1,3
          V(I,J+1) = DBVALU(BKPT,COEFF,N,NORD,J-1,XVAL,II,WORK)
   60   CONTINUE
        V(I,1) = XVAL
C
C       The variance function DCV is a companion subprogram to DFC.
C
        V(I,5) = SQRT(DCV(XVAL,NDATA,NCONST,NORD,NBKPT,BKPT,W))
   70 CONTINUE
C
      DIFF = 0.0D0
      DO 80 I = 1,NVAL
        DIFF = MAX(DIFF,ABS(V(I,2)-CHECK(I)))
   80 CONTINUE
      IF (DIFF .LE. TOL) THEN
        FATAL = .FALSE.
        IF (KPRINT .GE. 3) WRITE (LUN,9030)
      ELSE
        IPASS = 0
        FATAL = .TRUE.
        IF (KPRINT .GE. 2) WRITE (LUN,9040)
      ENDIF
C
      IF ((FATAL.AND.KPRINT.GE.2) .OR. KPRINT.GE.3) THEN
C
C       Print these values.
C
        CALL DMOUT(NVAL,5,NVAL,V,'(16X, ''X'', 10X, ''FNCN'', 8X,' //
     +             '''1ST D'', 7X, ''2ND D'', 7X, ''ERROR'')',IDIGIT)
        WRITE (LUN,9050)
      ENDIF
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
      IF (KPRINT .GE. 3) WRITE (LUN, 9060)
C
      CALL DFC(NDATA,XDATA,YDATA,SDDATA,0,NBKPT,BKPT,NCONST,XCONST,
     +        YCONST,NDERIV,MODE,COEFF,W,IW)
      IF (NUMXER(NERR) .NE. 2) FATAL = .TRUE.
      CALL XERCLR
C
      CALL DFC(NDATA,XDATA,YDATA,SDDATA,NORD,0,BKPT,NCONST,XCONST,
     +        YCONST,NDERIV,MODE,COEFF,W,IW)
      IF (NUMXER(NERR) .NE. 2) FATAL = .TRUE.
      CALL XERCLR
C
      CALL DFC(-1,XDATA,YDATA,SDDATA,NORD,NBKPT,BKPT,NCONST,XCONST,
     +        YCONST,NDERIV,MODE,COEFF,W,IW)
      IF (NUMXER(NERR) .NE. 2) FATAL = .TRUE.
      CALL XERCLR
C
      MODE = 0
      CALL DFC(NDATA,XDATA,YDATA,SDDATA,NORD,NBKPT,BKPT,NCONST,XCONST,
     +        YCONST,NDERIV,MODE,COEFF,W,IW)
      IF (NUMXER(NERR) .NE. 2) FATAL = .TRUE.
      CALL XERCLR
C
      IW(1) = 10
      CALL DFC(NDATA,XDATA,YDATA,SDDATA,NORD,NBKPT,BKPT,NCONST,XCONST,
     +        YCONST,NDERIV,MODE,COEFF,W,IW)
      IF (NUMXER(NERR) .NE. 2) FATAL = .TRUE.
      CALL XERCLR
C
      IW(1) = 529
      IW(2) = 2
      CALL DFC(NDATA,XDATA,YDATA,SDDATA,NORD,NBKPT,BKPT,NCONST,XCONST,
     +        YCONST,NDERIV,MODE,COEFF,W,IW)
      IF (NUMXER(NERR) .NE. 2) FATAL = .TRUE.
      CALL XERCLR
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
      RETURN
C
 9000 FORMAT ('1' / ' Test DFC')
 9010 FORMAT (/ ' DFC PASSED TEST 1')
 9020 FORMAT (/ ' DFC FAILED TEST 1')
 9030 FORMAT (/ ' DFC (AND DBVALU) PASSED TEST 2')
 9040 FORMAT (/ ' DFC (AND DBVALU) FAILED TEST 2')
 9050 FORMAT (/ ' VALUES SHOULD CORRESPOND TO THOSE IN ','SAND78-1291,',
     +        ' P. 26')
 9060 FORMAT (/ ' TRIGGER 6 ERROR MESSAGES',/)
 9070 FORMAT (' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
 9080 FORMAT (' ALL INCORRECT ARGUMENT TESTS PASSED')
 9100 FORMAT (/' ****************DFC PASSED ALL TESTS*****************')
 9110 FORMAT (/' ***************DFC FAILED SOME TESTS*****************')
      END
