*DECK FCQX
      SUBROUTINE FCQX (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  FCQX
C***PURPOSE  Quick check for FC.
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (FCQX-S, DFCQX-D)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C   Quick check subprogram for the subroutine FC.
C
C   Fit discrete data by an S-shaped curve.  Evaluate the fitted curve,
C   its first two derivatives, and probable error curve.
C
C   Use subprogram FC to obtain the constrained cubic B-spline
C   representation of the curve.
C
C   The values of the coefficients of the B-spline as computed by FC
C   and the values of the fitted curve as computed by BVALU in the
C   de Boor package are tested for accuracy with the expected values.
C   See the example program in the report sand78-1291, pp. 22-27.
C
C   The dimensions in the following arrays are as small as possible for
C   the problem being solved.
C
C***ROUTINES CALLED  BVALU, CV, FC, IVOUT, R1MACH, SCOPY, SMOUT, SVOUT
C***REVISION HISTORY  (YYMMDD)
C   780801  DATE WRITTEN
C   890718  Changed references from BVALUE to BVALU.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891004  Changed computation of XVAL.  (WRB)
C   891004  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   901010  Restructured using IF-THEN-ELSE-ENDIF, modified tolerances
C           to use R1MACH(4) rather than R1MACH(3) and cleaned up
C           FORMATs.  (RWC)
C   930214  Declarations sections added, code revised to test error
C           returns for all values of KPRINT and code polished.  (WRB)
C***END PROLOGUE  FCQX
C     .. Scalar Arguments ..
      INTEGER IPASS, KPRINT, LUN
C     .. Local Scalars ..
      REAL DIFF, ONE, T, TOL, XVAL, ZERO
      INTEGER KONTRL, I, IDIGIT, II, J, L, LAST, MODE, N, NBKPT,
     +        NCONST, NDATA, NDEG, NERR, NORD, NVAL
      LOGICAL FATAL
C     .. Local Arrays ..
      REAL BKPT(13), CHECK(51), COEFCK(9), COEFF(9), SDDATA(9), V(51,5),
     +     W(529), WORK(12), XCONST(11), XDATA(9), YCONST(11), YDATA(9)
      INTEGER IW(30), NDERIV(11)
C     .. External Functions ..
      REAL BVALU, CV, R1MACH
      INTEGER NUMXER
      EXTERNAL BVALU, CV, NUMXER, R1MACH
C     .. External Subroutines ..
      EXTERNAL FC, IVOUT, SCOPY, SMOUT, SVOUT, XGETF, XSETF
C     .. Intrinsic Functions ..
      INTRINSIC ABS, REAL, SQRT
C     .. Data statements ..
C
      DATA XDATA(1),XDATA(2),XDATA(3),XDATA(4),XDATA(5),XDATA(6),
     +     XDATA(7),XDATA(8),XDATA(9)
     +     /0.15E0,0.27E0,0.33E0,0.40E0,0.43E0,0.47E0,
     +      0.53E0,0.58E0,0.63E0/
      DATA YDATA(1),YDATA(2),YDATA(3),YDATA(4),YDATA(5),YDATA(6),
     +     YDATA(7),YDATA(8),YDATA(9)
     +     /0.025E0,0.05E0,0.13E0,0.27E0,0.37E0,0.47E0,
     +      0.64E0,0.77E0,0.87E0/
      DATA SDDATA(1)/0.015E0/, NDATA/9/, NORD/4/, NBKPT/13/, LAST/10/
      DATA BKPT(1),BKPT(2),BKPT(3),BKPT(4),BKPT(5),BKPT(6),BKPT(7),
     +     BKPT(8),BKPT(9),BKPT(10),BKPT(11),BKPT(12),BKPT(13)
     +     /-0.6E0,-0.4E0,-0.2E0,0.0E0,0.2E0,0.4E0,0.6E0,
     +      0.8E0,0.9E0,1.0E0,1.1E0,1.2E0,1.3E0/
C
C     Store the data to be used to check the accuracy of the computed
C     results.  See SAND78-1291, p.26.
C
      DATA COEFCK(1),COEFCK(2),COEFCK(3),COEFCK(4),COEFCK(5),COEFCK(6),
     +     COEFCK(7),COEFCK(8),COEFCK(9)
     +     /1.186380846E-13,-2.826166426E-14,-4.333929094E-15,
     +      1.722113311E-01, 9.421965984E-01, 9.684708719E-01,
     +      9.894902905E-01, 1.005254855E+00, 9.894902905E-01/
      DATA CHECK(1),CHECK(2),CHECK(3),CHECK(4),CHECK(5),CHECK(6),
     +     CHECK(7),CHECK(8),CHECK(9)
     +     /2.095830752E-16, 2.870188850E-05, 2.296151081E-04,
     +      7.749509897E-04, 1.836920865E-03, 3.587736064E-03,
     +      6.199607918E-03, 9.844747759E-03, 1.469536692E-02/
      DATA CHECK(10),CHECK(11),CHECK(12),CHECK(13),CHECK(14),CHECK(15),
     +     CHECK(16),CHECK(17),CHECK(18)
     +     /2.092367672E-02, 2.870188851E-02, 3.824443882E-02,
     +      4.993466504E-02, 6.419812979E-02, 8.146039566E-02,
     +      1.021470253E-01, 1.266835812E-01, 1.554956261E-01/
      DATA CHECK(19),CHECK(20),CHECK(21),CHECK(22),CHECK(23),CHECK(24),
     +     CHECK(25),CHECK(26),CHECK(27)
     +     /1.890087225E-01, 2.276484331E-01, 2.718403204E-01,
     +      3.217163150E-01, 3.762338189E-01, 4.340566020E-01,
     +      4.938484342E-01, 5.542730855E-01,6.139943258E-01/
      DATA CHECK(28),CHECK(29),CHECK(30),CHECK(31),CHECK(32),CHECK(33),
     +     CHECK(34),CHECK(35),CHECK(36)
     +     /6.716759250E-01, 7.259816530E-01, 7.755752797E-01,
     +      8.191205752E-01, 8.556270903E-01, 8.854875002E-01,
     +      9.094402609E-01, 9.282238286E-01, 9.425766596E-01/
      DATA CHECK(37),CHECK(38),CHECK(39),CHECK(40),CHECK(41),CHECK(42),
     +     CHECK(43),CHECK(44),CHECK(45)
     +     /9.532372098E-01, 9.609439355E-01, 9.664352927E-01,
     +      9.704497377E-01, 9.737257265E-01, 9.768786393E-01,
     +      9.800315521E-01, 9.831844649E-01, 9.863373777E-01/
      DATA CHECK(46),CHECK(47),CHECK(48),CHECK(49),CHECK(50),
     +     CHECK(51)
     +     /9.894902905E-01, 9.926011645E-01, 9.954598055E-01,
     +      9.978139804E-01, 9.994114563E-01, 1.000000000E+00/
C***FIRST EXECUTABLE STATEMENT  FCQX
      IF (KPRINT .GE. 2) WRITE (LUN,9000)
      IPASS = 1
C
C     Broadcast SDDATA(1) value to all of SDDATA(*).
C
      CALL SCOPY(NDATA,SDDATA,0,SDDATA,1)
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
        CALL SVOUT (NBKPT,BKPT,'('' ARRAY OF KNOTS.'')',IDIGIT)
        CALL SVOUT (NDATA,XDATA,
     +              '('' INDEPENDENT VARIABLE VALUES'')',IDIGIT)
        CALL SVOUT (NDATA,YDATA,'('' DEPENDENT VARIABLE VALUES'')',
     +              IDIGIT)
        CALL SVOUT (NDATA,SDDATA,
     +              '('' DEPENDENT VARIABLE UNCERTAINTY'')',IDIGIT)
        CALL SVOUT (NCONST,XCONST,
     +              '('' INDEPENDENT VARIABLE CONSTRAINT VALUES'')',
     +              IDIGIT)
        CALL SVOUT (NCONST,YCONST,'('' CONSTRAINT VALUES'')',IDIGIT)
        CALL IVOUT (NCONST,NDERIV,'('' CONSTRAINT INDICATOR'')',IDIGIT)
      ENDIF
C
C     Declare amount of working storage allocated to FC.
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
      CALL FC(NDATA,XDATA,YDATA,SDDATA,NORD,NBKPT,BKPT,NCONST,XCONST,
     +        YCONST,NDERIV,MODE,COEFF,W,IW)
C
C     Check coefficients.
C
      TOL = 7.0E0*SQRT(R1MACH(4))
      DIFF = 0.0E0
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
        CALL SVOUT (NDATA,COEFCK,
     +              '(/'' PREDICTED COEFFICIENTS OF THE B-SPLINE '//
     +              'FROM SAMPLE'')',IDIGIT)
        CALL SVOUT (NDATA,COEFF,
     +              '(/'' COEFFICIENTS OF THE B-SPLINE COMPUTED '//
     +              'BY FC'')',IDIGIT)
      ENDIF
C
C     Compute value, first two derivatives and probable uncertainty.
C
      N = NBKPT - NORD
      NVAL = 51
      DO 70 I = 1,NVAL
C
C       The function BVALU is in the de Boor B-spline package.
C
        XVAL = REAL(I-1)/(NVAL-1)
        II = 1
        DO 60 J = 1,3
          V(I,J+1) = BVALU(BKPT,COEFF,N,NORD,J-1,XVAL,II,WORK)
   60   CONTINUE
        V(I,1) = XVAL
C
C       The variance function CV is a companion subprogram to FC.
C
        V(I,5) = SQRT(CV(XVAL,NDATA,NCONST,NORD,NBKPT,BKPT,W))
   70 CONTINUE
C
      DIFF = 0.0E0
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
        CALL SMOUT(NVAL,5,NVAL,V,'(16X, ''X'', 10X, ''FNCN'', 8X,' //
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
      CALL FC(NDATA,XDATA,YDATA,SDDATA,0,NBKPT,BKPT,NCONST,XCONST,
     +        YCONST,NDERIV,MODE,COEFF,W,IW)
      IF (NUMXER(NERR) .NE. 2) FATAL = .TRUE.
      CALL XERCLR
C
      CALL FC(NDATA,XDATA,YDATA,SDDATA,NORD,0,BKPT,NCONST,XCONST,
     +        YCONST,NDERIV,MODE,COEFF,W,IW)
      IF (NUMXER(NERR) .NE. 2) FATAL = .TRUE.
      CALL XERCLR
C
      CALL FC(-1,XDATA,YDATA,SDDATA,NORD,NBKPT,BKPT,NCONST,XCONST,
     +        YCONST,NDERIV,MODE,COEFF,W,IW)
      IF (NUMXER(NERR) .NE. 2) FATAL = .TRUE.
      CALL XERCLR
C
      MODE = 0
      CALL FC(NDATA,XDATA,YDATA,SDDATA,NORD,NBKPT,BKPT,NCONST,XCONST,
     +        YCONST,NDERIV,MODE,COEFF,W,IW)
      IF (NUMXER(NERR) .NE. 2) FATAL = .TRUE.
      CALL XERCLR
C
      IW(1) = 10
      CALL FC(NDATA,XDATA,YDATA,SDDATA,NORD,NBKPT,BKPT,NCONST,XCONST,
     +        YCONST,NDERIV,MODE,COEFF,W,IW)
      IF (NUMXER(NERR) .NE. 2) FATAL = .TRUE.
      CALL XERCLR
C
      IW(1) = 529
      IW(2) = 2
      CALL FC(NDATA,XDATA,YDATA,SDDATA,NORD,NBKPT,BKPT,NCONST,XCONST,
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
 9000 FORMAT ('1' / ' Test FC')
 9010 FORMAT (/ ' FC PASSED TEST 1')
 9020 FORMAT (/ ' FC FAILED TEST 1')
 9030 FORMAT (/ ' FC (AND BVALU) PASSED TEST 2')
 9040 FORMAT (/ ' FC (AND BVALU) FAILED TEST 2')
 9050 FORMAT (/ ' VALUES SHOULD CORRESPOND TO THOSE IN ','SAND78-1291,',
     +        ' P. 26')
 9060 FORMAT (/ ' TRIGGER 6 ERROR MESSAGES',/)
 9070 FORMAT (' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
 9080 FORMAT (' ALL INCORRECT ARGUMENT TESTS PASSED')
 9100 FORMAT (/' *****************FC PASSED ALL TESTS*****************')
 9110 FORMAT (/' ****************FC FAILED SOME TESTS*****************')
      END
