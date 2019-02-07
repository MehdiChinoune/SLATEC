!*==DFCQX.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DFCQX
SUBROUTINE DFCQX(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--DFCQX5
  !***BEGIN PROLOGUE  DFCQX
  !***PURPOSE  Quick check for DFC.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (FCQX-S, DFCQX-D)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Hanson, R. J., (SNLA)
  !***DESCRIPTION
  !
  !   Quick check subprogram for the subroutine DFC.
  !
  !   Fit discrete data by an S-shaped curve.  Evaluate the fitted curve,
  !   its first two derivatives, and probable error curve.
  !
  !   Use subprogram DFC to obtain the constrained cubic B-spline
  !   representation of the curve.
  !
  !   The values of the coefficients of the B-spline as computed by DFC
  !   and the values of the fitted curve as computed by DBVALU in the
  !   de Boor package are tested for accuracy with the expected values.
  !   See the example program in the report sand78-1291, pp. 22-27.
  !
  !   The dimensions in the following arrays are as small as possible for
  !   the problem being solved.
  !
  !***ROUTINES CALLED  D1MACH, DBVALU, DCOPY, DCV, DFC, DMOUT, DVOUT,
  !                    IVOUT
  !***REVISION HISTORY  (YYMMDD)
  !   780801  DATE WRITTEN
  !   890718  Changed references from DBVLUE to DBVALU.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891004  Changed computation of XVAL.  (WRB)
  !   891004  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901010  Restructured using IF-THEN-ELSE-ENDIF, modified tolerances
  !           to use D1MACH(4) rather than D1MACH(3) and cleaned up
  !           FORMATs.  (RWC)
  !   930214  Declarations sections added, code revised to test error
  !           returns for all values of KPRINT and code polished.  (WRB)
  !***END PROLOGUE  DFCQX
  !     .. Scalar Arguments ..
  INTEGER Ipass, Kprint, Lun
  !     .. Local Scalars ..
  REAL(8) :: diff, one, t, tol, xval, zero
  INTEGER kontrl, i, idigit, ii, j, l, last, mode, n, nbkpt, &
    nconst, ndata, ndeg, nerr, nord, nval
  LOGICAL fatal
  !     .. Local Arrays ..
  REAL(8) :: bkpt(13), check(51), coefck(9), coeff(9), sddata(9), &
    v(51,5), w(529), work(12), xconst(11), xdata(9), &
    yconst(11), ydata(9)
  INTEGER iw(30), nderiv(11)
  !     .. External Functions ..
  REAL(8) :: D1MACH, DBVALU, DCV
  INTEGER NUMXER
  EXTERNAL DBVALU, DCV, NUMXER, D1MACH
  !     .. External Subroutines ..
  EXTERNAL DCOPY, DFC, DMOUT, DVOUT, IVOUT, XGETF, XSETF
  !     .. Intrinsic Functions ..
  INTRINSIC ABS, REAL, SQRT
  !     .. Data statements ..
  !
  DATA xdata(1), xdata(2), xdata(3), xdata(4), xdata(5), xdata(6), &
    xdata(7), xdata(8), xdata(9)/0.15D0, 0.27D0, 0.33D0, 0.40D0, &
    0.43D0, 0.47D0, 0.53D0, 0.58D0, 0.63D0/
  DATA ydata(1), ydata(2), ydata(3), ydata(4), ydata(5), ydata(6), &
    ydata(7), ydata(8), ydata(9)/0.025D0, 0.05D0, 0.13D0, 0.27D0, &
    0.37D0, 0.47D0, 0.64D0, 0.77D0, 0.87D0/
  DATA sddata(1)/0.015D0/, ndata/9/, nord/4/, nbkpt/13/, last/10/
  DATA bkpt(1), bkpt(2), bkpt(3), bkpt(4), bkpt(5), bkpt(6), bkpt(7), &
    bkpt(8), bkpt(9), bkpt(10), bkpt(11), bkpt(12), bkpt(13)&
    / - 0.6D0, -0.4D0, -0.2D0, 0.0D0, 0.2D0, 0.4D0, 0.6D0, 0.8D0, &
    0.9D0, 1.0D0, 1.1D0, 1.2D0, 1.3D0/
  !
  !     Store the data to be used to check the accuracy of the computed
  !     results.  See SAND78-1291, p.26.
  !
  DATA coefck(1), coefck(2), coefck(3), coefck(4), coefck(5), coefck(6)&
    , coefck(7), coefck(8), coefck(9)/1.186380846D-13, &
    -2.826166426D-14, -4.333929094D-15, 1.722113311D-01, &
    9.421965984D-01, 9.684708719D-01, 9.894902905D-01, &
    1.005254855D+00, 9.894902905D-01/
  DATA check(1), check(2), check(3), check(4), check(5), check(6), &
    check(7), check(8), check(9)/2.095830752D-16, 2.870188850D-05, &
    2.296151081D-04, 7.749509897D-04, 1.836920865D-03, &
    3.587736064D-03, 6.199607918D-03, 9.844747759D-03, &
    1.469536692D-02/
  DATA check(10), check(11), check(12), check(13), check(14), check(15)&
    , check(16), check(17), check(18)/2.092367672D-02, &
    2.870188851D-02, 3.824443882D-02, 4.993466504D-02, &
    6.419812979D-02, 8.146039566D-02, 1.021470253D-01, &
    1.266835812D-01, 1.554956261D-01/
  DATA check(19), check(20), check(21), check(22), check(23), check(24)&
    , check(25), check(26), check(27)/1.890087225D-01, &
    2.276484331D-01, 2.718403204D-01, 3.217163150D-01, &
    3.762338189D-01, 4.340566020D-01, 4.938484342D-01, &
    5.542730855D-01, 6.139943258D-01/
  DATA check(28), check(29), check(30), check(31), check(32), check(33)&
    , check(34), check(35), check(36)/6.716759250D-01, &
    7.259816530D-01, 7.755752797D-01, 8.191205752D-01, &
    8.556270903D-01, 8.854875002D-01, 9.094402609D-01, &
    9.282238286D-01, 9.425766596D-01/
  DATA check(37), check(38), check(39), check(40), check(41), check(42)&
    , check(43), check(44), check(45)/9.532372098D-01, &
    9.609439355D-01, 9.664352927D-01, 9.704497377D-01, &
    9.737257265D-01, 9.768786393D-01, 9.800315521D-01, &
    9.831844649D-01, 9.863373777D-01/
  DATA check(46), check(47), check(48), check(49), check(50), check(51)&
    /9.894902905D-01, 9.926011645D-01, 9.954598055D-01, &
    9.978139804D-01, 9.994114563D-01, 1.000000000D+00/
  !***FIRST EXECUTABLE STATEMENT  DFCQX
  IF ( Kprint>=2 ) WRITE (Lun,99001)
  !
  99001 FORMAT ('1'/' Test DFC')
  Ipass = 1
  !
  !     Broadcast SDDATA(1) value to all of SDDATA(*).
  !
  CALL DCOPY(ndata,sddata,0,sddata,1)
  zero = 0
  one = 1
  ndeg = nord - 1
  !
  !     Write the various constraints for the fitted curve.
  !
  nconst = 0
  t = bkpt(nord)
  !
  !     Constrain function to be zero at left-most breakpoint.
  !
  nconst = nconst + 1
  xconst(nconst) = t
  yconst(nconst) = zero
  nderiv(nconst) = 2 + 4*0
  !
  !     Constrain first derivative to be nonnegative at left-most
  !     breakpoint.
  !
  nconst = nconst + 1
  xconst(nconst) = t
  yconst(nconst) = zero
  nderiv(nconst) = 1 + 4*1
  !
  !     Constrain second derivatives to be nonnegative at left set of
  !     breakpoints.
  !
  DO i = 1, 3
    l = ndeg + i
    t = bkpt(l)
    nconst = nconst + 1
    xconst(nconst) = t
    yconst(nconst) = zero
    nderiv(nconst) = 1 + 4*2
  ENDDO
  !
  !     Constrain function value at right-most breakpoint to be one.
  !
  nconst = nconst + 1
  t = bkpt(last)
  xconst(nconst) = t
  yconst(nconst) = one
  nderiv(nconst) = 2 + 4*0
  !
  !     Constrain slope to agree at left- and right-most breakpoints.
  !
  nconst = nconst + 1
  xconst(nconst) = bkpt(nord)
  yconst(nconst) = bkpt(last)
  nderiv(nconst) = 3 + 4*1
  !
  !     Constrain second derivatives to be nonpositive at right set of
  !     breakpoints.
  !
  DO i = 1, 4
    nconst = nconst + 1
    l = last - 4 + i
    xconst(nconst) = bkpt(l)
    yconst(nconst) = zero
    nderiv(nconst) = 0 + 4*2
  ENDDO
  !
  idigit = -4
  !
  IF ( Kprint>=3 ) THEN
    CALL DVOUT(nbkpt,bkpt,'('' ARRAY OF KNOTS.'')',idigit)
    CALL DVOUT(ndata,xdata,'('' INDEPENDENT VARIABLE VALUES'')',idigit)
    CALL DVOUT(ndata,ydata,'('' DEPENDENT VARIABLE VALUES'')',idigit)
    CALL DVOUT(ndata,sddata,'('' DEPENDENT VARIABLE UNCERTAINTY'')',idigit)
    CALL DVOUT(nconst,xconst,'('' INDEPENDENT VARIABLE CONSTRAINT VALUES'')'&
      ,idigit)
    CALL DVOUT(nconst,yconst,'('' CONSTRAINT VALUES'')',idigit)
    CALL IVOUT(nconst,nderiv,'('' CONSTRAINT INDICATOR'')',idigit)
  ENDIF
  !
  !     Declare amount of working storage allocated to DFC.
  !
  iw(1) = 529
  iw(2) = 30
  !
  !     Set mode to indicate a new problem and request the variance
  !     function.
  !
  mode = 2
  !
  !     Obtain the coefficients of the B-spline.
  !
  CALL DFC(ndata,xdata,ydata,sddata,nord,nbkpt,bkpt,nconst,xconst,yconst,&
    nderiv,mode,coeff,w,iw)
  !
  !     Check coefficients.
  !
  tol = MAX(7.0D0*SQRT(D1MACH(4)),1.0D-8)
  diff = 0.0D0
  DO i = 1, ndata
    diff = MAX(diff,ABS(coeff(i)-coefck(i)))
  ENDDO
  IF ( diff<=tol ) THEN
    fatal = .FALSE.
    IF ( Kprint>=3 ) WRITE (Lun,99002)
    99002   FORMAT (/' DFC PASSED TEST 1')
  ELSE
    Ipass = 0
    fatal = .TRUE.
    IF ( Kprint>=2 ) WRITE (Lun,99003)
    99003   FORMAT (/' DFC FAILED TEST 1')
  ENDIF
  !
  IF ( (fatal.AND.Kprint>=2).OR.Kprint>=3 ) THEN
    CALL DVOUT(ndata,coefck,'(/'' PREDICTED COEFFICIENTS OF THE B-SPLINE '//&
      'FROM SAMPLE'')',idigit)
    CALL DVOUT(ndata,coeff,'(/'' COEFFICIENTS OF THE B-SPLINE COMPUTED '//&
      'BY DFC'')',idigit)
  ENDIF
  !
  !     Compute value, first two derivatives and probable uncertainty.
  !
  n = nbkpt - nord
  nval = 51
  DO i = 1, nval
    !
    !       The function DBVALU is in the de Boor B-spline package.
    !
    xval = REAL(i-1, 8)/(nval-1)
    ii = 1
    DO j = 1, 3
      v(i,j+1) = DBVALU(bkpt,coeff,n,nord,j-1,xval,ii,work)
    ENDDO
    v(i,1) = xval
    !
    !       The variance function DCV is a companion subprogram to DFC.
    !
    v(i,5) = SQRT(DCV(xval,ndata,nconst,nord,nbkpt,bkpt,w))
  ENDDO
  !
  diff = 0.0D0
  DO i = 1, nval
    diff = MAX(diff,ABS(v(i,2)-check(i)))
  ENDDO
  IF ( diff<=tol ) THEN
    fatal = .FALSE.
    IF ( Kprint>=3 ) WRITE (Lun,99004)
    99004   FORMAT (/' DFC (AND DBVALU) PASSED TEST 2')
  ELSE
    Ipass = 0
    fatal = .TRUE.
    IF ( Kprint>=2 ) WRITE (Lun,99005)
    99005   FORMAT (/' DFC (AND DBVALU) FAILED TEST 2')
  ENDIF
  !
  IF ( (fatal.AND.Kprint>=2).OR.Kprint>=3 ) THEN
    !
    !       Print these values.
    !
    CALL DMOUT(nval,5,nval,v,'(16X, ''X'', 10X, ''FNCN'', 8X,'//&
      '''1ST D'', 7X, ''2ND D'', 7X, ''ERROR'')',idigit)
    WRITE (Lun,99006)
    99006   FORMAT (/' VALUES SHOULD CORRESPOND TO THOSE IN ','SAND78-1291,',&
      ' P. 26')
  ENDIF
  !
  !     Trigger error conditions.
  !
  CALL XGETF(kontrl)
  IF ( Kprint<=2 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  ENDIF
  fatal = .FALSE.
  CALL XERCLR
  !
  IF ( Kprint>=3 ) WRITE (Lun,99007)
  99007 FORMAT (/' TRIGGER 6 ERROR MESSAGES',/)
  !
  CALL DFC(ndata,xdata,ydata,sddata,0,nbkpt,bkpt,nconst,xconst,yconst,&
    nderiv,mode,coeff,w,iw)
  IF ( NUMXER(nerr)/=2 ) fatal = .TRUE.
  CALL XERCLR
  !
  CALL DFC(ndata,xdata,ydata,sddata,nord,0,bkpt,nconst,xconst,yconst,nderiv,&
    mode,coeff,w,iw)
  IF ( NUMXER(nerr)/=2 ) fatal = .TRUE.
  CALL XERCLR
  !
  CALL DFC(-1,xdata,ydata,sddata,nord,nbkpt,bkpt,nconst,xconst,yconst,&
    nderiv,mode,coeff,w,iw)
  IF ( NUMXER(nerr)/=2 ) fatal = .TRUE.
  CALL XERCLR
  !
  mode = 0
  CALL DFC(ndata,xdata,ydata,sddata,nord,nbkpt,bkpt,nconst,xconst,yconst,&
    nderiv,mode,coeff,w,iw)
  IF ( NUMXER(nerr)/=2 ) fatal = .TRUE.
  CALL XERCLR
  !
  iw(1) = 10
  CALL DFC(ndata,xdata,ydata,sddata,nord,nbkpt,bkpt,nconst,xconst,yconst,&
    nderiv,mode,coeff,w,iw)
  IF ( NUMXER(nerr)/=2 ) fatal = .TRUE.
  CALL XERCLR
  !
  iw(1) = 529
  iw(2) = 2
  CALL DFC(ndata,xdata,ydata,sddata,nord,nbkpt,bkpt,nconst,xconst,yconst,&
    nderiv,mode,coeff,w,iw)
  IF ( NUMXER(nerr)/=2 ) fatal = .TRUE.
  CALL XERCLR
  !
  !     Restore KONTRL and check to see if the tests of error detection
  !     passed.
  !
  CALL XSETF(kontrl)
  IF ( fatal ) THEN
    Ipass = 0
    IF ( Kprint>=2 ) THEN
      WRITE (Lun,99008)
      99008     FORMAT (' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
    ENDIF
  ELSEIF ( Kprint>=3 ) THEN
    WRITE (Lun,99009)
    99009   FORMAT (' ALL INCORRECT ARGUMENT TESTS PASSED')
  ENDIF
  !
  !     Print PASS/FAIL message.
  !
  IF ( Ipass==1.AND.Kprint>=2 ) WRITE (Lun,99010)
  99010 FORMAT (/' ****************DFC PASSED ALL TESTS*****************')
  IF ( Ipass==0.AND.Kprint>=1 ) WRITE (Lun,99011)
  99011 FORMAT (/' ***************DFC FAILED SOME TESTS*****************')
  RETURN
END SUBROUTINE DFCQX
