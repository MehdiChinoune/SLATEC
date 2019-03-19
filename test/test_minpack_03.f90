MODULE TEST52_MOD
  IMPLICIT NONE

CONTAINS
  !DECK CMPARE
  SUBROUTINE CMPARE(Icnt,Itest)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  CMPARE
    !***PURPOSE  Compare values in COMMON block CHECK for quick check
    !            routine PFITQX.
    !***LIBRARY   SLATEC
    !***TYPE      SINGLE PRECISION (CMPARE-S, DCMPAR-D)
    !***AUTHOR  (UNKNOWN)
    !***ROUTINES CALLED  (NONE)
    !***COMMON BLOCKS    CHECK
    !***REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   890921  Realigned order of variables in the COMMON block.
    !           (WRB)
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   920214  Minor improvements to code for readability.  (WRB)
    !***END PROLOGUE  CMPARE
    !     .. Scalar Arguments ..
    INTEGER Icnt
    !     .. Array Arguments ..
    INTEGER Itest(9)
    !     .. Scalars in Common ..
    REAL EPS, RP, SVEps, TOL
    INTEGER IERp, IERr, NORd, NORdp
    !     .. Arrays in Common ..
    REAL R(11)
    !     .. Local Scalars ..
    REAL rpp, ss
    INTEGER ierpp, nrdp
    !     .. Local Arrays ..
    INTEGER itemp(4)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS
    !     .. Common blocks ..
    COMMON /CHECK / EPS, R, RP, SVEps, TOL, NORdp, NORd, IERp, IERr
    !***FIRST EXECUTABLE STATEMENT  CMPARE
    Icnt = Icnt + 1
    itemp(1) = 0
    itemp(2) = 0
    itemp(3) = 0
    itemp(4) = 0
    ss = SVEps - EPS
    nrdp = NORdp - NORd
    rpp = RP - R(11)
    ierpp = IERp - IERr
    IF ( ABS(ss)<=TOL.OR.Icnt<=2.OR.Icnt>=6 ) itemp(1) = 1
    IF ( ABS(nrdp)==0 ) itemp(2) = 1
    IF ( ABS(rpp)<=TOL ) itemp(3) = 1
    IF ( ABS(ierpp)==0 ) itemp(4) = 1
    !
    !     Check to see if all four tests were good.
    !     If so, set the test number equal to 1.
    !
    Itest(Icnt) = itemp(1)*itemp(2)*itemp(3)*itemp(4)
  END SUBROUTINE CMPARE
  !DECK PFITQX
  SUBROUTINE PFITQX(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  PFITQX
    !***PURPOSE  Quick check for POLFIT, PCOEF and PVALUE.
    !***LIBRARY   SLATEC
    !***TYPE      SINGLE PRECISION (PFITQX-S, DPFITT-D)
    !***AUTHOR  (UNKNOWN)
    !***ROUTINES CALLED  CMPARE, PASS, PCOEF, POLFIT, PVALUE, R1MACH,
    !                    XERCLR, XGETF, XSETF
    !***COMMON BLOCKS    CHECK
    !***REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   890921  Realigned order of variables in the COMMON block.
    !           (WRB)
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   900911  Test problem changed and cosmetic changes to code.  (WRB)
    !   901205  Changed usage of R1MACH(3) to R1MACH(4) and modified the
    !           FORMATs.  (RWC)
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    !   920214  Code restructured to test for all values of KPRINT and to
    !           provide more PASS/FAIL information.  (WRB)
    !***END PROLOGUE  PFITQX
    INTEGER kontrl
    !     .. Scalar Arguments ..
    INTEGER Ipass, Kprint, Lun
    !     .. Scalars in Common ..
    REAL EPS, RP, SVEps, TOL
    INTEGER IERp, IERr, NORd, NORdp
    !     .. Arrays in Common ..
    REAL R(11)
    !     .. Local Scalars ..
    REAL yfit
    INTEGER i, icnt, m, maxord
    !     .. Local Arrays ..
    REAL a(97), tc(5), w(11), x(11), y(11), yp(5)
    INTEGER itest(9)
    !     .. External Functions ..
    REAL R1MACH
    EXTERNAL R1MACH
    !     .. External Subroutines ..
    EXTERNAL PASS, PCOEF, POLFIT, PVALUE
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, SQRT
    !     .. Common blocks ..
    COMMON /CHECK / EPS, R, RP, SVEps, TOL, NORdp, NORd, IERp, IERr
    !***FIRST EXECUTABLE STATEMENT  PFITQX
    IF ( Kprint>=2 ) WRITE (Lun,FMT=99002)
    !
    !     Initialize variables for testing passage or failure of tests
    !
    DO i = 1, 9
      itest(i) = 0
    ENDDO
    icnt = 0
    TOL = SQRT(R1MACH(4))
    m = 11
    DO i = 1, m
      x(i) = i - 6
      y(i) = x(i)**4
    ENDDO
    !
    !     Test POLFIT
    !     Input EPS is negative - specified level
    !
    w(1) = -1.0E0
    EPS = -0.01E0
    SVEps = EPS
    maxord = 8
    NORdp = 4
    RP = 625.0E0
    IERp = 1
    CALL POLFIT(m,x,y,w,maxord,NORd,EPS,R,IERr,a)
    !
    !     See if test passed
    !
    CALL CMPARE(icnt,itest)
    !
    !     Check for suppression of printing.
    !
    IF ( Kprint/=0 ) THEN
      IF ( Kprint/=1.OR.itest(icnt)/=1 ) THEN
        WRITE (Lun,FMT=99003)
        WRITE (Lun,FMT=99004)
        IF ( Kprint>2.OR.itest(icnt)/=1 ) THEN
          WRITE (Lun,FMT=99005) SVEps, NORdp, RP, IERp
          WRITE (Lun,FMT=99006) EPS, NORd, R(11), IERr
        ENDIF
        !
        !     Send message indicating passage or failure of test
        !
        CALL PASS(Lun,icnt,itest(icnt))
      ENDIF
    ENDIF
    !
    !     Input EPS is negative - computed level
    !
    EPS = -1.0E0
    SVEps = EPS
    CALL POLFIT(m,x,y,w,maxord,NORd,EPS,R,IERr,a)
    !
    !     See if test passed
    !
    CALL CMPARE(icnt,itest)
    !
    !     Check for suppression of printing.
    !
    IF ( Kprint/=0 ) THEN
      IF ( Kprint/=1.OR.itest(icnt)/=1 ) THEN
        WRITE (Lun,FMT=99007)
        IF ( Kprint>2.OR.itest(icnt)/=1 ) THEN
          WRITE (Lun,FMT=99008) maxord
          WRITE (Lun,FMT=99005) SVEps, NORdp, RP, IERp
          WRITE (Lun,FMT=99006) EPS, NORd, R(11), IERr
        ENDIF
        !
        !     Send message indicating passage or failure of test
        !
        CALL PASS(Lun,icnt,itest(icnt))
      ENDIF
    ENDIF
    !
    !     Input EPS is zero
    !
    w(1) = -1.0E0
    EPS = 0.0E0
    SVEps = EPS
    NORdp = 5
    maxord = 5
    CALL POLFIT(m,x,y,w,maxord,NORd,EPS,R,IERr,a)
    !
    !     See if test passed
    !
    CALL CMPARE(icnt,itest)
    !
    !     Check for suppression of printing.
    !
    IF ( Kprint/=0 ) THEN
      IF ( Kprint/=1.OR.itest(icnt)/=1 ) THEN
        WRITE (Lun,FMT=99009)
        IF ( Kprint>2.OR.itest(icnt)/=1 ) THEN
          WRITE (Lun,FMT=99008) maxord
          WRITE (Lun,FMT=99005) SVEps, NORdp, RP, IERp
          WRITE (Lun,FMT=99006) EPS, NORd, R(11), IERr
        ENDIF
        !
        !     Send message indicating passage or failure of test
        !
        CALL PASS(Lun,icnt,itest(icnt))
      ENDIF
    ENDIF
    !
    !     Input EPS is positive
    !
    IERp = 1
    NORdp = 4
    EPS = 75.0E0*R1MACH(4)
    SVEps = EPS
    CALL POLFIT(m,x,y,w,maxord,NORd,EPS,R,IERr,a)
    !
    !     See if test passed
    !
    CALL CMPARE(icnt,itest)
    !
    !     Check for suppression of printing.
    !
    IF ( Kprint/=0 ) THEN
      IF ( Kprint/=1.OR.itest(icnt)/=1 ) THEN
        WRITE (Lun,FMT=99010)
        IF ( Kprint>2.OR.itest(icnt)/=1 ) THEN
          WRITE (Lun,FMT=99008) maxord
          WRITE (Lun,FMT=99005) SVEps, NORdp, RP, IERp
          WRITE (Lun,FMT=99006) EPS, NORd, R(11), IERr
        ENDIF
        !
        !     Send message indicating passage or failure of test
        !
        CALL PASS(Lun,icnt,itest(icnt))
      ENDIF
    ENDIF
    !
    !     Improper input
    !
    IERp = 2
    m = -2
    !
    !     Check for suppression of printing.
    !
    CALL XGETF(kontrl)
    IF ( Kprint<=2 ) THEN
      CALL XSETF(0)
    ELSE
      CALL XSETF(1)
    ENDIF
    CALL XERCLR
    !
    IF ( Kprint>=3 ) WRITE (Lun,99001)
    99001 FORMAT (/' Invalid input')
    CALL POLFIT(m,x,y,w,maxord,NORd,EPS,R,IERr,a)
    !
    !     See if test passed
    !
    icnt = icnt + 1
    IF ( IERr==2 ) THEN
      itest(icnt) = 1
      IF ( Kprint>=3 ) WRITE (Lun,99011) 'PASSED', IERr
    ELSEIF ( Kprint>=2 ) THEN
      WRITE (Lun,99011) 'FAILED', IERr
    ENDIF
    !
    !     Check for suppression of printing.
    !
    IF ( Kprint/=0 ) THEN
      IF ( Kprint/=1.OR.itest(icnt)/=1 ) THEN
        IF ( Kprint<=2.AND.itest(icnt)==1 ) THEN
        ENDIF
        !
        !     Send message indicating passage or failure of test
        !
        CALL PASS(Lun,icnt,itest(icnt))
        !
        CALL XERCLR
        CALL XSETF(kontrl)
      ENDIF
    ENDIF
    !
    !     MAXORD too small to meet RMS error
    !
    m = 11
    w(1) = -1.0E0
    EPS = 5.0E0*R1MACH(4)
    SVEps = EPS
    RP = 553.0E0
    maxord = 2
    IERp = 3
    NORdp = 2
    CALL POLFIT(m,x,y,w,maxord,NORd,EPS,R,IERr,a)
    !
    !     See if test passed
    !
    CALL CMPARE(icnt,itest)
    !
    !     Check for suppression of printing.
    !
    IF ( Kprint/=0 ) THEN
      IF ( Kprint/=1.OR.itest(icnt)/=1 ) THEN
        WRITE (Lun,FMT=99012)
        IF ( Kprint>2.OR.itest(icnt)/=1 ) THEN
          WRITE (Lun,FMT=99008) maxord
          WRITE (Lun,FMT=99005) SVEps, NORdp, RP, IERp
          WRITE (Lun,FMT=99006) EPS, NORd, R(11), IERr
        ENDIF
        !
        !     Send message indicating passage or failure of test
        !
        CALL PASS(Lun,icnt,itest(icnt))
      ENDIF
    ENDIF
    !
    !     MAXORD too small to meet statistical test
    !
    NORdp = 4
    IERp = 4
    RP = 625.0E0
    EPS = -0.01E0
    SVEps = EPS
    maxord = 5
    CALL POLFIT(m,x,y,w,maxord,NORd,EPS,R,IERr,a)
    !
    !     See if test passed
    !
    CALL CMPARE(icnt,itest)
    !
    !     Check for suppression of printing.
    !
    IF ( Kprint/=0 ) THEN
      IF ( Kprint/=1.OR.itest(icnt)/=1 ) THEN
        WRITE (Lun,FMT=99013)
        IF ( Kprint>2.OR.itest(icnt)/=1 ) THEN
          WRITE (Lun,FMT=99008) maxord
          WRITE (Lun,FMT=99005) SVEps, NORdp, RP, IERp
          WRITE (Lun,FMT=99006) EPS, NORd, R(11), IERr
        ENDIF
        !
        !     Send message indicating passage or failure of test
        !
        CALL PASS(Lun,icnt,itest(icnt))
      ENDIF
    ENDIF
    !
    !     Test PCOEF
    !
    maxord = 6
    EPS = 0.0E0
    SVEps = EPS
    y(6) = 1.0E0
    DO i = 1, m
      w(i) = 1.0E0/(y(i)**2)
    ENDDO
    y(6) = 0.0E0
    CALL POLFIT(m,x,y,w,maxord,NORd,EPS,R,IERr,a)
    CALL PCOEF(4,5.0E0,tc,a)
    !
    !     See if test passed
    !
    icnt = icnt + 1
    IF ( ABS(R(11)-tc(1))<=TOL ) itest(icnt) = 1
    !
    !     Check for suppression of printing
    !
    IF ( Kprint/=0 ) THEN
      IF ( Kprint/=1.OR.itest(icnt)/=1 ) THEN
        WRITE (Lun,FMT=99014)
        IF ( Kprint>2.OR.itest(icnt)/=1 ) WRITE (Lun,FMT=99015) R(11), tc(1)
        !
        !     Send message indicating passage or failure of test
        !
        CALL PASS(Lun,icnt,itest(icnt))
      ENDIF
    ENDIF
    !
    !     Test PVALUE
    !     Normal call
    !
    CALL PVALUE(6,0,x(8),yfit,yp,a)
    !
    !     See if test passed
    !
    icnt = icnt + 1
    IF ( ABS(R(8)-yfit)<=TOL ) itest(icnt) = 1
    !
    !     Check for suppression of printing
    !
    IF ( Kprint/=0 ) THEN
      IF ( Kprint/=1.OR.itest(icnt)/=1 ) THEN
        WRITE (Lun,FMT=99016)
        WRITE (Lun,FMT=99017)
        IF ( Kprint>2.OR.itest(icnt)/=1 ) WRITE (Lun,FMT=99018) x(8), R(8), &
          yfit
        !
        !     Send message indicating passage or failure of test
        !
        CALL PASS(Lun,icnt,itest(icnt))
      ENDIF
    ENDIF
    !
    !     Check to see if all tests passed
    !
    Ipass = 1
    DO i = 1, 9
      Ipass = Ipass*itest(i)
    ENDDO
    !
    IF ( Ipass==1.AND.Kprint>=3 ) WRITE (Lun,FMT=99019)
    IF ( Ipass==0.AND.Kprint>=2 ) WRITE (Lun,FMT=99020)
    RETURN
    !
    99002 FORMAT ('1'/' Test POLFIT, PCOEF and PVALUE')
    99003 FORMAT (' Exercise POLFIT')
    99004 FORMAT (' Input EPS is negative - specified significance level')
    99005 FORMAT (' Input EPS =  ',E15.8,'   correct order =  ',I3,'   R(1) = ',&
      E15.8,'   IERR = ',I1)
    99006 FORMAT (' Output EPS = ',E15.8,'   computed order = ',I3,'   R(1) = ',&
      E15.8,'   IERR = ',I1)
    99007 FORMAT (/' Input EPS is negative - computed significance level')
    99008 FORMAT (' Maximum order = ',I2)
    99009 FORMAT (/' Input EPS is zero')
    99010 FORMAT (/' Input EPS is positive')
    99011 FORMAT (' POLFIT incorrect argument test ',A/' IERR should be 2.  It is ',&
      I4)
    99012 FORMAT (/' Cannot meet RMS error requirement')
    99013 FORMAT (/' Cannot satisfy statistical test')
    99014 FORMAT (/' Exercise PCOEF')
    99015 FORMAT (/' For C=1.0, correct coefficient = ',E15.8,'   computed = ',&
      E15.8)
    99016 FORMAT (/' Exercise PVALUE')
    99017 FORMAT (' Normal execution')
    99018 FORMAT (' For X = ',F5.2,'   correct P(X) = ',E15.8,&
      '    P(X) from PVALUE = ',E15.8)
    99019 FORMAT (/' ***************POLFIT PASSED ALL TESTS***************')
    99020 FORMAT (/' ***************POLFIT FAILED SOME TESTS**************')
  END SUBROUTINE PFITQX
  !DECK SNLS1Q
  SUBROUTINE SNLS1Q(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  SNLS1Q
    !***PURPOSE  Quick check for SNLS1E, SNLS1 and SCOV.
    !***LIBRARY   SLATEC
    !***TYPE      SINGLE PRECISION (SNLS1Q-S, DNLS1Q-D)
    !***KEYWORDS  QUICK CHECK
    !***AUTHOR  (UNKNOWN)
    !***DESCRIPTION
    !
    !   This subroutine performs a quick check on the subroutines SNLS1E
    !   (and SNLS1) and SCOV.
    !
    !***ROUTINES CALLED  ENORM, FCN1, FCN2, FCN3, FDJAC3, PASS, R1MACH,
    !                    SCOV, SNLS1E
    !***REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   890911  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   930214  Declarations sections added, code revised to test error
    !           returns for all values of KPRINT and code polished.  (WRB)
    !***END PROLOGUE  SNLS1Q
    !     .. Scalar Arguments ..
    INTEGER Ipass, Kprint, Lun
    !     .. Local Scalars ..
    REAL fnorm, fnorms, one, sigma, temp1, temp2, temp3, tol, tol2, &
      zero
    INTEGER i, iflag, info, infos, iopt, kontrl, ldfjac, lwa, m, n, &
      nerr, nprint
    LOGICAL fatal
    !     .. Local Arrays ..
    REAL fjac(10,2), fjrow(2), fjtj(3), fvec(10), wa(40), x(2)
    INTEGER iw(2)
    !     .. External Functions ..
    REAL ENORM, R1MACH
    INTEGER NUMXER
    EXTERNAL ENORM, NUMXER, R1MACH
    !     .. External Subroutines ..
    EXTERNAL FDJAC3, PASS, SCOV, SNLS1E, XGETF, XSETF
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, SQRT
    !***FIRST EXECUTABLE STATEMENT  SNLS1Q
    IF ( Kprint>=2 ) WRITE (Lun,99001)
    !
    99001 FORMAT ('1'/' Test SNLS1E, SNLS1 and SCOV')
    !
    Ipass = 1
    infos = 1
    fnorms = 1.1151779E+01
    m = 10
    n = 2
    lwa = 40
    ldfjac = 10
    nprint = -1
    iflag = 1
    zero = 0.0E0
    one = 1.0E0
    tol = SQRT(40.0E0*R1MACH(4))
    tol2 = SQRT(tol)
    !
    !     OPTION=2, the full Jacobian is stored and the user provides the
    !     Jacobian.
    !
    iopt = 2
    x(1) = 3.0E-1
    x(2) = 4.0E-1
    CALL SNLS1E(FCN2,iopt,m,n,x,fvec,tol,nprint,info,iw,wa,lwa)
    fnorm = ENORM(m,fvec)
    IF ( info==infos.AND.ABS(fnorm-fnorms)/fnorms<=tol ) THEN
      fatal = .FALSE.
      IF ( Kprint>=3 ) CALL PASS(Lun,1,1)
    ELSE
      Ipass = 0
      fatal = .TRUE.
      IF ( Kprint>=2 ) CALL PASS(Lun,1,0)
    ENDIF
    IF ( (fatal.AND.Kprint>=2).OR.Kprint>=3 ) WRITE (Lun,99007) infos, &
      fnorms, info, fnorm
    !
    !     Form JAC-transpose*JAC.
    !
    sigma = fnorm*fnorm/(m-n)
    iflag = 2
    CALL FCN2(iflag,m,n,x,fvec,fjac,ldfjac)
    DO i = 1, 3
      fjtj(i) = zero
    ENDDO
    DO i = 1, m
      fjtj(1) = fjtj(1) + fjac(i,1)**2
      fjtj(2) = fjtj(2) + fjac(i,1)*fjac(i,2)
      fjtj(3) = fjtj(3) + fjac(i,2)**2
    ENDDO
    !
    !     Calculate the covariance matrix.
    !
    CALL SCOV(FCN2,iopt,m,n,x,fvec,fjac,ldfjac,info,wa(1),wa(n+1),wa(2*n+1),&
      wa(3*n+1))
    !
    !     Form JAC-transpose*JAC * covariance matrix (should = SIGMA*I).
    !
    temp1 = (fjtj(1)*fjac(1,1)+fjtj(2)*fjac(1,2))/sigma
    temp2 = (fjtj(1)*fjac(1,2)+fjtj(2)*fjac(2,2))/sigma
    temp3 = (fjtj(2)*fjac(1,2)+fjtj(3)*fjac(2,2))/sigma
    IF ( info==infos.AND.ABS(temp1-one)<tol2.AND.ABS(temp2)<tol2.AND.&
        ABS(temp3-one)<tol2 ) THEN
      fatal = .FALSE.
      IF ( Kprint>=3 ) CALL PASS(Lun,2,1)
    ELSE
      Ipass = 0
      fatal = .TRUE.
      IF ( Kprint>=2 ) CALL PASS(Lun,2,0)
    ENDIF
    IF ( (fatal.AND.Kprint>=2).OR.Kprint>=3 ) WRITE (Lun,99008) infos, info, &
      temp1, temp2, temp3
    !
    !     OPTION=1, the full Jacobian is stored and the code approximates
    !     the Jacobian.
    !
    iopt = 1
    x(1) = 3.0E-1
    x(2) = 4.0E-1
    CALL SNLS1E(FCN1,iopt,m,n,x,fvec,tol,nprint,info,iw,wa,lwa)
    fnorm = ENORM(m,fvec)
    IF ( info==infos.AND.ABS(fnorm-fnorms)/fnorms<=tol ) THEN
      fatal = .FALSE.
      IF ( Kprint>=3 ) CALL PASS(Lun,3,1)
    ELSE
      Ipass = 0
      fatal = .TRUE.
      IF ( Kprint>=2 ) CALL PASS(Lun,3,0)
    ENDIF
    IF ( (fatal.AND.Kprint>=2).OR.Kprint>=3 ) WRITE (Lun,99007) infos, &
      fnorms, info, fnorm
    !
    !     Form JAC-transpose*JAC.
    !
    sigma = fnorm*fnorm/(m-n)
    iflag = 1
    CALL FDJAC3(FCN1,m,n,x,fvec,fjac,ldfjac,iflag,zero,wa)
    DO i = 1, 3
      fjtj(i) = zero
    ENDDO
    DO i = 1, m
      fjtj(1) = fjtj(1) + fjac(i,1)**2
      fjtj(2) = fjtj(2) + fjac(i,1)*fjac(i,2)
      fjtj(3) = fjtj(3) + fjac(i,2)**2
    ENDDO
    !
    !     Calculate the covariance matrix.
    !
    CALL SCOV(FCN1,iopt,m,n,x,fvec,fjac,ldfjac,info,wa(1),wa(n+1),wa(2*n+1),&
      wa(3*n+1))
    !
    !     Form JAC-transpose*JAC * covariance matrix (should = SIGMA*I).
    !
    temp1 = (fjtj(1)*fjac(1,1)+fjtj(2)*fjac(1,2))/sigma
    temp2 = (fjtj(1)*fjac(1,2)+fjtj(2)*fjac(2,2))/sigma
    temp3 = (fjtj(2)*fjac(1,2)+fjtj(3)*fjac(2,2))/sigma
    IF ( info==infos.AND.ABS(temp1-one)<tol2.AND.ABS(temp2)<tol2.AND.&
        ABS(temp3-one)<tol2 ) THEN
      fatal = .FALSE.
      IF ( Kprint>=3 ) CALL PASS(Lun,4,1)
    ELSE
      Ipass = 0
      fatal = .TRUE.
      IF ( Kprint>=2 ) CALL PASS(Lun,4,0)
    ENDIF
    IF ( (fatal.AND.Kprint>=2).OR.Kprint>=3 ) WRITE (Lun,99008) infos, info, &
      temp1, temp2, temp3
    !
    !     OPTION=3, the full Jacobian is not stored.  Only the product of
    !     the Jacobian transpose and Jacobian is stored.  The user provides
    !     the Jacobian one row at a time.
    !
    iopt = 3
    x(1) = 3.0E-1
    x(2) = 4.0E-1
    CALL SNLS1E(FCN3,iopt,m,n,x,fvec,tol,nprint,info,iw,wa,lwa)
    fnorm = ENORM(m,fvec)
    IF ( info==infos.AND.ABS(fnorm-fnorms)/fnorms<=tol ) THEN
      fatal = .FALSE.
      IF ( Kprint>=3 ) CALL PASS(Lun,5,1)
    ELSE
      Ipass = 0
      fatal = .TRUE.
      IF ( Kprint>=2 ) CALL PASS(Lun,5,0)
    ENDIF
    IF ( (fatal.AND.Kprint>=2).OR.Kprint>=3 ) WRITE (Lun,99007) infos, &
      fnorms, info, fnorm
    !
    !     Form JAC-transpose*JAC.
    !
    sigma = fnorm*fnorm/(m-n)
    DO i = 1, 3
      fjtj(i) = zero
    ENDDO
    iflag = 3
    DO i = 1, m
      CALL FCN3(iflag,m,n,x,fvec,fjrow,i)
      fjtj(1) = fjtj(1) + fjrow(1)**2
      fjtj(2) = fjtj(2) + fjrow(1)*fjrow(2)
      fjtj(3) = fjtj(3) + fjrow(2)**2
    ENDDO
    !
    !     Calculate the covariance matrix.
    !
    CALL SCOV(FCN3,iopt,m,n,x,fvec,fjac,ldfjac,info,wa(1),wa(n+1),wa(2*n+1),&
      wa(3*n+1))
    !
    !     Form JAC-transpose*JAC * covariance matrix (should = SIGMA*I).
    !
    temp1 = (fjtj(1)*fjac(1,1)+fjtj(2)*fjac(1,2))/sigma
    temp2 = (fjtj(1)*fjac(1,2)+fjtj(2)*fjac(2,2))/sigma
    temp3 = (fjtj(2)*fjac(1,2)+fjtj(3)*fjac(2,2))/sigma
    IF ( info==infos.AND.ABS(temp1-one)<tol2.AND.ABS(temp2)<tol2.AND.&
        ABS(temp3-one)<tol2 ) THEN
      fatal = .FALSE.
      IF ( Kprint>=3 ) CALL PASS(Lun,6,1)
    ELSE
      Ipass = 0
      fatal = .TRUE.
      IF ( Kprint>=2 ) CALL PASS(Lun,6,0)
    ENDIF
    IF ( (fatal.AND.Kprint>=2).OR.Kprint>=3 ) WRITE (Lun,99008) infos, info, &
      temp1, temp2, temp3
    !
    !     Test improper input parameters.
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
    IF ( Kprint>=3 ) WRITE (Lun,99002)
    99002 FORMAT (/' TRIGGER 2 ERROR MESSAGES',/)
    !
    lwa = 35
    iopt = 2
    x(1) = 3.0E-1
    x(2) = 4.0E-1
    CALL SNLS1E(FCN2,iopt,m,n,x,fvec,tol,nprint,info,iw,wa,lwa)
    IF ( info/=0.OR.NUMXER(nerr)/=2 ) fatal = .TRUE.
    !
    m = 0
    CALL SCOV(FCN2,iopt,m,n,x,fvec,fjac,ldfjac,info,wa(1),wa(n+1),wa(2*n+1),&
      wa(3*n+1))
    IF ( info/=0.OR.NUMXER(nerr)/=2 ) fatal = .TRUE.
    !
    !     Restore KONTRL and check to see if the tests of error detection
    !     passed.
    !
    CALL XSETF(kontrl)
    IF ( fatal ) THEN
      Ipass = 0
      IF ( Kprint>=2 ) THEN
        WRITE (Lun,99003)
        99003 FORMAT (' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
      ENDIF
    ELSEIF ( Kprint>=3 ) THEN
      WRITE (Lun,99004)
      99004 FORMAT (' ALL INCORRECT ARGUMENT TESTS PASSED')
    ENDIF
    !
    !     Print PASS/FAIL message.
    !
    IF ( Ipass==1.AND.Kprint>=2 ) WRITE (Lun,99005)
    99005 FORMAT (/' *************SNLS1E PASSED ALL TESTS*****************')
    IF ( Ipass==0.AND.Kprint>=1 ) WRITE (Lun,99006)
    99006 FORMAT (/' ************SNLS1E FAILED SOME TESTS*****************')
    !
    RETURN
    99007 FORMAT (' EXPECTED VALUE OF INFO AND RESIDUAL NORM',I5,&
      E20.9/' RETURNED VALUE OF INFO AND RESIDUAL NORM',I5,E20.9/)
    99008 FORMAT (' EXPECTED AND RETURNED VALUE OF INFO',I5,10X,&
      I5/' RETURNED PRODUCT OF (J-TRANS*J)*COVARIANCE MATRIX/SIGMA'/&
      ' (SHOULD = THE IDENTITY -- 1.0, 0.0, 1.0)'/3E20.9/)
  END SUBROUTINE SNLS1Q
  !DECK FCQX
  SUBROUTINE FCQX(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  FCQX
    !***PURPOSE  Quick check for FC.
    !***LIBRARY   SLATEC
    !***TYPE      SINGLE PRECISION (FCQX-S, DFCQX-D)
    !***KEYWORDS  QUICK CHECK
    !***AUTHOR  Hanson, R. J., (SNLA)
    !***DESCRIPTION
    !
    !   Quick check subprogram for the subroutine FC.
    !
    !   Fit discrete data by an S-shaped curve.  Evaluate the fitted curve,
    !   its first two derivatives, and probable error curve.
    !
    !   Use subprogram FC to obtain the constrained cubic B-spline
    !   representation of the curve.
    !
    !   The values of the coefficients of the B-spline as computed by FC
    !   and the values of the fitted curve as computed by BVALU in the
    !   de Boor package are tested for accuracy with the expected values.
    !   See the example program in the report sand78-1291, pp. 22-27.
    !
    !   The dimensions in the following arrays are as small as possible for
    !   the problem being solved.
    !
    !***ROUTINES CALLED  BVALU, CV, FC, IVOUT, R1MACH, SCOPY, SMOUT, SVOUT
    !***REVISION HISTORY  (YYMMDD)
    !   780801  DATE WRITTEN
    !   890718  Changed references from BVALUE to BVALU.  (WRB)
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   891004  Changed computation of XVAL.  (WRB)
    !   891004  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Restructured using IF-THEN-ELSE-ENDIF, modified tolerances
    !           to use R1MACH(4) rather than R1MACH(3) and cleaned up
    !           FORMATs.  (RWC)
    !   930214  Declarations sections added, code revised to test error
    !           returns for all values of KPRINT and code polished.  (WRB)
    !***END PROLOGUE  FCQX
    !     .. Scalar Arguments ..
    INTEGER Ipass, Kprint, Lun
    !     .. Local Scalars ..
    REAL diff, one, t, tol, xval, zero
    INTEGER kontrl, i, idigit, ii, j, l, last, mode, n, nbkpt, &
      nconst, ndata, ndeg, nerr, nord, nval
    LOGICAL fatal
    !     .. Local Arrays ..
    REAL bkpt(13), check(51), coefck(9), coeff(9), sddata(9), v(51,5), &
      w(529), work(12), xconst(11), xdata(9), yconst(11), ydata(9)
    INTEGER iw(30), nderiv(11)
    !     .. External Functions ..
    REAL BVALU, CV, R1MACH
    INTEGER NUMXER
    EXTERNAL BVALU, CV, NUMXER, R1MACH
    !     .. External Subroutines ..
    EXTERNAL FC, IVOUT, SCOPY, SMOUT, SVOUT, XGETF, XSETF
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, REAL, SQRT
    !     .. Data statements ..
    !
    DATA xdata(1), xdata(2), xdata(3), xdata(4), xdata(5), xdata(6), &
      xdata(7), xdata(8), xdata(9)/0.15E0, 0.27E0, 0.33E0, 0.40E0, &
      0.43E0, 0.47E0, 0.53E0, 0.58E0, 0.63E0/
    DATA ydata(1), ydata(2), ydata(3), ydata(4), ydata(5), ydata(6), &
      ydata(7), ydata(8), ydata(9)/0.025E0, 0.05E0, 0.13E0, 0.27E0, &
      0.37E0, 0.47E0, 0.64E0, 0.77E0, 0.87E0/
    DATA sddata(1)/0.015E0/, ndata/9/, nord/4/, nbkpt/13/, last/10/
    DATA bkpt(1), bkpt(2), bkpt(3), bkpt(4), bkpt(5), bkpt(6), bkpt(7), &
      bkpt(8), bkpt(9), bkpt(10), bkpt(11), bkpt(12), bkpt(13)&
      / - 0.6E0, -0.4E0, -0.2E0, 0.0E0, 0.2E0, 0.4E0, 0.6E0, 0.8E0, &
      0.9E0, 1.0E0, 1.1E0, 1.2E0, 1.3E0/
    !
    !     Store the data to be used to check the accuracy of the computed
    !     results.  See SAND78-1291, p.26.
    !
    DATA coefck(1), coefck(2), coefck(3), coefck(4), coefck(5), coefck(6)&
      , coefck(7), coefck(8), coefck(9)/1.186380846E-13, &
      -2.826166426E-14, -4.333929094E-15, 1.722113311E-01, &
      9.421965984E-01, 9.684708719E-01, 9.894902905E-01, &
      1.005254855E+00, 9.894902905E-01/
    DATA check(1), check(2), check(3), check(4), check(5), check(6), &
      check(7), check(8), check(9)/2.095830752E-16, 2.870188850E-05, &
      2.296151081E-04, 7.749509897E-04, 1.836920865E-03, &
      3.587736064E-03, 6.199607918E-03, 9.844747759E-03, &
      1.469536692E-02/
    DATA check(10), check(11), check(12), check(13), check(14), check(15)&
      , check(16), check(17), check(18)/2.092367672E-02, &
      2.870188851E-02, 3.824443882E-02, 4.993466504E-02, &
      6.419812979E-02, 8.146039566E-02, 1.021470253E-01, &
      1.266835812E-01, 1.554956261E-01/
    DATA check(19), check(20), check(21), check(22), check(23), check(24)&
      , check(25), check(26), check(27)/1.890087225E-01, &
      2.276484331E-01, 2.718403204E-01, 3.217163150E-01, &
      3.762338189E-01, 4.340566020E-01, 4.938484342E-01, &
      5.542730855E-01, 6.139943258E-01/
    DATA check(28), check(29), check(30), check(31), check(32), check(33)&
      , check(34), check(35), check(36)/6.716759250E-01, &
      7.259816530E-01, 7.755752797E-01, 8.191205752E-01, &
      8.556270903E-01, 8.854875002E-01, 9.094402609E-01, &
      9.282238286E-01, 9.425766596E-01/
    DATA check(37), check(38), check(39), check(40), check(41), check(42)&
      , check(43), check(44), check(45)/9.532372098E-01, &
      9.609439355E-01, 9.664352927E-01, 9.704497377E-01, &
      9.737257265E-01, 9.768786393E-01, 9.800315521E-01, &
      9.831844649E-01, 9.863373777E-01/
    DATA check(46), check(47), check(48), check(49), check(50), check(51)&
      /9.894902905E-01, 9.926011645E-01, 9.954598055E-01, &
      9.978139804E-01, 9.994114563E-01, 1.000000000E+00/
    !***FIRST EXECUTABLE STATEMENT  FCQX
    IF ( Kprint>=2 ) WRITE (Lun,99001)
    !
    99001 FORMAT ('1'/' Test FC')
    Ipass = 1
    !
    !     Broadcast SDDATA(1) value to all of SDDATA(*).
    !
    CALL SCOPY(ndata,sddata,0,sddata,1)
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
      CALL SVOUT(nbkpt,bkpt,'('' ARRAY OF KNOTS.'')',idigit)
      CALL SVOUT(ndata,xdata,'('' INDEPENDENT VARIABLE VALUES'')',idigit)
      CALL SVOUT(ndata,ydata,'('' DEPENDENT VARIABLE VALUES'')',idigit)
      CALL SVOUT(ndata,sddata,'('' DEPENDENT VARIABLE UNCERTAINTY'')',idigit)
      CALL SVOUT(nconst,xconst,'('' INDEPENDENT VARIABLE CONSTRAINT VALUES'')'&
        ,idigit)
      CALL SVOUT(nconst,yconst,'('' CONSTRAINT VALUES'')',idigit)
      CALL IVOUT(nconst,nderiv,'('' CONSTRAINT INDICATOR'')',idigit)
    ENDIF
    !
    !     Declare amount of working storage allocated to FC.
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
    CALL FC(ndata,xdata,ydata,sddata,nord,nbkpt,bkpt,nconst,xconst,yconst,&
      nderiv,mode,coeff,w,iw)
    !
    !     Check coefficients.
    !
    tol = 7.0E0*SQRT(R1MACH(4))
    diff = 0.0E0
    DO i = 1, ndata
      diff = MAX(diff,ABS(coeff(i)-coefck(i)))
    ENDDO
    IF ( diff<=tol ) THEN
      fatal = .FALSE.
      IF ( Kprint>=3 ) WRITE (Lun,99002)
      99002 FORMAT (/' FC PASSED TEST 1')
    ELSE
      Ipass = 0
      fatal = .TRUE.
      IF ( Kprint>=2 ) WRITE (Lun,99003)
      99003 FORMAT (/' FC FAILED TEST 1')
    ENDIF
    !
    IF ( (fatal.AND.Kprint>=2).OR.Kprint>=3 ) THEN
      CALL SVOUT(ndata,coefck,'(/'' PREDICTED COEFFICIENTS OF THE B-SPLINE '//&
        'FROM SAMPLE'')',idigit)
      CALL SVOUT(ndata,coeff,'(/'' COEFFICIENTS OF THE B-SPLINE COMPUTED '//&
        'BY FC'')',idigit)
    ENDIF
    !
    !     Compute value, first two derivatives and probable uncertainty.
    !
    n = nbkpt - nord
    nval = 51
    DO i = 1, nval
      !
      !       The function BVALU is in the de Boor B-spline package.
      !
      xval = REAL(i-1)/(nval-1)
      ii = 1
      DO j = 1, 3
        v(i,j+1) = BVALU(bkpt,coeff,n,nord,j-1,xval,ii,work)
      ENDDO
      v(i,1) = xval
      !
      !       The variance function CV is a companion subprogram to FC.
      !
      v(i,5) = SQRT(CV(xval,ndata,nconst,nord,nbkpt,bkpt,w))
    ENDDO
    !
    diff = 0.0E0
    DO i = 1, nval
      diff = MAX(diff,ABS(v(i,2)-check(i)))
    ENDDO
    IF ( diff<=tol ) THEN
      fatal = .FALSE.
      IF ( Kprint>=3 ) WRITE (Lun,99004)
      99004 FORMAT (/' FC (AND BVALU) PASSED TEST 2')
    ELSE
      Ipass = 0
      fatal = .TRUE.
      IF ( Kprint>=2 ) WRITE (Lun,99005)
      99005 FORMAT (/' FC (AND BVALU) FAILED TEST 2')
    ENDIF
    !
    IF ( (fatal.AND.Kprint>=2).OR.Kprint>=3 ) THEN
      !
      !       Print these values.
      !
      CALL SMOUT(nval,5,nval,v,'(16X, ''X'', 10X, ''FNCN'', 8X,'//&
        '''1ST D'', 7X, ''2ND D'', 7X, ''ERROR'')',idigit)
      WRITE (Lun,99006)
      99006 FORMAT (/' VALUES SHOULD CORRESPOND TO THOSE IN ','SAND78-1291,',&
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
    CALL FC(ndata,xdata,ydata,sddata,0,nbkpt,bkpt,nconst,xconst,yconst,nderiv,&
      mode,coeff,w,iw)
    IF ( NUMXER(nerr)/=2 ) fatal = .TRUE.
    CALL XERCLR
    !
    CALL FC(ndata,xdata,ydata,sddata,nord,0,bkpt,nconst,xconst,yconst,nderiv,&
      mode,coeff,w,iw)
    IF ( NUMXER(nerr)/=2 ) fatal = .TRUE.
    CALL XERCLR
    !
    CALL FC(-1,xdata,ydata,sddata,nord,nbkpt,bkpt,nconst,xconst,yconst,nderiv,&
      mode,coeff,w,iw)
    IF ( NUMXER(nerr)/=2 ) fatal = .TRUE.
    CALL XERCLR
    !
    mode = 0
    CALL FC(ndata,xdata,ydata,sddata,nord,nbkpt,bkpt,nconst,xconst,yconst,&
      nderiv,mode,coeff,w,iw)
    IF ( NUMXER(nerr)/=2 ) fatal = .TRUE.
    CALL XERCLR
    !
    iw(1) = 10
    CALL FC(ndata,xdata,ydata,sddata,nord,nbkpt,bkpt,nconst,xconst,yconst,&
      nderiv,mode,coeff,w,iw)
    IF ( NUMXER(nerr)/=2 ) fatal = .TRUE.
    CALL XERCLR
    !
    iw(1) = 529
    iw(2) = 2
    CALL FC(ndata,xdata,ydata,sddata,nord,nbkpt,bkpt,nconst,xconst,yconst,&
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
        99008 FORMAT (' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
      ENDIF
    ELSEIF ( Kprint>=3 ) THEN
      WRITE (Lun,99009)
      99009 FORMAT (' ALL INCORRECT ARGUMENT TESTS PASSED')
    ENDIF
    !
    !     Print PASS/FAIL message.
    !
    IF ( Ipass==1.AND.Kprint>=2 ) WRITE (Lun,99010)
    99010 FORMAT (/' *****************FC PASSED ALL TESTS*****************')
    IF ( Ipass==0.AND.Kprint>=1 ) WRITE (Lun,99011)
    99011 FORMAT (/' ****************FC FAILED SOME TESTS*****************')
    RETURN
  END SUBROUTINE FCQX
  !DECK FCN1
  SUBROUTINE FCN1(Iflag,M,N,X,Fvec,Fjac,Ldfjac)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  FCN1
    !***PURPOSE  Subsidiary to SNLS1Q.
    !***LIBRARY   SLATEC
    !***TYPE      SINGLE PRECISION (FCN1-S, DFCN1-D)
    !***KEYWORDS  QUICK CHECK
    !***AUTHOR  (UNKNOWN)
    !***DESCRIPTION
    !
    !   Subroutine which evaluates the function for test program
    !   used in quick check of SNLS1E.
    !
    !   Numerical approximation of Jacobian is used.
    !
    !***ROUTINES CALLED  (NONE)
    !***REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   930214  TYPE and declarations sections added.  (WRB)
    !***END PROLOGUE  FCN1
    !     .. Scalar Arguments ..
    REAL Fjac
    INTEGER Iflag, Ldfjac, M, N
    !     .. Array Arguments ..
    REAL Fvec(*), X(*)
    !     .. Local Scalars ..
    REAL temp, two
    INTEGER i
    !     .. Intrinsic Functions ..
    INTRINSIC EXP
    !     .. Data statements ..
    DATA two/2.0E0/
    !***FIRST EXECUTABLE STATEMENT  FCN1
    IF ( Iflag/=1 ) RETURN
    DO i = 1, M
      temp = i
      Fvec(i) = two + two*temp - EXP(temp*X(1)) - EXP(temp*X(2))
    ENDDO
  END SUBROUTINE FCN1
  !DECK FCN2
  SUBROUTINE FCN2(Iflag,M,N,X,Fvec,Fjac,Ldfjac)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  FCN2
    !***PURPOSE  Subsidiary to SNLS1Q.
    !***LIBRARY   SLATEC
    !***TYPE      SINGLE PRECISION (FCN2-S, DFCN2-D)
    !***KEYWORDS  QUICK CHECK
    !***AUTHOR  (UNKNOWN)
    !***DESCRIPTION
    !
    !   Subroutine to evaluate function and full Jacobian for test
    !   problem in quick check of SNLS1E.
    !
    !***ROUTINES CALLED  (NONE)
    !***REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   890911  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   930214  TYPE and declarations sections added and code polished.
    !           (WRB)
    !***END PROLOGUE  FCN2
    !     .. Scalar Arguments ..
    INTEGER Iflag, Ldfjac, M, N
    !     .. Array Arguments ..
    REAL Fjac(Ldfjac,*), Fvec(*), X(*)
    !     .. Local Scalars ..
    REAL temp, two
    INTEGER i
    !     .. Intrinsic Functions ..
    INTRINSIC EXP
    !     .. Data statements ..
    DATA two/2.0E0/
    !***FIRST EXECUTABLE STATEMENT  FCN2
    IF ( Iflag==0 ) RETURN
    !
    !     Should we evaluate functions or Jacobian?
    !
    IF ( Iflag==1 ) THEN
      !
      !       Evaluate functions.
      !
      DO i = 1, M
        temp = i
        Fvec(i) = two + two*temp - EXP(temp*X(1)) - EXP(temp*X(2))
      ENDDO
    ELSE
      !
      !       Evaluate Jacobian.
      !
      IF ( Iflag/=2 ) RETURN
      DO i = 1, M
        temp = i
        Fjac(i,1) = -temp*EXP(temp*X(1))
        Fjac(i,2) = -temp*EXP(temp*X(2))
      ENDDO
    ENDIF
  END SUBROUTINE FCN2
  !DECK FCN3
  SUBROUTINE FCN3(Iflag,M,N,X,Fvec,Fjrow,Nrow)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  FCN3
    !***PURPOSE  Subsidiary to SNLS1Q.
    !***LIBRARY   SLATEC
    !***TYPE      SINGLE PRECISION (FCN3-S, DFCN3-D)
    !***KEYWORDS  QUICK CHECK
    !***AUTHOR  (UNKNOWN)
    !***DESCRIPTION
    !
    !   Subroutine to evaluate the Jacobian, one row at a time, for
    !   test problem used in quick check of SNLS1E.
    !
    !***ROUTINES CALLED  (NONE)
    !***REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   890911  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   930214  TYPE and declarations sections added and code polished.
    !           (WRB)
    !***END PROLOGUE  FCN3
    !     .. Scalar Arguments ..
    INTEGER Iflag, M, N, Nrow
    !     .. Array Arguments ..
    REAL Fjrow(*), Fvec(*), X(*)
    !     .. Local Scalars ..
    REAL temp, two
    INTEGER i
    !     .. Intrinsic Functions ..
    INTRINSIC EXP
    !     .. Data statements ..
    DATA two/2.0E0/
    !***FIRST EXECUTABLE STATEMENT  FCN3
    IF ( Iflag==0 ) RETURN
    !
    !     Should we evaluate functions or Jacobian?
    !
    IF ( Iflag==1 ) THEN
      !
      !       Evaluate functions.
      !
      DO i = 1, M
        temp = i
        Fvec(i) = two + two*temp - EXP(temp*X(1)) - EXP(temp*X(2))
      ENDDO
    ELSE
      !
      !       Evaluate one row of Jacobian.
      !
      IF ( Iflag/=3 ) RETURN
      temp = Nrow
      Fjrow(1) = -temp*EXP(temp*X(1))
      Fjrow(2) = -temp*EXP(temp*X(2))
    ENDIF
  END SUBROUTINE FCN3
END MODULE TEST52_MOD
!DECK TEST52
PROGRAM TEST52
  USE TEST52_MOD
  IMPLICIT NONE
  !***BEGIN PROLOGUE  TEST52
  !***PURPOSE  Driver for testing SLATEC subprograms
  !***LIBRARY   SLATEC
  !***CATEGORY  K1, E3, K6, L
  !***TYPE      SINGLE PRECISION (TEST52-S, TEST53-D)
  !***KEYWORDS  QUICK CHECK DRIVER
  !***AUTHOR  SLATEC Common Mathematical Library Committee
  !***DESCRIPTION
  !
  ! *Usage:
  !     One input data record is required
  !         READ (LIN, '(I1)') KPRINT
  !
  ! *Arguments:
  !     KPRINT = 0  Quick checks - No printing.
  !                 Driver       - Short pass or fail message printed.
  !              1  Quick checks - No message printed for passed tests,
  !                                short message printed for failed tests.
  !                 Driver       - Short pass or fail message printed.
  !              2  Quick checks - Print short message for passed tests,
  !                                fuller information for failed tests.
  !                 Driver       - Pass or fail message printed.
  !              3  Quick checks - Print complete quick check results.
  !                 Driver       - Pass or fail message printed.
  !
  ! *Description:
  !     Driver for testing SLATEC subprograms
  !        SNLS1E   SNLS1    SCOV
  !        BVALU    CV       FC
  !        POLFIT   PCOEF    PVALUE
  !
  !***REFERENCES  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***ROUTINES CALLED  FCQX, I1MACH, PFITQX, SNLS1Q, XERMAX, XSETF, XSETUN
  !***REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  !***END PROLOGUE  TEST52
  INTEGER I1MACH
  INTEGER ipass, kprint, lin, lun, nfail
  !***FIRST EXECUTABLE STATEMENT  TEST52
  lun = I1MACH(2)
  lin = I1MACH(1)
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  CALL XERMAX(1000)
  CALL XSETUN(lun)
  IF ( kprint<=1 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  ENDIF
  !
  !     Test SNLS1E and SNLS1
  !
  CALL SNLS1Q(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test FC (also BVALU and CV)
  !
  CALL FCQX(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test POLFIT (also PCOEF and PVALUE)
  !
  CALL PFITQX(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST52 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST52 *************')
  ENDIF
  STOP
END PROGRAM TEST52
