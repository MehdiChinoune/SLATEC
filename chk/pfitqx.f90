!*==PFITQX.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK PFITQX
SUBROUTINE PFITQX(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--PFITQX5
  !*** Start of declarations inserted by SPAG
  INTEGER kontrl
  !*** End of declarations inserted by SPAG
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
  !     .. Scalar Arguments ..
  INTEGER Ipass , Kprint , Lun
  !     .. Scalars in Common ..
  REAL EPS , RP , SVEps , TOL
  INTEGER IERp , IERr , NORd , NORdp
  !     .. Arrays in Common ..
  REAL R(11)
  !     .. Local Scalars ..
  REAL yfit
  INTEGER i , icnt , m , maxord
  !     .. Local Arrays ..
  REAL a(97) , tc(5) , w(11) , x(11) , y(11) , yp(5)
  INTEGER itest(9)
  !     .. External Functions ..
  REAL R1MACH
  EXTERNAL R1MACH
  !     .. External Subroutines ..
  EXTERNAL CMPARE , PASS , PCOEF , POLFIT , PVALUE
  !     .. Intrinsic Functions ..
  INTRINSIC ABS , SQRT
  !     .. Common blocks ..
  COMMON /CHECK / EPS , R , RP , SVEps , TOL , NORdp , NORd , IERp , IERr
  !***FIRST EXECUTABLE STATEMENT  PFITQX
  IF ( Kprint>=2 ) WRITE (Lun,FMT=99002)
  !
  !     Initialize variables for testing passage or failure of tests
  !
  DO i = 1 , 9
    itest(i) = 0
  ENDDO
  icnt = 0
  TOL = SQRT(R1MACH(4))
  m = 11
  DO i = 1 , m
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
        WRITE (Lun,FMT=99005) SVEps , NORdp , RP , IERp
        WRITE (Lun,FMT=99006) EPS , NORd , R(11) , IERr
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
        WRITE (Lun,FMT=99005) SVEps , NORdp , RP , IERp
        WRITE (Lun,FMT=99006) EPS , NORd , R(11) , IERr
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
        WRITE (Lun,FMT=99005) SVEps , NORdp , RP , IERp
        WRITE (Lun,FMT=99006) EPS , NORd , R(11) , IERr
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
        WRITE (Lun,FMT=99005) SVEps , NORdp , RP , IERp
        WRITE (Lun,FMT=99006) EPS , NORd , R(11) , IERr
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
    IF ( Kprint>=3 ) WRITE (Lun,99011) 'PASSED' , IERr
  ELSEIF ( Kprint>=2 ) THEN
    WRITE (Lun,99011) 'FAILED' , IERr
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
        WRITE (Lun,FMT=99005) SVEps , NORdp , RP , IERp
        WRITE (Lun,FMT=99006) EPS , NORd , R(11) , IERr
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
        WRITE (Lun,FMT=99005) SVEps , NORdp , RP , IERp
        WRITE (Lun,FMT=99006) EPS , NORd , R(11) , IERr
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
  DO i = 1 , m
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
      IF ( Kprint>2.OR.itest(icnt)/=1 ) WRITE (Lun,FMT=99015) R(11) , tc(1)
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
      IF ( Kprint>2.OR.itest(icnt)/=1 ) WRITE (Lun,FMT=99018) x(8) , R(8) , &
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
  DO i = 1 , 9
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
