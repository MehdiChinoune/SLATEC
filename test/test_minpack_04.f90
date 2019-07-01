MODULE TEST53_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE
  REAL(DP) :: EPS, RP, SVEps, TOL
  INTEGER :: IERp, IERr, NORd, NORdp
  REAL(DP) :: R(11)

CONTAINS
  !** DCMPAR
  SUBROUTINE DCMPAR(Icnt,Itest)
    !> Compare values in COMMON block DCHECK for quick check
    !            routine DPFITT.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (CMPARE-S, DCMPAR-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)
    !***
    ! COMMON BLOCKS    DCHECK

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   890921  Realigned order of variables in the COMMON block.
    !           (WRB)
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   920214  Minor improvements to code for readability.  (WRB)

    !     .. Scalar Arguments ..
    INTEGER :: Icnt
    !     .. Array Arguments ..
    INTEGER :: Itest(9)
    !     .. Local Scalars ..
    REAL(DP) :: rpp, ss
    INTEGER :: ierpp, nrdp
    !     .. Local Arrays ..
    INTEGER :: itemp(4)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS
    !* FIRST EXECUTABLE STATEMENT  DCMPAR
    Icnt = Icnt + 1
    itemp(1) = 0
    itemp(2) = 0
    itemp(3) = 0
    itemp(4) = 0
    ss = SVEps - EPS
    nrdp = NORdp - NORd
    rpp = RP - R(11)
    ierpp = IERp - IERr
    IF( ABS(ss)<=TOL .OR. Icnt<=2 .OR. Icnt>=6 ) itemp(1) = 1
    IF( ABS(nrdp)==0 ) itemp(2) = 1
    IF( Icnt==2 ) itemp(2) = 1
    IF( ABS(rpp)<=TOL ) itemp(3) = 1
    IF( ABS(ierpp)==0 ) itemp(4) = 1
    !
    !     Check to see if all four tests were good.
    !     If so, set the test number equal to 1.
    !
    Itest(Icnt) = itemp(1)*itemp(2)*itemp(3)*itemp(4)
  END SUBROUTINE DCMPAR
  !** DPFITT
  SUBROUTINE DPFITT(Lun,Kprint,Ipass)
    !> Quick check for DPOLFT, DPCOEF and DP1VLU.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (PFITQX-S, DPFITT-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  D1MACH, DCMPAR, DP1VLU, DPCOEF, DPOLFT, PASS,
    !                    XERCLR, XGETF, XSETF
    !***
    ! COMMON BLOCKS    DCHECK

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   890921  Realigned order of variables in the COMMON block.
    !           (WRB)
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   900911  Test problem changed and cosmetic changes to code.  (WRB)
    !   901205  Changed usage of D1MACH(3) to D1MACH(4) and modified the
    !           FORMATs.  (RWC)
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   900911  Test problem changed and cosmetic changes to code.  (WRB)
    !   920214  Code restructured to test for all values of KPRINT and to
    !           provide more PASS/FAIL information.  (WRB)
    USE slatec, ONLY : D1MACH, DP1VLU, DPCOEF, DPOLFT, num_xer, control_xer
    USE common_mod, ONLY : PASS
    INTEGER :: kontrl
    !     .. Scalar Arguments ..
    INTEGER :: Ipass, Kprint, Lun
    !     .. Local Scalars ..
    REAL(DP) :: yfit
    INTEGER :: i, icnt, m, maxord
    !     .. Local Arrays ..
    REAL(DP) :: a(97), tc(5), w(11), x(11), y(11), yp(5)
    INTEGER :: itest(9)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, SQRT
    !* FIRST EXECUTABLE STATEMENT  DPFITT
    IF( Kprint>=2 ) WRITE (Lun,FMT=99002)
    !
    !     Initialize variables for testing passage or failure of tests
    !
    DO i = 1, 9
      itest(i) = 0
    END DO
    icnt = 0
    TOL = SQRT(D1MACH(4))
    m = 11
    DO i = 1, m
      x(i) = i - 6
      y(i) = x(i)**4
    END DO
    !
    !     Test DPOLFT
    !     Input EPS is negative - specified level
    !
    w(1) = -1._DP
    EPS = -0.01_DP
    SVEps = EPS
    maxord = 8
    NORdp = 4
    RP = 625._DP
    IERp = 1
    CALL DPOLFT(m,x,y,w,maxord,NORd,EPS,R,IERr,a)
    !
    !     See if test passed
    !
    CALL DCMPAR(icnt,itest)
    !
    !     Check for suppression of printing.
    !
    IF( Kprint/=0 ) THEN
      IF( Kprint/=1 .OR. itest(icnt)/=1 ) THEN
        WRITE (Lun,FMT=99003)
        WRITE (Lun,FMT=99004)
        IF( Kprint>2 .OR. itest(icnt)/=1 ) THEN
          WRITE (Lun,FMT=99005) SVEps, NORdp, RP, IERp
          WRITE (Lun,FMT=99006) EPS, NORd, R(11), IERr
        END IF
        !
        !     Send message indicating passage or failure of test
        !
        CALL PASS(Lun,icnt,itest(icnt))
      END IF
    END IF
    !
    !     Input EPS is negative - computed level
    !
    EPS = -1._DP
    SVEps = EPS
    CALL DPOLFT(m,x,y,w,maxord,NORd,EPS,R,IERr,a)
    !
    !     See if test passed
    !
    CALL DCMPAR(icnt,itest)
    !
    !     Check for suppression of printing.
    !
    IF( Kprint/=0 ) THEN
      IF( Kprint/=1 .OR. itest(icnt)/=1 ) THEN
        WRITE (Lun,FMT=99007)
        IF( Kprint>2 .OR. itest(icnt)/=1 ) THEN
          WRITE (Lun,FMT=99008) maxord
          WRITE (Lun,FMT=99005) SVEps, NORdp, RP, IERp
          WRITE (Lun,FMT=99006) EPS, NORd, R(11), IERr
        END IF
        !
        !     Send message indicating passage or failure of test
        !
        CALL PASS(Lun,icnt,itest(icnt))
      END IF
    END IF
    !
    !     Input EPS is zero
    !
    w(1) = -1._DP
    EPS = 0._DP
    SVEps = EPS
    NORdp = 5
    maxord = 5
    CALL DPOLFT(m,x,y,w,maxord,NORd,EPS,R,IERr,a)
    !
    !     See if test passed
    !
    CALL DCMPAR(icnt,itest)
    !
    !     Check for suppression of printing.
    !
    IF( Kprint/=0 ) THEN
      IF( Kprint/=1 .OR. itest(icnt)/=1 ) THEN
        WRITE (Lun,FMT=99009)
        IF( Kprint>2 .OR. itest(icnt)/=1 ) THEN
          WRITE (Lun,FMT=99008) maxord
          WRITE (Lun,FMT=99005) SVEps, NORdp, RP, IERp
          WRITE (Lun,FMT=99006) EPS, NORd, R(11), IERr
        END IF
        !
        !     Send message indicating passage or failure of test
        !
        CALL PASS(Lun,icnt,itest(icnt))
      END IF
    END IF
    !
    !     Input EPS is positive
    !
    IERp = 1
    NORdp = 4
    EPS = 75._DP*D1MACH(4)
    SVEps = EPS
    CALL DPOLFT(m,x,y,w,maxord,NORd,EPS,R,IERr,a)
    !
    !     See if test passed
    !
    CALL DCMPAR(icnt,itest)
    !
    !     Check for suppression of printing.
    !
    IF( Kprint/=0 ) THEN
      IF( Kprint/=1 .OR. itest(icnt)/=1 ) THEN
        WRITE (Lun,FMT=99010)
        IF( Kprint>2 .OR. itest(icnt)/=1 ) THEN
          WRITE (Lun,FMT=99008) maxord
          WRITE (Lun,FMT=99005) SVEps, NORdp, RP, IERp
          WRITE (Lun,FMT=99006) EPS, NORd, R(11), IERr
        END IF
        !
        !     Send message indicating passage or failure of test
        !
        CALL PASS(Lun,icnt,itest(icnt))
      END IF
    END IF
    !
    !     Improper input
    !
    IERp = 2
    m = -2
    !
    !     Check for suppression of printing.
    !
    kontrl = control_xer
    IF( Kprint<=2 ) THEN
      control_xer = 0
    ELSE
      control_xer = 1
    END IF
    num_xer = 0
    !
    IF( Kprint>=3 ) WRITE (Lun,99001)
    99001 FORMAT (/' Invalid input')
    CALL DPOLFT(m,x,y,w,maxord,NORd,EPS,R,IERr,a)
    !
    !     See if test passed
    !
    icnt = icnt + 1
    IF( IERr==2 ) THEN
      itest(icnt) = 1
      IF( Kprint>=3 ) WRITE (Lun,99011) 'PASSED', IERr
    ELSEIF( Kprint>=2 ) THEN
      WRITE (Lun,99011) 'FAILED', IERr
    END IF
    !
    !     Check for suppression of printing.
    !
    IF( Kprint/=0 ) THEN
      IF( Kprint/=1 .OR. itest(icnt)/=1 ) THEN
        IF( Kprint<=2 .AND. itest(icnt)==1 ) THEN
        END IF
        !
        !     Send message indicating passage or failure of test
        !
        CALL PASS(Lun,icnt,itest(icnt))
        !
        num_xer = 0
        control_xer = kontrl
      END IF
    END IF
    !
    !     MAXORD too small to meet RMS error
    !
    m = 11
    w(1) = -1._DP
    EPS = 5._DP*D1MACH(4)
    SVEps = EPS
    RP = 553._DP
    maxord = 2
    IERp = 3
    NORdp = 2
    CALL DPOLFT(m,x,y,w,maxord,NORd,EPS,R,IERr,a)
    !
    !     See if test passed
    !
    CALL DCMPAR(icnt,itest)
    !
    !     Check for suppression of printing.
    !
    IF( Kprint/=0 ) THEN
      IF( Kprint/=1 .OR. itest(icnt)/=1 ) THEN
        WRITE (Lun,FMT=99012)
        IF( Kprint>2 .OR. itest(icnt)/=1 ) THEN
          WRITE (Lun,FMT=99008) maxord
          WRITE (Lun,FMT=99005) SVEps, NORdp, RP, IERp
          WRITE (Lun,FMT=99006) EPS, NORd, R(11), IERr
        END IF
        !
        !     Send message indicating passage or failure of test
        !
        CALL PASS(Lun,icnt,itest(icnt))
      END IF
    END IF
    !
    !     MAXORD too small to meet statistical test
    !
    NORdp = 4
    IERp = 4
    RP = 625._DP
    EPS = -0.01_DP
    SVEps = EPS
    maxord = 5
    CALL DPOLFT(m,x,y,w,maxord,NORd,EPS,R,IERr,a)
    !
    !     See if test passed
    !
    CALL DCMPAR(icnt,itest)
    !
    !     Check for suppression of printing.
    !
    IF( Kprint/=0 ) THEN
      IF( Kprint/=1 .OR. itest(icnt)/=1 ) THEN
        WRITE (Lun,FMT=99013)
        IF( Kprint>2 .OR. itest(icnt)/=1 ) THEN
          WRITE (Lun,FMT=99008) maxord
          WRITE (Lun,FMT=99005) SVEps, NORdp, RP, IERp
          WRITE (Lun,FMT=99006) EPS, NORd, R(11), IERr
        END IF
        !
        !     Send message indicating passage or failure of test
        !
        CALL PASS(Lun,icnt,itest(icnt))
      END IF
    END IF
    !
    !     Test DPCOEF
    !
    maxord = 6
    EPS = 0._DP
    SVEps = EPS
    y(6) = 1._DP
    DO i = 1, m
      w(i) = 1._DP/(y(i)**2)
    END DO
    y(6) = 0._DP
    CALL DPOLFT(m,x,y,w,maxord,NORd,EPS,R,IERr,a)
    CALL DPCOEF(4,5._DP,tc,a)
    !
    !     See if test passed
    !
    icnt = icnt + 1
    IF( ABS(R(11)-tc(1))<=TOL ) itest(icnt) = 1
    !
    !     Check for suppression of printing
    !
    IF( Kprint/=0 ) THEN
      IF( Kprint/=1 .OR. itest(icnt)/=1 ) THEN
        WRITE (Lun,FMT=99014)
        IF( Kprint>2 .OR. itest(icnt)/=1 ) WRITE (Lun,FMT=99015) R(11), tc(1)
        !
        !     Send message indicating passage or failure of test
        !
        CALL PASS(Lun,icnt,itest(icnt))
      END IF
    END IF
    !
    !     Test DP1VLU
    !     Normal call
    !
    CALL DP1VLU(6,0,x(8),yfit,yp,a)
    !
    !     See if test passed
    !
    icnt = icnt + 1
    IF( ABS(R(8)-yfit)<=TOL ) itest(icnt) = 1
    !
    !     Check for suppression of printing
    !
    IF( Kprint/=0 ) THEN
      IF( Kprint/=1 .OR. itest(icnt)/=1 ) THEN
        WRITE (Lun,FMT=99016)
        WRITE (Lun,FMT=99017)
        IF( Kprint>2 .OR. itest(icnt)/=1 ) WRITE (Lun,FMT=99018) x(8), R(8), yfit
        !
        !     Send message indicating passage or failure of test
        !
        CALL PASS(Lun,icnt,itest(icnt))
      END IF
    END IF
    !
    !     Check to see if all tests passed
    !
    Ipass = 1
    DO i = 1, 9
      Ipass = Ipass*itest(i)
    END DO
    !
    IF( Ipass==1 .AND. Kprint>=3 ) WRITE (Lun,FMT=99019)
    IF( Ipass==0 .AND. Kprint>=2 ) WRITE (Lun,FMT=99020)
    RETURN
    !
    99002 FORMAT ('1'/' Test DPOLFT, DPCOEF and DP1VLU')
    99003 FORMAT (' Exercise DPOLFT')
    99004 FORMAT (' Input EPS is negative - specified significance level')
    99005 FORMAT (' Input EPS =  ',E15.8,'   correct order =  ',I3,'   R(1) = ',&
      E15.8,'   IERR = ',I1)
    99006 FORMAT (' Output EPS = ',E15.8,'   computed order = ',I3,'   R(1) = ',&
      E15.8,'   IERR = ',I1)
    99007 FORMAT (/' Input EPS is negative - computed significance level')
    99008 FORMAT (' Maximum order = ',I2)
    99009 FORMAT (/' Input EPS is zero')
    99010 FORMAT (/' Input EPS is positive')
    99011 FORMAT (' DPOLFT incorrect argument test ',A/' IERR should be 2.  It is ',&
      I4)
    99012 FORMAT (/' Cannot meet RMS error requirement')
    99013 FORMAT (/' Cannot satisfy statistical test')
    99014 FORMAT (/' Exercise DPCOEF')
    99015 FORMAT (/' For C=1.0, correct coefficient = ',E15.8,'   computed = ',&
      E15.8)
    99016 FORMAT (/' Exercise DP1VLU')
    99017 FORMAT (' Normal execution')
    99018 FORMAT (' For X = ',F5.2,'   correct P(X) = ',E15.8,&
      '    P(X) from DP1VLU = ',E15.8)
    99019 FORMAT (/' ***************DPOLFT PASSED ALL TESTS***************')
    99020 FORMAT (/' ***************DPOLFT FAILED SOME TESTS**************')
  END SUBROUTINE DPFITT
  !** DNLS1Q
  SUBROUTINE DNLS1Q(Lun,Kprint,Ipass)
    !> Quick check for DNLS1E, DNLS1 and DCOV.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (SNLS1Q-S, DNLS1Q-D)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Description:**
    !
    !   This subroutine performs a quick check on the subroutines DNLS1E
    !   (and DNLS1) and DCOV.
    !
    !***
    ! **Routines called:**  DENORM, DFCN1, DFCN2, DFCN3, DFDJC3, PASS, D1MACH,
    !                    DCOV, DNLS1E

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   890911  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   930214  Declarations sections added, code revised to test error
    !           returns for all values of KPRINT and code polished.  (WRB)
    USE slatec, ONLY : DENORM, DFDJC3, D1MACH, DCOV, DNLS1E, control_xer, &
      num_xer
    USE common_mod, ONLY : PASS
    !     .. Scalar Arguments ..
    INTEGER :: Ipass, Kprint, Lun
    !     .. Local Scalars ..
    REAL(DP) :: fnorm, fnorms, one, sigma, temp1, temp2, temp3, tol, tol2, zero
    INTEGER :: i, iflag, info, infos, iopt, kontrl, ldfjac, lwa, m, n, nprint
    LOGICAL :: fatal
    !     .. Local Arrays ..
    REAL(DP) :: fjac(10,2), fjrow(2,1), fjtj(3), fvec(10), wa(40), x(2)
    INTEGER :: iw(2)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, SQRT
    !* FIRST EXECUTABLE STATEMENT  DNLS1Q
    IF( Kprint>=2 ) WRITE (Lun,99001)
    !
    99001 FORMAT ('1'/' Test DNLS1E, DNLS1 and DCOV')
    !
    Ipass = 1
    infos = 1
    fnorms = 1.1151779E+01_DP
    m = 10
    n = 2
    lwa = 40
    ldfjac = 10
    nprint = -1
    iflag = 1
    zero = 0._DP
    one = 1._DP
    tol = MAX(SQRT(40._DP*D1MACH(4)),1.E-12_DP)
    tol2 = SQRT(tol)
    !
    !     OPTION=2, the full Jacobian is stored and the user provides the
    !     Jacobian.
    !
    iopt = 2
    x(1) = 3.0E-1_DP
    x(2) = 4.0E-1_DP
    CALL DNLS1E(DFCN2,iopt,m,n,x,fvec,tol,nprint,info,iw,wa,lwa)
    fnorm = DENORM(m,fvec)
    IF( info==infos .AND. ABS(fnorm-fnorms)/fnorms<=tol2 ) THEN
      fatal = .FALSE.
      IF( Kprint>=3 ) CALL PASS(Lun,1,1)
    ELSE
      Ipass = 0
      fatal = .TRUE.
      IF( Kprint>=2 ) CALL PASS(Lun,1,0)
    END IF
    IF( (fatal .AND. Kprint>=2) .OR. Kprint>=3 ) WRITE (Lun,99007) infos, &
      fnorms, info, fnorm
    !
    !     Form JAC-transpose*JAC.
    !
    sigma = fnorm*fnorm/(m-n)
    iflag = 2
    CALL DFCN2(iflag,m,n,x,fvec,fjac,ldfjac)
    DO i = 1, 3
      fjtj(i) = zero
    END DO
    DO i = 1, m
      fjtj(1) = fjtj(1) + fjac(i,1)**2
      fjtj(2) = fjtj(2) + fjac(i,1)*fjac(i,2)
      fjtj(3) = fjtj(3) + fjac(i,2)**2
    END DO
    !
    !     Calculate the covariance matrix.
    !
    CALL DCOV(DFCN2,iopt,m,n,x,fvec,fjac,ldfjac,info,wa(1),wa(n+1),wa(2*n+1),&
      wa(3*n+1))
    !
    !     Form JAC-transpose*JAC * covariance matrix (should = SIGMA*I).
    !
    temp1 = (fjtj(1)*fjac(1,1)+fjtj(2)*fjac(1,2))/sigma
    temp2 = (fjtj(1)*fjac(1,2)+fjtj(2)*fjac(2,2))/sigma
    temp3 = (fjtj(2)*fjac(1,2)+fjtj(3)*fjac(2,2))/sigma
    IF( info==infos .AND. ABS(temp1-one)<tol2 .AND. ABS(temp2)<tol2 .AND. &
        ABS(temp3-one)<tol2 ) THEN
      fatal = .FALSE.
      IF( Kprint>=3 ) CALL PASS(Lun,2,1)
    ELSE
      Ipass = 0
      fatal = .TRUE.
      IF( Kprint>=2 ) CALL PASS(Lun,2,0)
    END IF
    IF( (fatal .AND. Kprint>=2) .OR. Kprint>=3 ) WRITE (Lun,99008) infos, info, &
      temp1, temp2, temp3
    !
    !     OPTION=1, the full Jacobian is stored and the code approximates
    !     the Jacobian.
    !
    iopt = 1
    x(1) = 3.0E-1_DP
    x(2) = 4.0E-1_DP
    CALL DNLS1E(DFCN1,iopt,m,n,x,fvec,tol,nprint,info,iw,wa,lwa)
    fnorm = DENORM(m,fvec)
    IF( info==infos .AND. ABS(fnorm-fnorms)/fnorms<=tol2 ) THEN
      fatal = .FALSE.
      IF( Kprint>=3 ) CALL PASS(Lun,3,1)
    ELSE
      Ipass = 0
      fatal = .TRUE.
      IF( Kprint>=2 ) CALL PASS(Lun,3,0)
    END IF
    IF( (fatal .AND. Kprint>=2) .OR. Kprint>=3 ) WRITE (Lun,99007) infos, &
      fnorms, info, fnorm
    !
    !     Form JAC-transpose*JAC.
    !
    sigma = fnorm*fnorm/(m-n)
    iflag = 1
    CALL DFDJC3(DFCN1,m,n,x,fvec,fjac,ldfjac,iflag,zero,wa)
    DO i = 1, 3
      fjtj(i) = zero
    END DO
    DO i = 1, m
      fjtj(1) = fjtj(1) + fjac(i,1)**2
      fjtj(2) = fjtj(2) + fjac(i,1)*fjac(i,2)
      fjtj(3) = fjtj(3) + fjac(i,2)**2
    END DO
    !
    !     Calculate the covariance matrix.
    !
    CALL DCOV(DFCN1,iopt,m,n,x,fvec,fjac,ldfjac,info,wa(1),wa(n+1),wa(2*n+1),&
      wa(3*n+1))
    !
    !     Form JAC-transpose*JAC * covariance matrix (should = SIGMA*I).
    !
    temp1 = (fjtj(1)*fjac(1,1)+fjtj(2)*fjac(1,2))/sigma
    temp2 = (fjtj(1)*fjac(1,2)+fjtj(2)*fjac(2,2))/sigma
    temp3 = (fjtj(2)*fjac(1,2)+fjtj(3)*fjac(2,2))/sigma
    IF( info==infos .AND. ABS(temp1-one)<tol2 .AND. ABS(temp2)<tol2 .AND. &
        ABS(temp3-one)<tol2 ) THEN
      fatal = .FALSE.
      IF( Kprint>=3 ) CALL PASS(Lun,4,1)
    ELSE
      Ipass = 0
      fatal = .TRUE.
      IF( Kprint>=2 ) CALL PASS(Lun,4,0)
    END IF
    IF( (fatal .AND. Kprint>=2) .OR. Kprint>=3 ) WRITE (Lun,99008) infos, info, &
      temp1, temp2, temp3
    !
    !     OPTION=3, the full Jacobian is not stored.  Only the product of
    !     the Jacobian transpose and Jacobian is stored.  The user provides
    !     the Jacobian one row at a time.
    !
    iopt = 3
    x(1) = 3.0E-1_DP
    x(2) = 4.0E-1_DP
    CALL DNLS1E(DFCN3,iopt,m,n,x,fvec,tol,nprint,info,iw,wa,lwa)
    fnorm = DENORM(m,fvec)
    IF( info==infos .AND. ABS(fnorm-fnorms)/fnorms<=tol2 ) THEN
      fatal = .FALSE.
      IF( Kprint>=3 ) CALL PASS(Lun,5,1)
    ELSE
      Ipass = 0
      fatal = .TRUE.
      IF( Kprint>=2 ) CALL PASS(Lun,5,0)
    END IF
    IF( (fatal .AND. Kprint>=2) .OR. Kprint>=3 ) WRITE (Lun,99007) infos, &
      fnorms, info, fnorm
    !
    !     Form JAC-transpose*JAC.
    !
    sigma = fnorm*fnorm/(m-n)
    DO i = 1, 3
      fjtj(i) = zero
    END DO
    iflag = 3
    DO i = 1, m
      CALL DFCN3(iflag,m,n,x,fvec,fjrow,i)
      fjtj(1) = fjtj(1) + fjrow(1,1)**2
      fjtj(2) = fjtj(2) + fjrow(1,1)*fjrow(2,1)
      fjtj(3) = fjtj(3) + fjrow(2,1)**2
    END DO
    !
    !     Calculate the covariance matrix.
    !
    CALL DCOV(DFCN3,iopt,m,n,x,fvec,fjac,ldfjac,info,wa(1),wa(n+1),wa(2*n+1),&
      wa(3*n+1))
    !
    !     Form JAC-transpose*JAC * covariance matrix (should = SIGMA*I).
    !
    temp1 = (fjtj(1)*fjac(1,1)+fjtj(2)*fjac(1,2))/sigma
    temp2 = (fjtj(1)*fjac(1,2)+fjtj(2)*fjac(2,2))/sigma
    temp3 = (fjtj(2)*fjac(1,2)+fjtj(3)*fjac(2,2))/sigma
    IF( info==infos .AND. ABS(temp1-one)<tol2 .AND. ABS(temp2)<tol2 .AND. &
        ABS(temp3-one)<tol2 ) THEN
      fatal = .FALSE.
      IF( Kprint>=3 ) CALL PASS(Lun,6,1)
    ELSE
      Ipass = 0
      fatal = .TRUE.
      IF( Kprint>=2 ) CALL PASS(Lun,6,0)
    END IF
    IF( (fatal .AND. Kprint>=2) .OR. Kprint>=3 ) WRITE (Lun,99008) infos, info, &
      temp1, temp2, temp3
    !
    !     Test improper input parameters.
    !
    kontrl = control_xer
    IF( Kprint<=2 ) THEN
      control_xer = 0
    ELSE
      control_xer = 1
    END IF
    fatal = .FALSE.
    num_xer = 0
    !
    IF( Kprint>=3 ) WRITE (Lun,99002)
    99002 FORMAT (/' TRIGGER 2 ERROR MESSAGES',/)
    !
    lwa = 35
    iopt = 2
    x(1) = 3.0E-1_DP
    x(2) = 4.0E-1_DP
    CALL DNLS1E(DFCN2,iopt,m,n,x,fvec,tol,nprint,info,iw,wa,lwa)
    IF( info/=0 .OR. num_xer/=2 ) fatal = .TRUE.
    !
    m = 0
    CALL DCOV(DFCN2,iopt,m,n,x,fvec,fjac,ldfjac,info,wa(1),wa(n+1),wa(2*n+1),&
      wa(3*n+1))
    IF( info/=0 .OR. num_xer/=2 ) fatal = .TRUE.
    !
    !     Restore KONTRL and check to see if the tests of error detection
    !     passed.
    !
    control_xer = kontrl
    IF( fatal ) THEN
      Ipass = 0
      IF( Kprint>=2 ) THEN
        WRITE (Lun,99003)
        99003 FORMAT (' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
      END IF
    ELSEIF( Kprint>=3 ) THEN
      WRITE (Lun,99004)
      99004 FORMAT (' ALL INCORRECT ARGUMENT TESTS PASSED')
    END IF
    !
    !     Print PASS/FAIL message.
    !
    IF( Ipass==1 .AND. Kprint>=2 ) WRITE (Lun,99005)
    99005 FORMAT (/' *************DNLS1E PASSED ALL TESTS*****************')
    IF( Ipass==0 .AND. Kprint>=1 ) WRITE (Lun,99006)
    99006 FORMAT (/' ************DNLS1E FAILED SOME TESTS*****************')
    !
    RETURN
    99007 FORMAT (' EXPECTED VALUE OF INFO AND RESIDUAL NORM',I5,&
      D20.9/' RETURNED VALUE OF INFO AND RESIDUAL NORM',I5,D20.9/)
    99008 FORMAT (' EXPECTED AND RETURNED VALUE OF INFO',I5,10X,&
      I5/' RETURNED PRODUCT OF (J-TRANS*J)*COVARIANCE MATRIX/SIGMA'/&
      ' (SHOULD = THE IDENTITY -- 1.0, 0.0, 1.0)'/3D20.9/)
  END SUBROUTINE DNLS1Q
  !** DFCQX
  SUBROUTINE DFCQX(Lun,Kprint,Ipass)
    !> Quick check for DFC.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (FCQX-S, DFCQX-D)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Hanson, R. J., (SNLA)
    !***
    ! **Description:**
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
    !***
    ! **Routines called:**  D1MACH, DBVALU, DCOPY, DCV, DFC, DMOUT, DVOUT,
    !                    IVOUT

    !* REVISION HISTORY  (YYMMDD)
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
    USE slatec, ONLY : D1MACH, DBVALU, DCV, DFC, DMOUT, DVOUT, IVOUT, &
      control_xer, num_xer
    !     .. Scalar Arguments ..
    INTEGER :: Ipass, Kprint, Lun
    !     .. Local Scalars ..
    REAL(DP) :: diff, one, t, tol, xval, zero
    INTEGER :: kontrl, i, idigit, ii, j, l, mode, n, nconst, ndeg, nval
    LOGICAL :: fatal
    !     .. Local Arrays ..
    REAL(DP) :: coeff(9), v(51,5), w(529), xconst(11), yconst(11)
    INTEGER :: iw(30), nderiv(11)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, REAL, SQRT
    !     .. Data statements ..
    !
    INTEGER, PARAMETER :: ndata = 9, nord = 4, nbkpt = 13, last = 10
    REAL(DP), PARAMETER :: xdata(9) = [ 0.15_DP, 0.27_DP, 0.33_DP, 0.40_DP, 0.43_DP, &
      0.47_DP, 0.53_DP, 0.58_DP, 0.63_DP ]
    REAL(DP), PARAMETER :: ydata(9) = [ 0.025_DP, 0.05_DP, 0.13_DP, 0.27_DP, 0.37_DP, &
      0.47_DP, 0.64_DP, 0.77_DP, 0.87_DP ]
    REAL(DP), PARAMETER :: sddata(9)  = 0.015_DP
    REAL(DP), PARAMETER :: bkpt(13) = [ -0.6_DP, -0.4_DP, -0.2_DP, 0._DP, 0.2_DP, &
      0.4_DP, 0.6_DP, 0.8_DP, 0.9_DP, 1._DP, 1.1_DP, 1.2_DP, 1.3_DP ]
    !
    !     Store the data to be used to check the accuracy of the computed
    !     results.  See SAND78-1291, p.26.
    !
    REAL(DP), PARAMETER :: coefck(9) = [ 1.186380846E-13_DP, -2.826166426E-14_DP, &
      -4.333929094E-15_DP, 1.722113311E-01_DP, 9.421965984E-01_DP, 9.684708719E-01_DP, &
      9.894902905E-01_DP,  1.005254855E+00_DP, 9.894902905E-01_DP ]
    REAL(DP), PARAMETER :: check(51) = [ 2.095830752E-16_DP, 2.870188850E-05_DP, &
      2.296151081E-04_DP, 7.749509897E-04_DP, 1.836920865E-03_DP, 3.587736064E-03_DP, &
      6.199607918E-03_DP, 9.844747759E-03_DP, 1.469536692E-02_DP, 2.092367672E-02_DP, &
      2.870188851E-02_DP, 3.824443882E-02_DP, 4.993466504E-02_DP, 6.419812979E-02_DP, &
      8.146039566E-02_DP, 1.021470253E-01_DP, 1.266835812E-01_DP, 1.554956261E-01_DP, &
      1.890087225E-01_DP, 2.276484331E-01_DP, 2.718403204E-01_DP, 3.217163150E-01_DP, &
      3.762338189E-01_DP, 4.340566020E-01_DP, 4.938484342E-01_DP, 5.542730855E-01_DP, &
      6.139943258E-01_DP, 6.716759250E-01_DP, 7.259816530E-01_DP, 7.755752797E-01_DP, &
      8.191205752E-01_DP, 8.556270903E-01_DP, 8.854875002E-01_DP, 9.094402609E-01_DP, &
      9.282238286E-01_DP, 9.425766596E-01_DP, 9.532372098E-01_DP, 9.609439355E-01_DP, &
      9.664352927E-01_DP, 9.704497377E-01_DP, 9.737257265E-01_DP, 9.768786393E-01_DP, &
      9.800315521E-01_DP, 9.831844649E-01_DP, 9.863373777E-01_DP, 9.894902905E-01_DP, &
      9.926011645E-01_DP, 9.954598055E-01_DP, 9.978139804E-01_DP, 9.994114563E-01_DP, &
      1.000000000E+00_DP ]
    !* FIRST EXECUTABLE STATEMENT  DFCQX
    IF( Kprint>=2 ) WRITE (Lun,99001)
    !
    99001 FORMAT ('1'/' Test DFC')
    Ipass = 1
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
    END DO
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
    END DO
    !
    idigit = -4
    !
    IF( Kprint>=3 ) THEN
      CALL DVOUT(nbkpt,bkpt,'('' ARRAY OF KNOTS.'')',idigit)
      CALL DVOUT(ndata,xdata,'('' INDEPENDENT VARIABLE VALUES'')',idigit)
      CALL DVOUT(ndata,ydata,'('' DEPENDENT VARIABLE VALUES'')',idigit)
      CALL DVOUT(ndata,sddata,'('' DEPENDENT VARIABLE UNCERTAINTY'')',idigit)
      CALL DVOUT(nconst,xconst,'('' INDEPENDENT VARIABLE CONSTRAINT VALUES'')'&
        ,idigit)
      CALL DVOUT(nconst,yconst,'('' CONSTRAINT VALUES'')',idigit)
      CALL IVOUT(nconst,nderiv,'('' CONSTRAINT INDICATOR'')',idigit)
    END IF
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
    tol = MAX(7._DP*SQRT(D1MACH(4)),1.E-8_DP)
    diff = 0._DP
    DO i = 1, ndata
      diff = MAX(diff,ABS(coeff(i)-coefck(i)))
    END DO
    IF( diff<=tol ) THEN
      fatal = .FALSE.
      IF( Kprint>=3 ) WRITE (Lun,99002)
      99002 FORMAT (/' DFC PASSED TEST 1')
    ELSE
      Ipass = 0
      fatal = .TRUE.
      IF( Kprint>=2 ) WRITE (Lun,99003)
      99003 FORMAT (/' DFC FAILED TEST 1')
    END IF
    !
    IF( (fatal .AND. Kprint>=2) .OR. Kprint>=3 ) THEN
      CALL DVOUT(ndata,coefck,'(/'' PREDICTED COEFFICIENTS OF THE B-SPLINE '//&
        'FROM SAMPLE'')',idigit)
      CALL DVOUT(ndata,coeff,'(/'' COEFFICIENTS OF THE B-SPLINE COMPUTED '//&
        'BY DFC'')',idigit)
    END IF
    !
    !     Compute value, first two derivatives and probable uncertainty.
    !
    n = nbkpt - nord
    nval = 51
    DO i = 1, nval
      !
      !       The function DBVALU is in the de Boor B-spline package.
      !
      xval = REAL( i-1, DP)/(nval-1)
      ii = 1
      DO j = 1, 3
        v(i,j+1) = DBVALU(bkpt,coeff,n,nord,j-1,xval)
      END DO
      v(i,1) = xval
      !
      !       The variance function DCV is a companion subprogram to DFC.
      !
      v(i,5) = SQRT(DCV(xval,ndata,nconst,nord,nbkpt,bkpt,w))
    END DO
    !
    diff = 0._DP
    DO i = 1, nval
      diff = MAX(diff,ABS(v(i,2)-check(i)))
    END DO
    IF( diff<=tol ) THEN
      fatal = .FALSE.
      IF( Kprint>=3 ) WRITE (Lun,99004)
      99004 FORMAT (/' DFC (AND DBVALU) PASSED TEST 2')
    ELSE
      Ipass = 0
      fatal = .TRUE.
      IF( Kprint>=2 ) WRITE (Lun,99005)
      99005 FORMAT (/' DFC (AND DBVALU) FAILED TEST 2')
    END IF
    !
    IF( (fatal .AND. Kprint>=2) .OR. Kprint>=3 ) THEN
      !
      !       Print these values.
      !
      CALL DMOUT(nval,5,nval,v,'(16X, ''X'', 10X, ''FNCN'', 8X,'//&
        '''1ST D'', 7X, ''2ND D'', 7X, ''ERROR'')',idigit)
      WRITE (Lun,99006)
      99006 FORMAT (/' VALUES SHOULD CORRESPOND TO THOSE IN ','SAND78-1291,',&
        ' P. 26')
    END IF
    !
    !     Trigger error conditions.
    !
    kontrl = control_xer
    IF( Kprint<=2 ) THEN
      control_xer = 0
    ELSE
      control_xer = 1
    END IF
    fatal = .FALSE.
    num_xer = 0
    !
    IF( Kprint>=3 ) WRITE (Lun,99007)
    99007 FORMAT (/' TRIGGER 6 ERROR MESSAGES',/)
    !
    CALL DFC(ndata,xdata,ydata,sddata,0,nbkpt,bkpt,nconst,xconst,yconst,&
      nderiv,mode,coeff,w,iw)
    IF( num_xer/=2 ) fatal = .TRUE.
    num_xer = 0
    !
    CALL DFC(ndata,xdata,ydata,sddata,nord,0,bkpt,nconst,xconst,yconst,nderiv,&
      mode,coeff,w,iw)
    IF( num_xer/=2 ) fatal = .TRUE.
    num_xer = 0
    !
    CALL DFC(-1,xdata,ydata,sddata,nord,nbkpt,bkpt,nconst,xconst,yconst,&
      nderiv,mode,coeff,w,iw)
    IF( num_xer/=2 ) fatal = .TRUE.
    num_xer = 0
    !
    mode = 0
    CALL DFC(ndata,xdata,ydata,sddata,nord,nbkpt,bkpt,nconst,xconst,yconst,&
      nderiv,mode,coeff,w,iw)
    IF( num_xer/=2 ) fatal = .TRUE.
    num_xer = 0
    !
    iw(1) = 10
    CALL DFC(ndata,xdata,ydata,sddata,nord,nbkpt,bkpt,nconst,xconst,yconst,&
      nderiv,mode,coeff,w,iw)
    IF( num_xer/=2 ) fatal = .TRUE.
    num_xer = 0
    !
    iw(1) = 529
    iw(2) = 2
    CALL DFC(ndata,xdata,ydata,sddata,nord,nbkpt,bkpt,nconst,xconst,yconst,&
      nderiv,mode,coeff,w,iw)
    IF( num_xer/=2 ) fatal = .TRUE.
    num_xer = 0
    !
    !     Restore KONTRL and check to see if the tests of error detection
    !     passed.
    !
    control_xer = kontrl
    IF( fatal ) THEN
      Ipass = 0
      IF( Kprint>=2 ) THEN
        WRITE (Lun,99008)
        99008 FORMAT (' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
      END IF
    ELSEIF( Kprint>=3 ) THEN
      WRITE (Lun,99009)
      99009 FORMAT (' ALL INCORRECT ARGUMENT TESTS PASSED')
    END IF
    !
    !     Print PASS/FAIL message.
    !
    IF( Ipass==1 .AND. Kprint>=2 ) WRITE (Lun,99010)
    99010 FORMAT (/' ****************DFC PASSED ALL TESTS*****************')
    IF( Ipass==0 .AND. Kprint>=1 ) WRITE (Lun,99011)
    99011 FORMAT (/' ***************DFC FAILED SOME TESTS*****************')
    RETURN
  END SUBROUTINE DFCQX
  !** DFCN1
  SUBROUTINE DFCN1(Iflag,M,N,X,Fvec,Fjac,Ldfjac)
    !> Subsidiary to DNLS1Q.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (FCN1-S, DFCN1-D)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Description:**
    !
    !   Subroutine which evaluates the function for test program
    !   used in quick check of DNLS1E.
    !
    !   Numerical approximation of Jacobian is used.
    !
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   890911  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   930214  TYPE and declarations sections added.  (WRB)

    !     .. Scalar Arguments ..
    INTEGER :: Iflag, Ldfjac, M, N
    !     .. Array Arguments ..
    REAL(DP) :: Fvec(M), X(N), Fjac(:,:)
    !     .. Local Scalars ..
    INTEGER :: i
    REAL(DP) :: temp
    !     .. Intrinsic Functions ..
    INTRINSIC EXP
    !     .. Data statements ..
    REAL(DP), PARAMETER :: two = 2._DP
    !* FIRST EXECUTABLE STATEMENT  DFCN1
    IF( Iflag/=1 ) RETURN
    DO i = 1, M
      temp = i
      Fvec(i) = two + two*temp - EXP(temp*X(1)) - EXP(temp*X(2))
    END DO
  END SUBROUTINE DFCN1
  !** DFCN2
  SUBROUTINE DFCN2(Iflag,M,N,X,Fvec,Fjac,Ldfjac)
    !> Subsidiary to DNLS1Q.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (FCN2-S, DFCN2-D)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Description:**
    !
    !   Subroutine to evaluate function and full Jacobian for test
    !   problem in quick check of DNLS1E.
    !
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   890911  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   930214  TYPE and declarations sections added and code polished.
    !           (WRB)

    !     .. Scalar Arguments ..
    INTEGER :: Iflag, Ldfjac, M, N
    !     .. Array Arguments ..
    REAL(DP) :: Fjac(:,:), Fvec(M), X(N)
    !     .. Local Scalars ..
    REAL(DP) :: temp
    INTEGER :: i
    !     .. Intrinsic Functions ..
    INTRINSIC EXP
    !     .. Data statements ..
    REAL(DP), PARAMETER :: two = 2._DP
    !* FIRST EXECUTABLE STATEMENT  DFCN2
    IF( Iflag==0 ) RETURN
    !
    !     Should we evaluate function or Jacobian?
    !
    IF( Iflag==1 ) THEN
      !
      !       Evaluate functions.
      !
      DO i = 1, M
        temp = i
        Fvec(i) = two + two*temp - EXP(temp*X(1)) - EXP(temp*X(2))
      END DO
    ELSE
      !
      !       Evaluate Jacobian.
      !
      IF( Iflag/=2 ) RETURN
      DO i = 1, M
        temp = i
        Fjac(i,1) = -temp*EXP(temp*X(1))
        Fjac(i,2) = -temp*EXP(temp*X(2))
      END DO
    END IF
  END SUBROUTINE DFCN2
  !** DFCN3
  SUBROUTINE DFCN3(Iflag,M,N,X,Fvec,Fjrow,Nrow)
    !> Subsidiary to DNLS1Q.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (FCN3-S, DFCN3-D)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Description:**
    !
    !   Subroutine to evaluate the Jacobian, one row at a time, for
    !   test problem used in quick check of DNLS1E.
    !
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   890911  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   930214  TYPE and declarations sections added and code polished.
    !           (WRB)

    !     .. Scalar Arguments ..
    INTEGER :: Iflag, M, N, Nrow
    !     .. Array Arguments ..
    REAL(DP) :: Fjrow(:,:), Fvec(M), X(N)
    !     .. Local Scalars ..
    REAL(DP) :: temp
    INTEGER :: i
    !     .. Intrinsic Functions ..
    INTRINSIC EXP
    !     .. Data statements ..
    REAL(DP), PARAMETER :: two = 2._DP
    !* FIRST EXECUTABLE STATEMENT  DFCN3
    IF( Iflag==0 ) RETURN
    !
    !     Should we evaluate functions or Jacobian?
    !
    IF( Iflag==1 ) THEN
      !
      !       Evaluate functions.
      !
      DO i = 1, M
        temp = i
        Fvec(i) = two + two*temp - EXP(temp*X(1)) - EXP(temp*X(2))
      END DO
    ELSE
      !
      !       Evaluate one row of Jacobian.
      !
      IF( Iflag/=3 ) RETURN
      temp = Nrow
      Fjrow(1,1) = -temp*EXP(temp*X(1))
      Fjrow(2,1) = -temp*EXP(temp*X(2))
    END IF
  END SUBROUTINE DFCN3
END MODULE TEST53_MOD
!** TEST53
PROGRAM TEST53
  USE TEST53_MOD, ONLY : DFCQX, DNLS1Q, DPFITT
  USE slatec, ONLY : I1MACH, control_xer, max_xer
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !> Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  K1, E3, K6, L
  !***
  ! **Type:**      DOUBLE PRECISION (TEST52-S, TEST53-D)
  !***
  ! **Keywords:**  QUICK CHECK DRIVER
  !***
  ! **Author:**  SLATEC Common Mathematical Library Committee
  !***
  ! **Description:**
  !
  !- Usage:
  !     One input data record is required
  !         READ (LIN, '(I1)') KPRINT
  !
  !- Arguments:
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
  !- Description:
  !     Driver for testing SLATEC subprograms
  !        DNLS1E   DNLS1    DCOV
  !        DBVALU   DCV      DFC
  !        DPOLFT   DPCOEF   DP1VLU
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  DFCQX, DNLS1Q, DPFITT, I1MACH, XERMAX, XSETF,
  !                    XSETUN

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  INTEGER :: ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST53
  lun = I1MACH(2)
  lin = I1MACH(1)
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  max_xer = 1000
  IF( kprint<=1 ) THEN
    control_xer = 0
  ELSE
    control_xer = 1
  END IF
  !
  !     Test DNLS1E and DNLS1
  !
  CALL DNLS1Q(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Test DFC (also DBVALU and DCV)
  !
  CALL DFCQX(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Test DPOLFT (also DPCOEF and DPLVlU)
  !
  CALL DPFITT(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST53 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST53 *************')
  END IF
  STOP
END PROGRAM TEST53
