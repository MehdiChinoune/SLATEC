!*==DNLS1Q.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DNLS1Q
SUBROUTINE DNLS1Q(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--DNLS1Q5
  !*** Start of declarations inserted by SPAG
  REAL DFCN1
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DNLS1Q
  !***PURPOSE  Quick check for DNLS1E, DNLS1 and DCOV.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (SNLS1Q-S, DNLS1Q-D)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !   This subroutine performs a quick check on the subroutines DNLS1E
  !   (and DNLS1) and DCOV.
  !
  !***ROUTINES CALLED  DENORM, DFCN1, DFCN2, DFCN3, DFDJC3, PASS, D1MACH,
  !                    DCOV, DNLS1E
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   930214  Declarations sections added, code revised to test error
  !           returns for all values of KPRINT and code polished.  (WRB)
  !***END PROLOGUE  DNLS1Q
  !     .. Scalar Arguments ..
  INTEGER Ipass, Kprint, Lun
  !     .. Local Scalars ..
  REAL(8) :: fnorm, fnorms, one, sigma, temp1, temp2, temp3, &
    tol, tol2, zero
  INTEGER i, iflag, info, infos, iopt, kontrl, ldfjac, lwa, m, n, &
    nerr, nprint
  LOGICAL fatal
  !     .. Local Arrays ..
  REAL(8) :: fjac(10,2), fjrow(2), fjtj(3), fvec(10), wa(40), &
    x(2)
  INTEGER iw(2)
  !     .. External Functions ..
  REAL(8) :: D1MACH, DENORM
  INTEGER NUMXER
  EXTERNAL D1MACH, DENORM, NUMXER
  !     .. External Subroutines ..
  EXTERNAL DFCN1, DFCN2, DFCN3, DFDJC3, PASS, DCOV, DNLS1E, XGETF, &
    XSETF
  !     .. Intrinsic Functions ..
  INTRINSIC ABS, SQRT
  !***FIRST EXECUTABLE STATEMENT  DNLS1Q
  IF ( Kprint>=2 ) WRITE (Lun,99001)
  !
  99001 FORMAT ('1'/' Test DNLS1E, DNLS1 and DCOV')
  !
  Ipass = 1
  infos = 1
  fnorms = 1.1151779D+01
  m = 10
  n = 2
  lwa = 40
  ldfjac = 10
  nprint = -1
  iflag = 1
  zero = 0.0D0
  one = 1.0D0
  tol = MAX(SQRT(40.0D0*D1MACH(4)),1.0D-12)
  tol2 = SQRT(tol)
  !
  !     OPTION=2, the full Jacobian is stored and the user provides the
  !     Jacobian.
  !
  iopt = 2
  x(1) = 3.0D-1
  x(2) = 4.0D-1
  CALL DNLS1E(DFCN2,iopt,m,n,x,fvec,tol,nprint,info,iw,wa,lwa)
  fnorm = DENORM(m,fvec)
  IF ( info==infos.AND.ABS(fnorm-fnorms)/fnorms<=tol2 ) THEN
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
  CALL DFCN2(iflag,m,n,x,fvec,fjac,ldfjac)
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
  CALL DCOV(DFCN2,iopt,m,n,x,fvec,fjac,ldfjac,info,wa(1),wa(n+1),wa(2*n+1),&
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
  x(1) = 3.0D-1
  x(2) = 4.0D-1
  CALL DNLS1E(DFCN1,iopt,m,n,x,fvec,tol,nprint,info,iw,wa,lwa)
  fnorm = DENORM(m,fvec)
  IF ( info==infos.AND.ABS(fnorm-fnorms)/fnorms<=tol2 ) THEN
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
  CALL DFDJC3(DFCN1,m,n,x,fvec,fjac,ldfjac,iflag,zero,wa)
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
  CALL DCOV(DFCN1,iopt,m,n,x,fvec,fjac,ldfjac,info,wa(1),wa(n+1),wa(2*n+1),&
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
  x(1) = 3.0D-1
  x(2) = 4.0D-1
  CALL DNLS1E(DFCN3,iopt,m,n,x,fvec,tol,nprint,info,iw,wa,lwa)
  fnorm = DENORM(m,fvec)
  IF ( info==infos.AND.ABS(fnorm-fnorms)/fnorms<=tol2 ) THEN
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
    CALL DFCN3(iflag,m,n,x,fvec,fjrow,i)
    fjtj(1) = fjtj(1) + fjrow(1)**2
    fjtj(2) = fjtj(2) + fjrow(1)*fjrow(2)
    fjtj(3) = fjtj(3) + fjrow(2)**2
  ENDDO
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
  x(1) = 3.0D-1
  x(2) = 4.0D-1
  CALL DNLS1E(DFCN2,iopt,m,n,x,fvec,tol,nprint,info,iw,wa,lwa)
  IF ( info/=0.OR.NUMXER(nerr)/=2 ) fatal = .TRUE.
  !
  m = 0
  CALL DCOV(DFCN2,iopt,m,n,x,fvec,fjac,ldfjac,info,wa(1),wa(n+1),wa(2*n+1),&
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
      99003     FORMAT (' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
    ENDIF
  ELSEIF ( Kprint>=3 ) THEN
    WRITE (Lun,99004)
    99004   FORMAT (' ALL INCORRECT ARGUMENT TESTS PASSED')
  ENDIF
  !
  !     Print PASS/FAIL message.
  !
  IF ( Ipass==1.AND.Kprint>=2 ) WRITE (Lun,99005)
  99005 FORMAT (/' *************DNLS1E PASSED ALL TESTS*****************')
  IF ( Ipass==0.AND.Kprint>=1 ) WRITE (Lun,99006)
  99006 FORMAT (/' ************DNLS1E FAILED SOME TESTS*****************')
  !
  RETURN
  99007 FORMAT (' EXPECTED VALUE OF INFO AND RESIDUAL NORM',I5,&
    D20.9/' RETURNED VALUE OF INFO AND RESIDUAL NORM',I5,D20.9/)
  99008 FORMAT (' EXPECTED AND RETURNED VALUE OF INFO',I5,10X,&
    I5/' RETURNED PRODUCT OF (J-TRANS*J)*COVARIANCE MATRIX/SIGMA'/&
    ' (SHOULD = THE IDENTITY -- 1.0, 0.0, 1.0)'/3D20.9/)
END SUBROUTINE DNLS1Q
