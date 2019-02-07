!*==DNSQQK.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DNSQQK
SUBROUTINE DNSQQK(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--DNSQQK5
  !*** Start of declarations inserted by SPAG
  REAL DQFCN2, DQJAC2
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DNSQQK
  !***PURPOSE  Quick check for DNSQE and DNSQ.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (SNSQQK-S, DNSQQK-D)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !   This subroutine performs a quick check on the subroutine DNSQE
  !   (and DNSQ).
  !
  !***ROUTINES CALLED  D1MACH, DENORM, DNSQE, DQFCN2, DQJAC2, PASS
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920310  Code cleaned up and TYPE section added.  (RWC, WRB)
  !***END PROLOGUE  DNSQQK
  !     .. Scalar Arguments ..
  INTEGER Ipass, Kprint, Lun
  !     .. Local Scalars ..
  REAL(8) :: fnorm, fnorms, tol
  INTEGER icnt, info, infos, iopt, lwa, n, nprint
  !     .. Local Arrays ..
  REAL(8) :: fvec(2), wa(19), x(2)
  INTEGER itest(3)
  !     .. External Functions ..
  REAL(8) :: D1MACH, DENORM
  EXTERNAL D1MACH, DENORM
  !     .. External Subroutines ..
  EXTERNAL DNSQE, DQFCN2, DQJAC2, PASS
  !     .. Intrinsic Functions ..
  INTRINSIC SQRT
  !***FIRST EXECUTABLE STATEMENT  DNSQQK
  infos = 1
  fnorms = 0.0D0
  n = 2
  lwa = 19
  nprint = -1
  tol = SQRT(D1MACH(4))
  IF ( Kprint>=2 ) WRITE (Lun,99001)
  99001 FORMAT ('1'/'  DNSQE QUICK CHECK'/)
  !
  !     Option 1, the user provides the Jacobian.
  !
  iopt = 1
  x(1) = -1.2D0
  x(2) = 1.0D0
  CALL DNSQE(DQFCN2,DQJAC2,iopt,n,x,fvec,tol,nprint,info,wa,lwa)
  icnt = 1
  fnorm = DENORM(n,fvec)
  itest(icnt) = 0
  IF ( (info==infos).AND.(fnorm-fnorms<=tol) ) itest(icnt) = 1
  !
  IF ( Kprint/=0 ) THEN
    IF ( (Kprint>=2.AND.itest(icnt)/=1).OR.Kprint>=3 ) WRITE (Lun,99004)&
      infos, fnorms, info, fnorm
    IF ( (Kprint>=2).OR.(Kprint==1.AND.itest(icnt)/=1) )&
      CALL PASS(Lun,icnt,itest(icnt))
  ENDIF
  !
  !     Option 2, the code approximates the Jacobian.
  !
  iopt = 2
  x(1) = -1.2D0
  x(2) = 1.0D0
  CALL DNSQE(DQFCN2,DQJAC2,iopt,n,x,fvec,tol,nprint,info,wa,lwa)
  icnt = 2
  fnorm = DENORM(n,fvec)
  itest(icnt) = 0
  IF ( (info==infos).AND.(fnorm-fnorms<=tol) ) itest(icnt) = 1
  !
  IF ( Kprint/=0 ) THEN
    IF ( Kprint>=3.OR.(Kprint>=2.AND.itest(icnt)/=1) ) WRITE (Lun,99004)&
      infos, fnorms, info, fnorm
    IF ( Kprint>=2.OR.(Kprint==1.AND.itest(icnt)/=1) )&
      CALL PASS(Lun,icnt,itest(icnt))
  ENDIF
  !
  !     Test improper input parameters.
  !
  lwa = 15
  iopt = 1
  x(1) = -1.2D0
  x(2) = 1.0D0
  CALL DNSQE(DQFCN2,DQJAC2,iopt,n,x,fvec,tol,nprint,info,wa,lwa)
  icnt = 3
  itest(icnt) = 0
  IF ( info==0 ) itest(icnt) = 1
  IF ( Kprint>=2.OR.(Kprint==1.AND.itest(icnt)/=1) )&
    CALL PASS(Lun,icnt,itest(icnt))
  !
  !     Set IPASS.
  !
  Ipass = itest(1)*itest(2)*itest(3)
  IF ( Kprint>=1.AND.Ipass/=1 ) WRITE (Lun,99002)
  99002 FORMAT (/' **********WARNING -- DNSQE/DNSQ FAILED SOME TESTS****',&
    '******')
  IF ( Kprint>=2.AND.Ipass==1 ) WRITE (Lun,99003)
  99003 FORMAT (/' ----------DNSQE/DNSQ PASSED ALL TESTS----------')
  RETURN
  99004 FORMAT (' EXPECTED VALUE OF INFO AND RESIDUAL NORM',I5,&
    D20.5/' RETURNED VALUE OF INFO AND RESIDUAL NORM',I5,D20.5/)
END SUBROUTINE DNSQQK
