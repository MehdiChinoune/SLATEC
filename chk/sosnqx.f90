!*==SOSNQX.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK SOSNQX
SUBROUTINE SOSNQX(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--SOSNQX5
  !***BEGIN PROLOGUE  SOSNQX
  !***PURPOSE  Quick check for SOS.
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (SOSNQX-S, DSOSQX-D)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Watts, H. A., (SNLA)
  !***DESCRIPTION
  !
  !   This subroutine performs a quick check on the subroutine SOS.
  !
  !***ROUTINES CALLED  PASS, R1MACH, SNRM2, SOS, SOSFNC
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920310  Code cleaned up and TYPE section added.  (RWC, WRB)
  !***END PROLOGUE  SOSNQX
  !     .. Scalar Arguments ..
  INTEGER Ipass, Kprint, Lun
  !     .. Local Scalars ..
  REAL aer, fnorm, fnorms, rer, tolf
  INTEGER icnt, iflag, iflags, liw, lwa, n
  !     .. Local Arrays ..
  REAL fvec(2), wa(17), x(2)
  INTEGER itest(2), iw(6)
  !     .. External Functions ..
  REAL R1MACH, SNRM2, SOSFNC
  EXTERNAL R1MACH, SNRM2, SOSFNC
  !     .. External Subroutines ..
  EXTERNAL PASS, SOS
  !     .. Intrinsic Functions ..
  INTRINSIC SQRT
  !***FIRST EXECUTABLE STATEMENT  SOSNQX
  iflags = 3
  fnorms = 0.0E0
  n = 2
  lwa = 17
  liw = 6
  tolf = SQRT(R1MACH(4))
  rer = SQRT(R1MACH(4))
  aer = 0.0E0
  IF ( Kprint>=2 ) WRITE (Lun,99001)
  99001 FORMAT ('1'/'  SOS QUICK CHECK'/)
  !
  !     Test the code with proper input values.
  !
  iflag = 0
  x(1) = -1.2E0
  x(2) = 1.0E0
  CALL SOS(SOSFNC,n,x,rer,aer,tolf,iflag,wa,lwa,iw,liw)
  icnt = 1
  fvec(1) = SOSFNC(x,1)
  fvec(2) = SOSFNC(x,2)
  fnorm = SNRM2(n,fvec,1)
  itest(icnt) = 0
  IF ( iflag<=iflags.AND.fnorm-fnorms<=rer ) itest(icnt) = 1
  !
  IF ( Kprint/=0 ) THEN
    IF ( Kprint>=3.OR.(Kprint>=2.AND.itest(icnt)/=1) ) WRITE (Lun,99002)&
      iflags, fnorms, iflag, fnorm
    99002   FORMAT (' EXPECTED VALUE OF IFLAG AND RESIDUAL NORM',I5,&
      E20.5/' RETURNED VALUE OF IFLAG AND RESIDUAL NORM',I5,E20.5/)
    IF ( Kprint>=2.OR.(Kprint==1.AND.itest(icnt)/=1) )&
      CALL PASS(Lun,icnt,itest(icnt))
  ENDIF
  !
  !     Test improper input parameters.
  !
  lwa = 15
  iflag = 0
  x(1) = -1.2E0
  x(2) = 1.0E0
  CALL SOS(SOSFNC,n,x,rer,aer,tolf,iflag,wa,lwa,iw,liw)
  icnt = 2
  itest(icnt) = 0
  IF ( iflag==9 ) itest(icnt) = 1
  IF ( Kprint>=2.OR.(Kprint==1.AND.itest(icnt)/=1) )&
    CALL PASS(Lun,icnt,itest(icnt))
  !
  !     Set IPASS.
  !
  Ipass = itest(1)*itest(2)
  IF ( Kprint>=1.AND.Ipass/=1 ) WRITE (Lun,99003)
  99003 FORMAT (/' **********WARNING -- SOS FAILED SOME TESTS**********')
  IF ( Kprint>=2.AND.Ipass==1 ) WRITE (Lun,99004)
  99004 FORMAT (/' ----------SOS PASSED ALL TESTS----------')
  RETURN
END SUBROUTINE SOSNQX
