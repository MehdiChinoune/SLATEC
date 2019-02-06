!*==DSOSQX.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DSOSQX
      SUBROUTINE DSOSQX(Lun,Kprint,Ipass)
      IMPLICIT NONE
!*--DSOSQX5
!***BEGIN PROLOGUE  DSOSQX
!***PURPOSE  Quick check for DSOS.
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (SOSNQX-S, DSOSQX-D)
!***KEYWORDS  QUICK CHECK
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!   This subroutine performs a quick check on the subroutine DSOS.
!
!***ROUTINES CALLED  D1MACH, DNRM2, DSOS, DSOSFN, PASS
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890618  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920310  Code cleaned up and TYPE section added.  (RWC, WRB)
!***END PROLOGUE  DSOSQX
!     .. Scalar Arguments ..
      INTEGER Ipass , Kprint , Lun
!     .. Local Scalars ..
      DOUBLE PRECISION aer , fnorm , fnorms , rer , tolf
      INTEGER icnt , iflag , iflags , liw , lwa , n
!     .. Local Arrays ..
      DOUBLE PRECISION fvec(2) , wa(17) , x(2)
      INTEGER itest(2) , iw(6)
!     .. External Functions ..
      DOUBLE PRECISION D1MACH , DNRM2 , DSOSFN
      EXTERNAL D1MACH , DNRM2 , DSOSFN
!     .. External Subroutines ..
      EXTERNAL DSOS , PASS
!     .. Intrinsic Functions ..
      INTRINSIC SQRT
!***FIRST EXECUTABLE STATEMENT  DSOSQX
      iflags = 3
      fnorms = 0.0D0
      n = 2
      lwa = 17
      liw = 6
      tolf = SQRT(D1MACH(4))
      rer = SQRT(D1MACH(4))
      aer = 0.0D0
      IF ( Kprint>=2 ) WRITE (Lun,99001)
99001 FORMAT ('1'/'  DSOS QUICK CHECK'/)
!
!     Test the code with proper input values.
!
      iflag = 0
      x(1) = -1.2D0
      x(2) = 1.0D0
      CALL DSOS(DSOSFN,n,x,rer,aer,tolf,iflag,wa,lwa,iw,liw)
      icnt = 1
      fvec(1) = DSOSFN(x,1)
      fvec(2) = DSOSFN(x,2)
      fnorm = DNRM2(n,fvec,1)
      itest(icnt) = 0
      IF ( iflag<=iflags.AND.fnorm-fnorms<=rer ) itest(icnt) = 1
!
      IF ( Kprint/=0 ) THEN
        IF ( Kprint>=3.OR.(Kprint>=2.AND.itest(icnt)/=1) ) WRITE (Lun,99002)
     &       iflags , fnorms , iflag , fnorm
99002   FORMAT (' EXPECTED VALUE OF IFLAG AND RESIDUAL NORM',I5,
     &          D20.5/' RETURNED VALUE OF IFLAG AND RESIDUAL NORM',I5,D20.5/)
        IF ( Kprint>=2.OR.(Kprint==1.AND.itest(icnt)/=1) )
     &       CALL PASS(Lun,icnt,itest(icnt))
      ENDIF
!
!     Test improper input parameters.
!
      lwa = 15
      iflag = 0
      x(1) = -1.2D0
      x(2) = 1.0D0
      CALL DSOS(DSOSFN,n,x,rer,aer,tolf,iflag,wa,lwa,iw,liw)
      icnt = 2
      itest(icnt) = 0
      IF ( iflag==9 ) itest(icnt) = 1
      IF ( Kprint>=2.OR.(Kprint==1.AND.itest(icnt)/=1) )
     &     CALL PASS(Lun,icnt,itest(icnt))
!
!     Set IPASS.
!
      Ipass = itest(1)*itest(2)
      IF ( Kprint>=1.AND.Ipass/=1 ) WRITE (Lun,99003)
99003 FORMAT (/' **********WARNING -- DSOS FAILED SOME TESTS**********')
      IF ( Kprint>=2.AND.Ipass==1 ) WRITE (Lun,99004)
99004 FORMAT (/' ----------DSOS PASSED ALL TESTS----------')
      RETURN
      END SUBROUTINE DSOSQX
