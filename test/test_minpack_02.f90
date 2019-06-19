MODULE TEST36_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  !** DNSQQK
  SUBROUTINE DNSQQK(Lun,Kprint,Ipass)
    !> Quick check for DNSQE and DNSQ.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (SNSQQK-S, DNSQQK-D)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Description:**
    !
    !   This subroutine performs a quick check on the subroutine DNSQE
    !   (and DNSQ).
    !
    !***
    ! **Routines called:**  D1MACH, DENORM, DNSQE, DQFCN2, DQJAC2, PASS

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   920310  Code cleaned up and TYPE section added.  (RWC, WRB)
    USE slatec, ONLY : D1MACH, DENORM, DNSQE
    USE common_mod, ONLY : PASS
    !     .. Scalar Arguments ..
    INTEGER :: Ipass, Kprint, Lun
    !     .. Local Scalars ..
    REAL(DP) :: fnorm, fnorms, tol
    INTEGER :: icnt, info, infos, iopt, lwa, n, nprint
    !     .. Local Arrays ..
    REAL(DP) :: fvec(2), wa(19), x(2)
    INTEGER :: itest(3)
    !     .. Intrinsic Functions ..
    INTRINSIC SQRT
    !* FIRST EXECUTABLE STATEMENT  DNSQQK
    infos = 1
    fnorms = 0.0D0
    n = 2
    lwa = 19
    nprint = -1
    tol = SQRT(D1MACH(4))
    IF( Kprint>=2 ) WRITE (Lun,99001)
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
    IF( (info==infos) .AND. (fnorm-fnorms<=tol) ) itest(icnt) = 1
    !
    IF( Kprint/=0 ) THEN
      IF( (Kprint>=2 .AND. itest(icnt)/=1) .OR. Kprint>=3 ) WRITE (Lun,99004)&
        infos, fnorms, info, fnorm
      IF( (Kprint>=2) .OR. (Kprint==1 .AND. itest(icnt)/=1) )&
        CALL PASS(Lun,icnt,itest(icnt))
    END IF
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
    IF( (info==infos) .AND. (fnorm-fnorms<=tol) ) itest(icnt) = 1
    !
    IF( Kprint/=0 ) THEN
      IF( Kprint>=3 .OR. (Kprint>=2 .AND. itest(icnt)/=1) ) WRITE (Lun,99004)&
        infos, fnorms, info, fnorm
      IF( Kprint>=2 .OR. (Kprint==1 .AND. itest(icnt)/=1) )&
        CALL PASS(Lun,icnt,itest(icnt))
    END IF
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
    IF( info==0 ) itest(icnt) = 1
    IF( Kprint>=2 .OR. (Kprint==1 .AND. itest(icnt)/=1) )&
      CALL PASS(Lun,icnt,itest(icnt))
    !
    !     Set IPASS.
    !
    Ipass = itest(1)*itest(2)*itest(3)
    IF( Kprint>=1 .AND. Ipass/=1 ) WRITE (Lun,99002)
    99002 FORMAT (/' **********WARNING -- DNSQE/DNSQ FAILED SOME TESTS****',&
      '******')
    IF( Kprint>=2 .AND. Ipass==1 ) WRITE (Lun,99003)
    99003 FORMAT (/' ----------DNSQE/DNSQ PASSED ALL TESTS----------')
    RETURN
    99004 FORMAT (' EXPECTED VALUE OF INFO AND RESIDUAL NORM',I5,&
      D20.5/' RETURNED VALUE OF INFO AND RESIDUAL NORM',I5,D20.5/)
  END SUBROUTINE DNSQQK
  !** DSOSFN
  REAL(DP) FUNCTION DSOSFN(X,K)
    !> Function evaluator for DSOS quick check.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Watts, H. A., (SNLA)
    !***
    ! **Description:**
    !
    !     FUNCTION WHICH EVALUATES THE FUNCTIONS, ONE AT A TIME,
    !     FOR TEST PROGRAM USED IN QUICK CHECK OF DSOS.
    !
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   801001  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    INTEGER :: K
    REAL(DP) :: X(:)
    !* FIRST EXECUTABLE STATEMENT  DSOSFN
    IF( K==1 ) THEN
      DSOSFN = 1.0D0 - X(1)
    ELSEIF( K==2 ) THEN
      DSOSFN = 1.0D1*(X(2)-X(1)**2)
    ELSE
      DSOSFN = 0.D0
    END IF
  END FUNCTION DSOSFN
  !** DSOSQX
  SUBROUTINE DSOSQX(Lun,Kprint,Ipass)
    !> Quick check for DSOS.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (SOSNQX-S, DSOSQX-D)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Watts, H. A., (SNLA)
    !***
    ! **Description:**
    !
    !   This subroutine performs a quick check on the subroutine DSOS.
    !
    !***
    ! **Routines called:**  D1MACH, DNRM2, DSOS, DSOSFN, PASS

    !* REVISION HISTORY  (YYMMDD)
    !   801001  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   920310  Code cleaned up and TYPE section added.  (RWC, WRB)
    USE slatec, ONLY : D1MACH, DSOS
    USE common_mod, ONLY : PASS
    !     .. Scalar Arguments ..
    INTEGER :: Ipass, Kprint, Lun
    !     .. Local Scalars ..
    REAL(DP) :: aer, fnorm, fnorms, rer, tolf
    INTEGER :: icnt, iflag, iflags, liw, lwa, n
    !     .. Local Arrays ..
    REAL(DP) :: fvec(2), wa(17), x(2)
    INTEGER :: itest(2), iw(6)
    !     .. Intrinsic Functions ..
    INTRINSIC SQRT
    !* FIRST EXECUTABLE STATEMENT  DSOSQX
    iflags = 3
    fnorms = 0.0D0
    n = 2
    lwa = 17
    liw = 6
    tolf = SQRT(D1MACH(4))
    rer = SQRT(D1MACH(4))
    aer = 0.0D0
    IF( Kprint>=2 ) WRITE (Lun,99001)
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
    fnorm = NORM2(fvec(1:n))
    itest(icnt) = 0
    IF( iflag<=iflags .AND. fnorm-fnorms<=rer ) itest(icnt) = 1
    !
    IF( Kprint/=0 ) THEN
      IF( Kprint>=3 .OR. (Kprint>=2 .AND. itest(icnt)/=1) ) WRITE (Lun,99002)&
        iflags, fnorms, iflag, fnorm
      99002 FORMAT (' EXPECTED VALUE OF IFLAG AND RESIDUAL NORM',I5,&
        D20.5/' RETURNED VALUE OF IFLAG AND RESIDUAL NORM',I5,D20.5/)
      IF( Kprint>=2 .OR. (Kprint==1 .AND. itest(icnt)/=1) )&
        CALL PASS(Lun,icnt,itest(icnt))
    END IF
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
    IF( iflag==9 ) itest(icnt) = 1
    IF( Kprint>=2 .OR. (Kprint==1 .AND. itest(icnt)/=1) )&
      CALL PASS(Lun,icnt,itest(icnt))
    !
    !     Set IPASS.
    !
    Ipass = itest(1)*itest(2)
    IF( Kprint>=1 .AND. Ipass/=1 ) WRITE (Lun,99003)
    99003 FORMAT (/' **********WARNING -- DSOS FAILED SOME TESTS**********')
    IF( Kprint>=2 .AND. Ipass==1 ) WRITE (Lun,99004)
    99004 FORMAT (/' ----------DSOS PASSED ALL TESTS----------')
    RETURN
  END SUBROUTINE DSOSQX
  !** DQFCN2
  SUBROUTINE DQFCN2(N,X,Fvec,Iflag)
    !> Evaluate function used in DNSQE.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (SQFCN2-S, DQFCN2-D)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Description:**
    !
    !   Subroutine which evaluates the function for test program
    !   used in quick check of DNSQE.
    !
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   890831  Modified array declarations.  (WRB)
    !   890831  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   930214  TYPE and declarations sections added.  (WRB)

    !     .. Scalar Arguments ..
    INTEGER :: Iflag, N
    !     .. Array Arguments ..
    REAL(DP) :: Fvec(N), X(N)
    !* FIRST EXECUTABLE STATEMENT  DQFCN2
    Fvec(1) = 1.0D0 - X(1)
    Fvec(2) = 10.0D0*(X(2)-X(1)**2)
  END SUBROUTINE DQFCN2
  !** DQJAC2
  SUBROUTINE DQJAC2(N,X,Fvec,Fjac,Ldfjac,Iflag)
    !> **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Description:**
    !
    !     SUBROUTINE TO EVALUATE THE FULL JACOBIAN FOR TEST PROBLEM USED
    !     IN QUICK CHECK OF DNSQE.
    !
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   890831  Modified array declarations.  (WRB)
    !   890831  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    INTEGER :: Iflag, Ldfjac, N
    REAL(DP) :: Fjac(Ldfjac,N), Fvec(N), X(N)
    !* FIRST EXECUTABLE STATEMENT  DQJAC2
    Fjac(1,1) = -1.0D0
    Fjac(1,2) = 0.0D0
    Fjac(2,1) = -2.0D1*X(1)
    Fjac(2,2) = 1.0D1
  END SUBROUTINE DQJAC2
END MODULE TEST36_MOD
!** TEST36
PROGRAM TEST36
  USE TEST36_MOD, ONLY : DNSQQK, DSOSQX
  USE slatec, ONLY : I1MACH, XSETF, XSETUN, XERMAX
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !> Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  F2
  !***
  ! **Type:**      DOUBLE PRECISION (TEST35-S, TEST36-D)
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
  !        DNSQE    DNSQ     DSOS
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  DNSQQK, DSOSQX, I1MACH, XERMAX, XSETF, XSETUN

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  INTEGER :: ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST36
  lun = I1MACH(2)
  lin = I1MACH(1)
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  CALL XERMAX(1000)
  CALL XSETUN(lun)
  IF( kprint<=1 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  END IF
  !
  !     Test DNSQE and DNSQ
  !
  CALL DNSQQK(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Test DSOS
  !
  CALL DSOSQX(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST36 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST36 *************')
  END IF
  STOP
END PROGRAM TEST36
