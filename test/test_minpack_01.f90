MODULE TEST35_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  !** SNSQQK
  SUBROUTINE SNSQQK(Lun,Kprint,Ipass)
    !> Quick check for SNSQE and SNSQ.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      SINGLE PRECISION (SNSQQK-S, DNSQQK-D)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Description:**
    !
    !   This subroutine performs a quick check on the subroutine SNSQE
    !   (and SNSQ).
    !
    !***
    ! **Routines called:**  ENORM, PASS, R1MACH, SNSQE, SQFCN2, SQJAC2

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891009  Removed unreferenced variable.  (WRB)
    !   891009  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   920310  Code cleaned up and TYPE section added.  (RWC, WRB)
    USE slatec, ONLY : ENORM, R1MACH, SNSQE
    USE common_mod, ONLY : PASS
    !     .. Scalar Arguments ..
    INTEGER :: Ipass, Kprint, Lun
    !     .. Local Scalars ..
    REAL(SP) :: fnorm, fnorms, tol
    INTEGER :: icnt, info, infos, iopt, lwa, n, nprint
    !     .. Local Arrays ..
    REAL(SP) :: fvec(2), wa(19), x(2)
    INTEGER :: itest(3)
    !     .. Intrinsic Functions ..
    INTRINSIC SQRT
    !* FIRST EXECUTABLE STATEMENT  SNSQQK
    infos = 1
    fnorms = 0._SP
    n = 2
    lwa = 19
    nprint = -1
    tol = SQRT(R1MACH(4))
    IF( Kprint>=2 ) WRITE (Lun,99001)
    99001 FORMAT ('1'/'  SNSQE QUICK CHECK'/)
    !
    !     Option 1, the user provides the Jacobian.
    !
    iopt = 1
    x(1) = -1.2_SP
    x(2) = 1._SP
    CALL SNSQE(SQFCN2,SQJAC2,iopt,n,x,fvec,tol,nprint,info,wa,lwa)
    icnt = 1
    fnorm = ENORM(n,fvec)
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
    x(1) = -1.2_SP
    x(2) = 1._SP
    CALL SNSQE(SQFCN2,SQJAC2,iopt,n,x,fvec,tol,nprint,info,wa,lwa)
    icnt = 2
    fnorm = ENORM(n,fvec)
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
    x(1) = -1.2_SP
    x(2) = 1._SP
    CALL SNSQE(SQFCN2,SQJAC2,iopt,n,x,fvec,tol,nprint,info,wa,lwa)
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
    99002 FORMAT (/' **********WARNING -- SNSQE/SNSQ FAILED SOME TESTS****',&
      '******')
    IF( Kprint>=2 .AND. Ipass==1 ) WRITE (Lun,99003)
    99003 FORMAT (/' ----------SNSQE/SNSQ PASSED ALL TESTS----------')
    RETURN
    99004 FORMAT (' EXPECTED VALUE OF INFO AND RESIDUAL NORM',I5,&
      E20.5/' RETURNED VALUE OF INFO AND RESIDUAL NORM',I5,E20.5/)
  END SUBROUTINE SNSQQK
  !** SOSFNC
  REAL(SP) FUNCTION SOSFNC(X,K)
    !> Function evaluator for SOS quick check.
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
    !     FOR TEST PROGRAM USED IN QUICK CHECK OF SOS.
    !
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   801001  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    INTEGER :: K
    REAL(SP) :: X(:)
    !* FIRST EXECUTABLE STATEMENT  SOSFNC
    IF( K==1 ) THEN
      SOSFNC = 1._SP - X(1)
    ELSEIF( K==2 ) THEN
      SOSFNC = 1.E1_SP*(X(2)-X(1)**2)
    ELSE
      SOSFNC = 0._SP
    END IF
  END FUNCTION SOSFNC
  !** SOSNQX
  SUBROUTINE SOSNQX(Lun,Kprint,Ipass)
    !> Quick check for SOS.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      SINGLE PRECISION (SOSNQX-S, DSOSQX-D)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Watts, H. A., (SNLA)
    !***
    ! **Description:**
    !
    !   This subroutine performs a quick check on the subroutine SOS.
    !
    !***
    ! **Routines called:**  PASS, R1MACH, SNRM2, SOS, SOSFNC

    !* REVISION HISTORY  (YYMMDD)
    !   801001  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   920310  Code cleaned up and TYPE section added.  (RWC, WRB)
    USE slatec, ONLY : R1MACH, SOS
    USE common_mod, ONLY : PASS
    !     .. Scalar Arguments ..
    INTEGER :: Ipass, Kprint, Lun
    !     .. Local Scalars ..
    REAL(SP) :: aer, fnorm, fnorms, rer, tolf
    INTEGER :: icnt, iflag, iflags, liw, lwa, n
    !     .. Local Arrays ..
    REAL(SP) :: fvec(2), wa(17), x(2)
    INTEGER :: itest(2), iw(6)
    !     .. Intrinsic Functions ..
    INTRINSIC SQRT
    !* FIRST EXECUTABLE STATEMENT  SOSNQX
    iflags = 3
    fnorms = 0._SP
    n = 2
    lwa = 17
    liw = 6
    tolf = SQRT(R1MACH(4))
    rer = SQRT(R1MACH(4))
    aer = 0._SP
    IF( Kprint>=2 ) WRITE (Lun,99001)
    99001 FORMAT ('1'/'  SOS QUICK CHECK'/)
    !
    !     Test the code with proper input values.
    !
    iflag = 0
    x(1) = -1.2_SP
    x(2) = 1._SP
    CALL SOS(SOSFNC,n,x,rer,aer,tolf,iflag,wa,lwa,iw,liw)
    icnt = 1
    fvec(1) = SOSFNC(x,1)
    fvec(2) = SOSFNC(x,2)
    fnorm = NORM2(fvec(1:n))
    itest(icnt) = 0
    IF( iflag<=iflags .AND. fnorm-fnorms<=rer ) itest(icnt) = 1
    !
    IF( Kprint/=0 ) THEN
      IF( Kprint>=3 .OR. (Kprint>=2 .AND. itest(icnt)/=1) ) WRITE (Lun,99002)&
        iflags, fnorms, iflag, fnorm
      99002 FORMAT (' EXPECTED VALUE OF IFLAG AND RESIDUAL NORM',I5,&
        E20.5/' RETURNED VALUE OF IFLAG AND RESIDUAL NORM',I5,E20.5/)
      IF( Kprint>=2 .OR. (Kprint==1 .AND. itest(icnt)/=1) )&
        CALL PASS(Lun,icnt,itest(icnt))
    END IF
    !
    !     Test improper input parameters.
    !
    lwa = 15
    iflag = 0
    x(1) = -1.2_SP
    x(2) = 1._SP
    CALL SOS(SOSFNC,n,x,rer,aer,tolf,iflag,wa,lwa,iw,liw)
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
    99003 FORMAT (/' **********WARNING -- SOS FAILED SOME TESTS**********')
    IF( Kprint>=2 .AND. Ipass==1 ) WRITE (Lun,99004)
    99004 FORMAT (/' ----------SOS PASSED ALL TESTS----------')
    RETURN
  END SUBROUTINE SOSNQX
  !** SQFCN2
  SUBROUTINE SQFCN2(N,X,Fvec,Iflag)
    !> Evaluate function used in SNSQE.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      SINGLE PRECISION (SQFCN2-S, DQFCN2-D)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Description:**
    !
    !   Subroutine which evaluates the function for test program
    !   used in quick check of SNSQE.
    !
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   930214  TYPE and declarations sections added.  (WRB)

    !     .. Scalar Arguments ..
    INTEGER :: Iflag, N
    !     .. Array Arguments ..
    REAL(SP) :: Fvec(N), X(N)
    !* FIRST EXECUTABLE STATEMENT  SQFCN2
    Fvec(1) = 1._SP - X(1)
    Fvec(2) = 10._SP*(X(2)-X(1)**2)
  END SUBROUTINE SQFCN2
  !** SQJAC2
  SUBROUTINE SQJAC2(N,X,Fvec,Fjac,Ldfjac,Iflag)
    !> Evaluate full Jacobian for SNSQE test.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Description:**
    !
    !     SUBROUTINE TO EVALUATE THE FULL JACOBIAN FOR TEST PROBLEM USED
    !     IN QUICK CHECK OF SNSQE.
    !
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   890831  Modified array declarations.  (WRB)
    !   890831  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    INTEGER :: Iflag, Ldfjac, N
    REAL(SP) :: Fjac(Ldfjac,N), Fvec(N), X(N)
    !* FIRST EXECUTABLE STATEMENT  SQJAC2
    Fjac(1,1) = -1._SP
    Fjac(1,2) = 0._SP
    Fjac(2,1) = -2.E1_SP*X(1)
    Fjac(2,2) = 1.E1_SP
  END SUBROUTINE SQJAC2
END MODULE TEST35_MOD
!** TEST35
PROGRAM TEST35
  USE TEST35_MOD, ONLY : SNSQQK, SOSNQX
  USE slatec, ONLY : I1MACH, XSETF, XSETUN, XERMAX
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !> Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  F2
  !***
  ! **Type:**      SINGLE PRECISION (TEST35-S, TEST36-D)
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
  !        SNSQE    SNSQ     SOS
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  I1MACH, SNSQQK, SOSNQX, XERMAX, XSETF, XSETUN

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  INTEGER :: ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST35
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
  !     Test SNSQE and SNSQ
  !
  CALL SNSQQK(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Test SOS
  !
  CALL SOSNQX(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST35 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST35 *************')
  END IF
  STOP
END PROGRAM TEST35
