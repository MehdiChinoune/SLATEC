MODULE TEST37_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  !** SPLPQX
  SUBROUTINE SPLPQX(Lun,Kprint,Ipass)
    !> Quick check for SPLP.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      SINGLE PRECISION (SPLPQX-S, DPLPQX-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  PASS, SCOPY, SPLP, USRMAT

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901013  Added additional printout on failure.  (RWC)
    USE optimization, ONLY : USRMAT, SPLP
    USE common_mod, ONLY : PASS
    !
    REAL(SP) :: bl(60), bu(60), d(14,37), dattrv(210), duals(60), prgopt(50), &
      primal(60), work(800)
    INTEGER :: i, ibasis(60), ic, icnt, ind(60), info, Ipass, isoln(14), iv, ivv, &
      iwork(900), j, kk, kount, Kprint, liw, Lun, lw, mm, mrelas
    INTEGER :: nvars
    REAL(SP) :: costs(37)
    !* FIRST EXECUTABLE STATEMENT  SPLPQX
    IF( Kprint>=2 ) WRITE (Lun,99001)
    99001 FORMAT ('1 SPLP QUICK CHECK')
    icnt = 1
    !
    !     DEFINE WORKING ARRAY LENGTHS
    !
    liw = 900
    lw = 800
    mrelas = 14
    nvars = 37
    !
    !     DEFINE THE ARRAY COSTS(*) FOR THE OBJECTIVE FUNCTION
    !
    costs(1) = 1.030_SP
    costs(2) = 0.985_SP
    costs(3) = 0.997_SP
    costs(4) = 1.036_SP
    costs(5) = 1.005_SP
    costs(6) = 0.980_SP
    costs(7) = 1.004_SP
    costs(8) = 0.993_SP
    costs(9) = 1.018_SP
    costs(10) = 0.947_SP
    costs(11) = 0.910_SP
    costs(12) = 1.028_SP
    costs(13) = 0.957_SP
    costs(14) = 1.025_SP
    costs(15) = 1.036_SP
    costs(16) = 1.060_SP
    costs(17) = 0.954_SP
    costs(18) = 0.891_SP
    costs(19) = 0.921_SP
    costs(20) = 1.040_SP
    costs(21) = 0.912_SP
    costs(22) = 0.926_SP
    costs(23) = 1.000_SP
    costs(24) = 0.000_SP
    costs(25) = 0.000_SP
    costs(26) = 0.000_SP
    costs(27) = 0.000_SP
    costs(28) = 0.000_SP
    costs(29) = 0.000_SP
    costs(30) = 0.000_SP
    costs(31) = 0.000_SP
    costs(32) = 0.000_SP
    costs(33) = 0.000_SP
    costs(34) = 0.000_SP
    costs(35) = 0.000_SP
    costs(36) = 0.000_SP
    costs(37) = 0.000_SP
    !
    !     PLACE THE NONZERO INFORMATION ABOUT THE MATRIX IN DATTRV(*)
    !
    d = 0._SP
    d(1,1) = 1.04000_SP
    d(1,23) = 1.00000_SP
    d(1,24) = -1.00000_SP
    d(2,6) = 0.04125_SP
    d(2,7) = 0.05250_SP
    d(2,17) = 0.04875_SP
    d(2,24) = 1.00000_SP
    d(2,25) = -1.00000_SP
    d(3,8) = 0.05625_SP
    d(3,9) = 0.06875_SP
    d(3,11) = 0.02250_SP
    d(3,25) = 1.00000_SP
    d(3,26) = -1.00000_SP
    d(4,2) = 1.04000_SP
    d(4,3) = 1.05375_SP
    d(4,5) = 1.06125_SP
    d(4,12) = 0.08000_SP
    d(4,16) = 0.09375_SP
    d(4,18) = 0.03750_SP
    d(4,19) = 0.04625_SP
    d(4,20) = 0.08125_SP
    d(4,22) = 0.05250_SP
    d(4,26) = 1.00000_SP
    d(4,27) = -1.00000_SP
    d(5,10) = 0.04375_SP
    d(5,27) = 1.00000_SP
    d(5,28) = -1.00000_SP
    d(6,4) = 1.05875_SP
    d(6,13) = 0.04500_SP
    d(6,14) = 0.06375_SP
    d(6,15) = 0.06625_SP
    d(6,21) = 0.05000_SP
    d(6,28) = 1.00000_SP
    d(6,29) = -1.00000_SP
    d(7,6) = 1.04125_SP
    d(7,7) = 1.05250_SP
    d(7,8) = 1.05625_SP
    d(7,9) = 1.06875_SP
    d(7,11) = 0.02250_SP
    d(7,17) = 0.04875_SP
    d(7,29) = 1.00000_SP
    d(7,30) = -1.00000_SP
    d(8,10) = 1.04375_SP
    d(8,12) = 0.08000_SP
    d(8,13) = 0.04500_SP
    d(8,14) = 0.06375_SP
    d(8,15) = 0.06625_SP
    d(8,16) = 0.09375_SP
    d(8,18) = 0.03750_SP
    d(8,19) = 0.04625_SP
    d(8,20) = 0.08125_SP
    d(8,21) = 0.05000_SP
    d(8,22) = 0.05250_SP
    d(8,30) = 1.00000_SP
    d(8,31) = -1.00000_SP
    d(9,11) = 1.02250_SP
    d(9,17) = 0.04875_SP
    d(9,31) = 1.00000_SP
    d(9,32) = -1.00000_SP
    d(10,12) = 1.08000_SP
    d(10,13) = 1.04500_SP
    d(10,14) = 1.06375_SP
    d(10,15) = 1.06625_SP
    d(10,16) = 1.09375_SP
    d(10,18) = 0.03750_SP
    d(10,19) = 0.04625_SP
    d(10,20) = 0.08125_SP
    d(10,21) = 0.05000_SP
    d(10,22) = 0.05250_SP
    d(10,32) = 1.00000_SP
    d(10,33) = -1.00000_SP
    d(11,17) = 1.04875_SP
    d(11,33) = 1.00000_SP
    d(11,34) = -1.00000_SP
    d(12,18) = 1.03750_SP
    d(12,19) = 1.04625_SP
    d(12,20) = 1.08125_SP
    d(12,21) = 1.05000_SP
    d(12,22) = 0.05250_SP
    d(12,34) = 1.00000_SP
    d(12,35) = -1.00000_SP
    d(13,35) = 1.00000_SP
    d(13,36) = -1.00000_SP
    d(14,22) = 1.05250_SP
    d(14,36) = 1.00000_SP
    d(14,37) = -1.00000_SP
    kount = 1
    DO mm = 1, nvars
      dattrv(kount) = -mm
      DO kk = 1, mrelas
        IF( d(kk,mm)/=0._SP ) THEN
          kount = kount + 1
          dattrv(kount) = kk
          kount = kount + 1
          dattrv(kount) = d(kk,mm)
        END IF
      END DO
      kount = kount + 1
    END DO
    dattrv(kount) = 0._SP
    !
    !     NON-NEGATIVITY CONSTRAINT
    !
    DO ic = 1, nvars
      bl(ic) = 0._SP
      ind(ic) = 3
      bu(ic) = 10000000.000
    END DO
    !
    !     LE CONSTRAINTS
    !
    DO iv = 1, mrelas
      ivv = iv + nvars
      ind(ivv) = 3
      bl(ivv) = 100.00000
      bu(ivv) = 100000000.00000
    END DO
    prgopt(01) = 18
    prgopt(02) = 59
    prgopt(03) = 0
    prgopt(04) = 1
    prgopt(05) = 3
    prgopt(06) = 8
    prgopt(07) = 10
    prgopt(08) = 11
    prgopt(09) = 16
    prgopt(10) = 17
    prgopt(11) = 21
    prgopt(12) = 22
    prgopt(13) = 24
    prgopt(14) = 25
    prgopt(15) = 27
    prgopt(16) = 28
    prgopt(17) = 35
    prgopt(18) = 21
    prgopt(19) = 51
    prgopt(20) = 0
    prgopt(21) = 1
    CALL SPLP(USRMAT,mrelas,nvars,costs,prgopt,dattrv,bl,bu,ind,info,primal,&
      duals,ibasis,work,lw,iwork,liw)
    !
    !     LOOK FOR THE KNOWN BASIS AT THE SOLN., NOW IS ISOLN(*).
    !
    DO i = 1, mrelas
      isoln(i) = INT( prgopt(i+3) )
    END DO
    !
    Ipass = 1
    DO j = 1, mrelas
      DO i = 1, mrelas
        IF( isoln(i)==ibasis(j) ) GOTO 100
      END DO
      Ipass = 0
      EXIT
      100 CONTINUE
    END DO
    !
    IF( Kprint>=2 ) WRITE (Lun,99002) (isoln(i),ibasis(i),i=1,mrelas)
    !
    99002 FORMAT (/'     ISOLN    IBASIS'/(2I10))
    !
    IF( Kprint>=2 .OR. (Kprint==1 .AND. Ipass/=1) ) CALL PASS(Lun,icnt,Ipass)
    !
    !     HERE IPASS=0 IF CODE FAILED QUICK CHECK;
    !               =1 IF CODE PASSED QUICK CHECK.
    !
    IF( Kprint>=1 .AND. Ipass/=1 ) WRITE (Lun,99003)
    99003 FORMAT (/' ************ SPLP FAILED SOME TESTS ****************')
    IF( Kprint>=2 .AND. Ipass==1 ) WRITE (Lun,99004)
    99004 FORMAT (/' ************ SPLP PASSED ALL TESTS *****************')
    RETURN
  END SUBROUTINE SPLPQX
  !** SBOCQX
  SUBROUTINE SBOCQX(Lun,Kprint,Ipass)
    !> Quick check for SBOCLS.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      SINGLE PRECISION (SBOCQX-S, DBOCQX-D)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Description:**
    !
    !     MINIMAL TEST DRIVER FOR SBOCLS, BOUNDED CONSTRAINED LEAST
    !     SQUARES SOLVER.  DELIVERS THE VALUE IPASS=1 IF 8 TESTS WERE
    !     PASSED.  DELIVER THE VALUE IPASS=0 IF ANY ONE OF THEM FAILED.
    !
    !     RUN FOUR BOUNDED LEAST SQUARES PROBLEMS THAT COME FROM THE
    !     DIPLOME WORK OF P. ZIMMERMANN.
    !
    !***
    ! **Routines called:**  R1MACH, SBOCLS, SBOLS, SCOPY, SNRM2

    !* REVISION HISTORY  (YYMMDD)
    !   850310  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901013  Added PASS/FAIL message and cleaned up FORMATs.  (RWC)
    USE service, ONLY : eps_sp
    USE approximation, ONLY : SBOLS, SBOCLS
    !
    INTEGER :: ib, Ipass, irhs, itest, j, Kprint, Lun, mcon, mdw, &
      mode, mpass, mrows, ncols
    REAL(SP) :: rnorm, rnormc, sr
    REAL(SP) :: w(11,11), x(30), rw(55), bl1(10), bu1(10)
    INTEGER :: ind(10), iw(20), iopt(40)
    CHARACTER(4) :: msg
    !
    REAL(SP), PARAMETER :: c(5,5) = RESHAPE( [ 1._SP, 10._SP, 4._SP, 8._SP, 1._SP, &
      1._SP, 10._SP, 2._SP, -1._SP, 1._SP,    1._SP, -3._SP, -3._SP, 2._SP, 1._SP, &
      1._SP, 5._SP, 5._SP, 5._SP, 1._SP,      1._SP,  4._SP, -1._SP, -3._SP, 1._SP ], [5,5] )
    REAL(SP), PARAMETER :: d(6,5) = RESHAPE( [ &
      -74._SP, 14._SP, 66._SP, -12._SP, 3._SP, 4._SP, &
      80._SP, -69._SP, -72._SP, 66._SP, 8._SP, -12._SP, &
      18._SP, 21._SP, -5._SP, -30._SP, -7._SP, 4._SP, &
      -11._SP, 28._SP, 7._SP, -23._SP, -4._SP, 4._SP, &
      -4._SP, 0._SP, 1._SP, 3._SP, 1._SP, 0._SP ], [6,5] )
    REAL(SP), PARAMETER :: bl(5,2) = RESHAPE( [ 1._SP, 0._SP, -1._SP, 1._SP, -4._SP, &
      -1._SP, 0._SP, -3._SP, 1._SP, -6._SP ], [5,2] )
    REAL(SP), PARAMETER :: bu(5,2) = RESHAPE( [ 3._SP, 2._SP, 1._SP, 3._SP, -2._SP, &
      3._SP, 4._SP, 1._SP, 5._SP, -2._SP ], [5,2] )
    REAL(SP), PARAMETER :: rhs(6,2) = RESHAPE( [ &
      51._SP, -61._SP, -56._SP, 69._SP, 10._SP, -12._SP, &
      -5._SP, -9._SP, 708._SP, 4165._SP, -13266._SP, 8409._SP ], [6,2] )
    REAL(SP), PARAMETER :: xtrue(9) = [ 1._SP, 2._SP, -1._SP, 3._SP, -4._SP, &
      1._SP, 32._SP, 30._SP, 31._SP ]
    !* FIRST EXECUTABLE STATEMENT  SBOCQX
    mdw = 11
    mrows = 6
    ncols = 5
    mcon = 4
    iopt(1) = 99
    Ipass = 1
    itest = 0
    ! Initialize W
    w = 0._SP
    !
    IF( Kprint>=2 ) WRITE (Lun,99001)
    99001 FORMAT (' TEST   IB IRHS             SR')
    !
    DO ib = 1, 2
      DO irhs = 1, 2
        !
        !           TRANSFER DATA TO WORKING ARRAY W(*,*).
        !
        w(1:mrows,1:ncols) = d(1:mrows,1:ncols)
        !
        w(1:mrows,ncols+1) = rhs(1:mrows,irhs)
        !
        !             SET BOUND INDICATOR FLAGS.
        !
        DO j = 1, ncols
          ind(j) = 3
        END DO
        !
        CALL SBOLS(w,mdw,mrows,ncols,bl(1,ib),bu(1,ib),ind,iopt,x,rnorm,mode,&
          rw,iw)
        DO j = 1, ncols
          x(j) = x(j) - xtrue(j)
        END DO
        !
        sr = NORM2(x(1:ncols))
        mpass = 1
        IF( sr>10.E3*SQRT(eps_sp) ) mpass = 0
        Ipass = Ipass*mpass
        IF( Kprint>=2 ) THEN
          msg = 'PASS'
          IF( mpass==0 ) msg = 'FAIL'
          itest = itest + 1
          WRITE (Lun,99003) itest, ib, irhs, sr, msg
        END IF
      END DO
    END DO
    !
    !     RUN STOER'S PROBLEM FROM 1971 SIAM J. N. ANAL. PAPER.
    !
    DO ib = 1, 2
      DO irhs = 1, 2
        w(:,1:10) = 0._SP
        bl1(1:ncols) = bl(1:ncols,ib)
        bu1(1:ncols) = bu(1:ncols,ib)
        ind(ncols+1) = 2
        ind(ncols+2) = 1
        ind(ncols+3) = 2
        ind(ncols+4) = 3
        bu1(ncols+1) = 5.
        bl1(ncols+2) = 20._SP
        bu1(ncols+3) = 30._SP
        bl1(ncols+4) = 11._SP
        bu1(ncols+4) = 40._SP
        w(1:mcon,1:ncols) = c(1:mcon,1:ncols)
        w(mcon+1:mcon+mrows,1:ncols) = d(1:mrows,1:ncols)
        !
        w(mcon+1:mcon+mrows,ncols+1) = rhs(1:mrows,irhs)
        !
        !           CHECK LENGTHS OF REQD. ARRAYS.
        !
        iopt(01) = 2
        iopt(02) = 11
        iopt(03) = 11
        iopt(04) = 10
        iopt(05) = 30
        iopt(06) = 55
        iopt(07) = 20
        iopt(08) = 40
        iopt(09) = 99
        CALL SBOCLS(w,mdw,mcon,mrows,ncols,bl1,bu1,ind,iopt,x,rnormc,rnorm,&
          mode,rw,iw)
        DO j = 1, ncols + mcon
          x(j) = x(j) - xtrue(j)
        END DO
        !
        sr = NORM2(x(1:ncols+mcon))
        mpass = 1
        IF( sr>10.E3*SQRT(eps_sp) ) mpass = 0
        Ipass = Ipass*mpass
        IF( Kprint>=2 ) THEN
          msg = 'PASS'
          IF( mpass==0 ) msg = 'FAIL'
          itest = itest + 1
          WRITE (Lun,99003) itest, ib, irhs, sr, msg
        END IF
      END DO
    END DO
    !
    !     HERE THE VALUE OF IPASS=1 SAYS THAT SBOCLS() HAS PASSED ITS TESTS.
    !          THE VALUE OF IPASS=0 SAYS THAT SBOCLS() HAS NOT PASSED.
    !
    IF( Kprint>=3 ) WRITE (Lun,&
      '('' IPASS VALUE. (A 1 IS GOOD, 0 IS BAD.)'',I4)') Ipass
    IF( Kprint>=2 .AND. Ipass==0 ) WRITE (Lun,99002)
    !
    99002 FORMAT (' ERROR IN SBOCLS OR SBOLS')
    RETURN
    99003 FORMAT (3I5,1P,E20.6,' TEST ',A,'ED.')
  END SUBROUTINE SBOCQX
END MODULE TEST37_MOD
!** TEST37
PROGRAM TEST37
  USE TEST37_MOD, ONLY : SBOCQX, SPLPQX
  USE ISO_FORTRAN_ENV, ONLY : INPUT_UNIT, OUTPUT_UNIT
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !> Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  G2
  !***
  ! **Type:**      SINGLE PRECISION (TEST37-S, TEST38-D)
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
  !        SPLP     SBOCLS
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  I1MACH, SBOCQX, SPLPQX, XERMAX, XSETF

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  INTEGER :: ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST37
  lun = OUTPUT_UNIT
  lin = INPUT_UNIT
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  !
  !     Test SPLP package
  !
  CALL SPLPQX(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Test SBOCLS package
  !
  CALL SBOCQX(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST37 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST37 *************')
  END IF
  STOP
END PROGRAM TEST37
