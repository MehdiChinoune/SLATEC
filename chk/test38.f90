MODULE TEST38_MOD
  IMPLICIT NONE

CONTAINS
  !DECK DPLPQX
  SUBROUTINE DPLPQX(Lun,Kprint,Ipass)
    !***BEGIN PROLOGUE  DPLPQX
    !***PURPOSE  Quick check for DSPLP.
    !***LIBRARY   SLATEC
    !***TYPE      DOUBLE PRECISION (SPLPQX-S, DPLPQX-D)
    !***AUTHOR  (UNKNOWN)
    !***ROUTINES CALLED  DCOPY, DSPLP, DUSRMT, PASS
    !***REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901013  Added additional printout on failure.  (RWC)
    !***END PROLOGUE  DPLPQX
    IMPLICIT NONE
    REAL(8) :: DUSRMT
    INTEGER i, ic, iv, ivv, j, kk, kount, Kprint, Lun, mm
    EXTERNAL DUSRMT
    INTEGER icnt, ind(60), ibasis(60), Ipass, iwork(900), isoln(14)
    REAL(8) :: costs(37)
    REAL(8) :: prgopt(50), dattrv(210), bl(60), bu(60)
    REAL(8) :: primal(60), duals(60)
    REAL(8) :: work(800)
    REAL(8) :: d(14,37)
    REAL(8) :: zero
    INTEGER mrelas, nvars, info, lw, liw
    !***FIRST EXECUTABLE STATEMENT  DPLPQX
    IF ( Kprint>=2 ) WRITE (Lun,99001)
    99001 FORMAT ('1 DSPLP QUICK CHECK')
    icnt = 1
    zero = 0.0D0
    Ipass = 0
    !     DEFINE WORKING ARRAY LENGTHS
    liw = 900
    lw = 800
    mrelas = 14
    nvars = 37
    !     DEFINE THE ARRAY COSTS(*) FOR THE OBJECTIVE FUNCTION
    costs(1) = 1.030D0
    costs(2) = 0.985D0
    costs(3) = 0.997D0
    costs(4) = 1.036D0
    costs(5) = 1.005D0
    costs(6) = 0.980D0
    costs(7) = 1.004D0
    costs(8) = 0.993D0
    costs(9) = 1.018D0
    costs(10) = 0.947D0
    costs(11) = 0.910D0
    costs(12) = 1.028D0
    costs(13) = 0.957D0
    costs(14) = 1.025D0
    costs(15) = 1.036D0
    costs(16) = 1.060D0
    costs(17) = 0.954D0
    costs(18) = 0.891D0
    costs(19) = 0.921D0
    costs(20) = 1.040D0
    costs(21) = 0.912D0
    costs(22) = 0.926D0
    costs(23) = 1.000D0
    costs(24) = 0.000D0
    costs(25) = 0.000D0
    costs(26) = 0.000D0
    costs(27) = 0.000D0
    costs(28) = 0.000D0
    costs(29) = 0.000D0
    costs(30) = 0.000D0
    costs(31) = 0.000D0
    costs(32) = 0.000D0
    costs(33) = 0.000D0
    costs(34) = 0.000D0
    costs(35) = 0.000D0
    costs(36) = 0.000D0
    costs(37) = 0.000D0
    !     PLACE THE NONZERO INFORMATION ABOUT THE MATRIX IN DATTRV(*)
    CALL DCOPY(14*37,zero,0,d,1)
    d(1,1) = 1.04000D0
    d(1,23) = 1.00000D0
    d(1,24) = -1.00000D0
    d(2,6) = 0.04125D0
    d(2,7) = 0.05250D0
    d(2,17) = 0.04875D0
    d(2,24) = 1.00000D0
    d(2,25) = -1.00000D0
    d(3,8) = 0.05625D0
    d(3,9) = 0.06875D0
    d(3,11) = 0.02250D0
    d(3,25) = 1.00000D0
    d(3,26) = -1.00000D0
    d(4,2) = 1.04000D0
    d(4,3) = 1.05375D0
    d(4,5) = 1.06125D0
    d(4,12) = 0.08000D0
    d(4,16) = 0.09375D0
    d(4,18) = 0.03750D0
    d(4,19) = 0.04625D0
    d(4,20) = 0.08125D0
    d(4,22) = 0.05250D0
    d(4,26) = 1.00000D0
    d(4,27) = -1.00000D0
    d(5,10) = 0.04375D0
    d(5,27) = 1.00000D0
    d(5,28) = -1.00000D0
    d(6,4) = 1.05875D0
    d(6,13) = 0.04500D0
    d(6,14) = 0.06375D0
    d(6,15) = 0.06625D0
    d(6,21) = 0.05000D0
    d(6,28) = 1.00000D0
    d(6,29) = -1.00000D0
    d(7,6) = 1.04125D0
    d(7,7) = 1.05250D0
    d(7,8) = 1.05625D0
    d(7,9) = 1.06875D0
    d(7,11) = 0.02250D0
    d(7,17) = 0.04875D0
    d(7,29) = 1.00000D0
    d(7,30) = -1.00000D0
    d(8,10) = 1.04375D0
    d(8,12) = 0.08000D0
    d(8,13) = 0.04500D0
    d(8,14) = 0.06375D0
    d(8,15) = 0.06625D0
    d(8,16) = 0.09375D0
    d(8,18) = 0.03750D0
    d(8,19) = 0.04625D0
    d(8,20) = 0.08125D0
    d(8,21) = 0.05000D0
    d(8,22) = 0.05250D0
    d(8,30) = 1.00000D0
    d(8,31) = -1.00000D0
    d(9,11) = 1.02250D0
    d(9,17) = 0.04875D0
    d(9,31) = 1.00000D0
    d(9,32) = -1.00000D0
    d(10,12) = 1.08000D0
    d(10,13) = 1.04500D0
    d(10,14) = 1.06375D0
    d(10,15) = 1.06625D0
    d(10,16) = 1.09375D0
    d(10,18) = 0.03750D0
    d(10,19) = 0.04625D0
    d(10,20) = 0.08125D0
    d(10,21) = 0.05000D0
    d(10,22) = 0.05250D0
    d(10,32) = 1.00000D0
    d(10,33) = -1.00000D0
    d(11,17) = 1.04875D0
    d(11,33) = 1.00000D0
    d(11,34) = -1.00000D0
    d(12,18) = 1.03750D0
    d(12,19) = 1.04625D0
    d(12,20) = 1.08125D0
    d(12,21) = 1.05000D0
    d(12,22) = 0.05250D0
    d(12,34) = 1.00000D0
    d(12,35) = -1.00000D0
    d(13,35) = 1.00000D0
    d(13,36) = -1.00000D0
    d(14,22) = 1.05250D0
    d(14,36) = 1.00000D0
    d(14,37) = -1.00000D0
    kount = 1
    DO mm = 1, nvars
      dattrv(kount) = -mm
      DO kk = 1, mrelas
        IF ( d(kk,mm)/=zero ) THEN
          kount = kount + 1
          dattrv(kount) = kk
          kount = kount + 1
          dattrv(kount) = d(kk,mm)
        ENDIF
      ENDDO
      kount = kount + 1
    ENDDO
    dattrv(kount) = zero
    !     NON-NEGATIVITY CONSTRAINT
    DO ic = 1, nvars
      bl(ic) = zero
      ind(ic) = 3
      bu(ic) = 10000000.000D0
    ENDDO
    !     LE CONSTRAINTS
    DO iv = 1, mrelas
      ivv = iv + nvars
      ind(ivv) = 3
      bl(ivv) = 100.00000D0
      bu(ivv) = 100000000.00000D0
    ENDDO
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
    CALL DSPLP(DUSRMT,mrelas,nvars,costs,prgopt,dattrv,bl,bu,ind,info,primal,&
      duals,ibasis,work,lw,iwork,liw)
    !
    !     LOOK FOR THE KNOWN BASIS AT THE SOLN., NOW IS ISOLN(*).
    !
    DO i = 1, mrelas
      isoln(i) = prgopt(i+3)
    ENDDO
    !
    Ipass = 1
    DO j = 1, mrelas
      DO i = 1, mrelas
        IF ( isoln(i)==ibasis(j) ) GOTO 100
      ENDDO
      Ipass = 0
      EXIT
      100  ENDDO
      !
      IF ( Kprint>=2 ) WRITE (Lun,99002) (isoln(i),ibasis(i),i=1,mrelas)
      !
      99002 FORMAT (/'     ISOLN    IBASIS'/(2I10))
      !
      IF ( Kprint>=2.OR.(Kprint==1.AND.Ipass/=1) ) CALL PASS(Lun,icnt,Ipass)
      !
      !     HERE IPASS=0 IF CODE FAILED QUICK CHECK;
      !               =1 IF CODE PASSED QUICK CHECK.
      !
      IF ( Kprint>=1.AND.Ipass/=1 ) WRITE (Lun,99003)
      99003 FORMAT (/' ************ DSPLP FAILED SOME TESTS ***************')
      IF ( Kprint>=2.AND.Ipass==1 ) WRITE (Lun,99004)
      99004 FORMAT (/' ************ DSPLP PASSED ALL TESTS ****************')
      RETURN
  END SUBROUTINE DPLPQX
  !DECK DBOCQX
  SUBROUTINE DBOCQX(Lun,Kprint,Ipass)
    !***BEGIN PROLOGUE  DBOCQX
    !***PURPOSE  Quick check for DBOCLS.
    !***LIBRARY   SLATEC
    !***TYPE      DOUBLE PRECISION (SBOCQX-S, DBOCQX-D)
    !***KEYWORDS  QUICK CHECK
    !***AUTHOR  (UNKNOWN)
    !***DESCRIPTION
    !
    !     MINIMAL TEST DRIVER FOR DBOCLS, BOUNDED CONSTRAINED LEAST
    !     SQUARES SOLVER.  DELIVERS THE VALUE IPASS=1 IF 8 TESTS WERE
    !     PASSED.  DELIVER THE VALUE IPASS=0 IF ANY ONE OF THEM FAILED.
    !
    !     RUN FOUR BOUNDED LEAST SQUARES PROBLEMS THAT COME FROM THE
    !     DIPLOME WORK OF P. ZIMMERMANN.
    !
    !***ROUTINES CALLED  D1MACH, DBOCLS, DBOLS, DCOPY, DNRM2
    !***REVISION HISTORY  (YYMMDD)
    !   850310  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Added PASS/FAIL message.  (RWC)
    !***END PROLOGUE  DBOCQX
    IMPLICIT NONE
    REAL(8) :: D1MACH, DNRM2, rnorm, rnormc, sr
    INTEGER i, ib, Ipass, irhs, itest, j, Kprint, Lun, mcon, mdw ,&
      mode, mpass, mrows, ncols
    REAL(8) :: d(6,5), w(11,11), bl(5,2), bu(5,2), x(30), rw(55) ,&
      xtrue(9)
    REAL(8) :: c(5,5)
    REAL(8) :: bl1(10), bu1(10)
    INTEGER ind(10), iw(20), iopt(40)
    REAL(8) :: rhs(6,2)
    CHARACTER(4) :: msg
    !
    DATA ((c(i,j),i=1,5),j=1,5)/1.D0, 10.D0, 4.D0, 8.D0, 1.D0, 1.D0 ,&
      10.D0, 2.D0, -1.D0, 1.D0, 1.D0, -3.D0, -3.D0, 2.D0, 1.D0 ,&
      1.D0, 5.D0, 5.D0, 5.D0, 1.D0, 1.D0, 4.D0, -1.D0, -3.D0 ,&
      1.D0/
    DATA ((d(i,j),i=1,6),j=1,5)/ - 74.D0, 14.D0, 66.D0, -12.D0, 3.D0 ,&
      4.D0, 80.D0, -69.D0, -72.D0, 66.D0, 8.D0, -12.D0, 18.D0 ,&
      21.D0, -5.D0, -30.D0, -7.D0, 4.D0, -11.D0, 28.D0, 7.D0 ,&
      -23.D0, -4.D0, 4.D0, -4.D0, 0.D0, 1.D0, 3.D0, 1.D0, 0.D0/
    DATA ((bl(i,j),i=1,5),j=1,2)/1.D0, 0.D0, -1.D0, 1.D0, -4.D0, -1.D0 ,&
      0.D0, -3.D0, 1.D0, -6.D0/
    DATA ((bu(i,j),i=1,5),j=1,2)/3.D0, 2.D0, 1.D0, 3.D0, -2.D0, 3.D0 ,&
      4.D0, 1.D0, 5.D0, -2.D0/
    DATA ((rhs(i,j),i=1,6),j=1,2)/51.D0, -61.D0, -56.D0, 69.D0, 10.D0 ,&
      -12.D0, -5.D0, -9.D0, 708.D0, 4165.D0, -13266.D0, 8409.D0/
    DATA (xtrue(j),j=1,9)/1.D0, 2.D0, -1.D0, 3.D0, -4.D0, 1.D0, 32.D0 ,&
      30.D0, 31.D0/
    !***FIRST EXECUTABLE STATEMENT  DBOCQX
    mdw = 11
    mrows = 6
    ncols = 5
    mcon = 4
    iopt(1) = 99
    Ipass = 1
    itest = 0
    !
    IF ( Kprint>=2 ) WRITE (Lun,99001)
    99001 FORMAT (' TEST   IB IRHS             SR')
    !
    DO ib = 1, 2
      DO irhs = 1, 2
        !
        !           TRANSFER DATA TO WORKING ARRAY W(*,*).
        !
        DO j = 1, ncols
          CALL DCOPY(mrows,d(1,j),1,w(1,j),1)
        ENDDO
        !
        CALL DCOPY(mrows,rhs(1,irhs),1,w(1,ncols+1),1)
        !
        !             SET BOUND INDICATOR FLAGS.
        !
        DO j = 1, ncols
          ind(j) = 3
        ENDDO
        !
        CALL DBOLS(w,mdw,mrows,ncols,bl(1,ib),bu(1,ib),ind,iopt,x,rnorm,mode,&
          rw,iw)
        DO j = 1, ncols
          x(j) = x(j) - xtrue(j)
        ENDDO
        !
        sr = DNRM2(ncols,x,1)
        mpass = 1
        IF ( sr>10.D2*SQRT(D1MACH(4)) ) mpass = 0
        Ipass = Ipass*mpass
        IF ( Kprint>=2 ) THEN
          msg = 'PASS'
          IF ( mpass==0 ) msg = 'FAIL'
          itest = itest + 1
          WRITE (Lun,99003) itest, ib, irhs, sr, msg
        ENDIF
      ENDDO
    ENDDO
    !
    !     RUN STOER'S PROBLEM FROM 1971 SIAM J. N. ANAL. PAPER.
    !
    DO ib = 1, 2
      DO irhs = 1, 2
        CALL DCOPY(11*10,0.D0,0,w,1)
        CALL DCOPY(ncols,bl(1,ib),1,bl1,1)
        CALL DCOPY(ncols,bu(1,ib),1,bu1,1)
        ind(ncols+1) = 2
        ind(ncols+2) = 1
        ind(ncols+3) = 2
        ind(ncols+4) = 3
        bu1(ncols+1) = 5.
        bl1(ncols+2) = 20.
        bu1(ncols+3) = 30.
        bl1(ncols+4) = 11.
        bu1(ncols+4) = 40.
        DO j = 1, ncols
          CALL DCOPY(mcon,c(1,j),1,w(1,j),1)
          CALL DCOPY(mrows,d(1,j),1,w(mcon+1,j),1)
        ENDDO
        !
        CALL DCOPY(mrows,rhs(1,irhs),1,w(mcon+1,ncols+1),1)
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
        CALL DBOCLS(w,mdw,mcon,mrows,ncols,bl1,bu1,ind,iopt,x,rnormc,rnorm,&
          mode,rw,iw)
        DO j = 1, ncols + mcon
          x(j) = x(j) - xtrue(j)
        ENDDO
        !
        sr = DNRM2(ncols+mcon,x,1)
        mpass = 1
        IF ( sr>10.D2*SQRT(D1MACH(4)) ) mpass = 0
        Ipass = Ipass*mpass
        IF ( Kprint>=2 ) THEN
          msg = 'PASS'
          IF ( mpass==0 ) msg = 'FAIL'
          itest = itest + 1
          WRITE (Lun,99003) itest, ib, irhs, sr, msg
        ENDIF
      ENDDO
    ENDDO
    !
    !     HERE THE VALUE OF IPASS=1 SAYS THAT DBOCLS HAS PASSED ITS TESTS.
    !          THE VALUE OF IPASS=0 SAYS THAT DBOCLS HAS NOT PASSED.
    !
    IF ( Kprint>=3 ) WRITE (Lun,&
      '('' IPASS VALUE. (A 1 IS GOOD, 0 IS BAD.)'',I4)')&
      Ipass
    IF ( Kprint>=2.AND.Ipass==0 ) WRITE (Lun,99002)
    !
    99002 FORMAT (' ERROR IN DBOCLS OR DBOLS')
    RETURN
    99003 FORMAT (3I5,1P,E20.6,' TEST ',A,'ED.')
  END SUBROUTINE DBOCQX
END MODULE TEST38_MOD
!DECK TEST38
PROGRAM TEST38
  USE TEST38_MOD
  IMPLICIT NONE
  INTEGER I1MACH
  !***BEGIN PROLOGUE  TEST38
  !***PURPOSE  Driver for testing SLATEC subprograms
  !***LIBRARY   SLATEC
  !***CATEGORY  G2
  !***TYPE      DOUBLE PRECISION (TEST37-S, TEST38-D)
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
  !        DSPLP    DBOCLS
  !
  !***REFERENCES  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***ROUTINES CALLED  DBOCQX, DPLPQX, I1MACH, XERMAX, XSETF, XSETUN
  !***REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  !***END PROLOGUE  TEST38
  INTEGER ipass, kprint, lin, lun, nfail
  !***FIRST EXECUTABLE STATEMENT  TEST38
  lun = I1MACH(2)
  lin = I1MACH(1)
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  READ (lin,'(I1)') kprint
  CALL XERMAX(1000)
  CALL XSETUN(lun)
  IF ( kprint<=1 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  ENDIF
  !
  !     Test DSPLP package
  !
  CALL DPLPQX(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test DBOCLS package
  !
  CALL DBOCQX(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001   FORMAT (/' --------------TEST38 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002   FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST38 *************')
  ENDIF
  STOP
END PROGRAM TEST38
