MODULE TEST30_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  !** BSPCK
  SUBROUTINE BSPCK(Lun,Kprint,Ipass)
    !> Quick check for the B-Spline package.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      SINGLE PRECISION (BSPCK-S, DBSPCK-D)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Description:**
    !
    !   BSPCK is a quick check routine for the B-Spline package which
    !   tests consistency between results from higher level routines.
    !   Those routines not explicitly called are exercised at some lower
    !   level.  The routines exercised are BFQAD, BINT4, BINTK, BNFAC,
    !   BNSLV, BSGQ8, BSPDR, BSPEV, BSPPP, BSPVD, BSPVN, BSQAD, BVALU,
    !   INTRV, PFQAD, PPGQ8, PPQAD and PPVAL.
    !
    !***
    ! **Routines called:**  BFQAD, BINT4, BINTK, BSPDR, BSPEV, BSPPP, BSPVD,
    !                    BSPVN, BSQAD, BVALU, FB, INTRV, PFQAD, PPQAD,
    !                    PPVAL, R1MACH

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   891004  Removed unreachable code.  (WRB)
    !   891009  Removed unreferenced variables.  (WRB)
    !   891009  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   930214  Declarations sections added, code revised to test error
    !           returns for all values of KPRINT and code polished.  (WRB)
    USE slatec, ONLY : BFQAD, BINT4, BINTK, BSPDR, BSPEV, BSPPP, BSPVD, BSPVN, &
      BSQAD, BVALU, INTRV, PFQAD, PPQAD, PPVAL, R1MACH
    !     .. Scalar Arguments ..
    INTEGER :: Ipass, Kprint, Lun
    !     .. Local Scalars ..
    REAL(SP) :: atol, bquad, bv, den, dn, er, fbcl, fbcr, pi, pquad, quad, &
      spv, tol, x1, x2, xl, xx
    INTEGER :: i, ibcl, ibcr, id, ierr, iknt, ileft, ilo, inbv, inev, &
      inppv, iwork, j, jhigh, jj, k, kk, knt, kntopt, kontrl, &
      ldc, ldcc, lxi, mflag, n, ndata, nmk, nn
    LOGICAL :: fatal
    !     .. Local Arrays ..
    REAL(SP) :: adif(52), bc(13), c(4,10), cc(4,4), q(3), qq(77), qsave(2), &
      sv(4), t(17), w(65), x(11), xi(11), y(11)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, SIN
    !* FIRST EXECUTABLE STATEMENT  BSPCK
    IF( Kprint>=2 ) WRITE (Lun,99001)
    !
    99001 FORMAT ('1 QUICK CHECK FOR SPLINE ROUTINES',//)
    !
    Ipass = 1
    pi = 3.14159265358979324_SP
    tol = 1000._SP*R1MACH(4)
    !
    !     Generate data.
    !
    ndata = 11
    den = ndata - 1
    DO i = 1, ndata
      x(i) = (i-1)/den
      y(i) = SIN(pi*x(i))
    END DO
    x(3) = 2._SP/den
    y(3) = SIN(pi*x(3))
    !
    !     Compute splines for two knot arrays.
    !
    DO iknt = 1, 2
      knt = 3 - iknt
      ibcl = 1
      ibcr = 2
      fbcl = pi
      fbcr = 0._SP
      CALL BINT4(x,y,ndata,ibcl,ibcr,fbcl,fbcr,knt,t,bc,n,k,w)
      !
      !       Error test on BINT4.
      !
      inbv = 1
      DO i = 1, ndata
        xx = x(i)
        bv = BVALU(t,bc,n,k,0,xx)
        er = ABS(y(i)-bv)
        IF( er>tol ) THEN
          Ipass = 0
          IF( Kprint>=2 ) WRITE (Lun,99002)
          99002 FORMAT (' ERROR TEST FOR INTERPOLATION BY BINT4 NOT SATISFIED')
        END IF
      END DO
      inbv = 1
      bv = BVALU(t,bc,n,k,1,x(1))
      er = ABS(pi-bv)
      IF( er>tol ) THEN
        Ipass = 0
        IF( Kprint>=2 ) WRITE (Lun,99003)
        99003 FORMAT (' ERROR TEST FOR INTERPOLATION BY BINT4 NOT SATISFIED ',&
          'BY FIRST DERIVATIVE')
      END IF
      bv = BVALU(t,bc,n,k,2,x(ndata))
      er = ABS(bv)
      IF( er>tol ) THEN
        Ipass = 0
        IF( Kprint>=2 ) WRITE (Lun,99004)
        99004 FORMAT (' ERROR TEST FOR INTERPOLATION BY BINT4 NOT SATISFIED ',&
          'BY SECOND DERIVATIVE')
      END IF
      !
      !       Test for equality of area from 4 routines.
      !
      x1 = x(1)
      x2 = x(ndata)
      CALL BSQAD(t,bc,n,k,x1,x2,bquad)
      ldc = 4
      CALL BSPPP(t,bc,n,k,ldc,c,xi,lxi,w)
      CALL PPQAD(ldc,c,xi,lxi,k,x1,x2,q(1))
      CALL BFQAD(FB,t,bc,n,k,0,x1,x2,tol,q(2),ierr)
      CALL PFQAD(FB,ldc,c,xi,lxi,k,0,x1,x2,tol,q(3),ierr)
      !
      !       Error test for quadratures.
      !
      DO i = 1, 3
        er = ABS(bquad-q(i))
        IF( er>tol ) THEN
          Ipass = 0
          IF( Kprint>=2 ) WRITE (Lun,99005)
          99005 FORMAT (' ERROR IN QUADRATURE CHECKS')
        END IF
      END DO
      qsave(knt) = bquad
    END DO
    er = ABS(qsave(1)-qsave(2))
    IF( er>tol ) THEN
      Ipass = 0
      IF( Kprint>=2 ) WRITE (Lun,99006)
      99006 FORMAT (' ERROR IN QUADRATURE CHECK USING TWO SETS OF KNOTS')
    END IF
    !
    !     Check BSPDR and BSPEV against BVALU, PPVAL and BSPVD.
    !
    CALL BSPDR(t,bc,n,k,k,adif)
    inev = 1
    inbv = 1
    inppv = 1
    ilo = 1
    DO i = 1, 6
      xx = x(i+i-1)
      CALL BSPEV(t,adif,n,k,k,xx,inev,sv,w)
      atol = tol
      DO j = 1, k
        spv = BVALU(t,bc,n,k,j-1,xx)
        er = ABS(spv-sv(j))
        x2 = ABS(sv(j))
        IF( x2>1._SP ) er = er/x2
        IF( er>atol ) THEN
          Ipass = 0
          IF( Kprint>=2 ) WRITE (Lun,99007)
          99007 FORMAT (' COMPARISONS FROM BSPEV AND BVALU DO NOT AGREE')
        END IF
        atol = 10._SP*atol
      END DO
      atol = tol
      DO j = 1, k
        spv = PPVAL(ldc,c,xi,lxi,k,j-1,xx)
        er = ABS(spv-sv(j))
        x2 = ABS(sv(j))
        IF( x2>1._SP ) er = er/x2
        IF( er>atol ) THEN
          Ipass = 0
          IF( Kprint>=2 ) WRITE (Lun,99008)
          99008 FORMAT (' COMPARISONS FROM BSPEV AND PPVAL DO NOT AGREE')
        END IF
        atol = 10._SP*atol
      END DO
      atol = tol
      ldcc = 4
      x1 = xx
      IF( i+i-1==ndata ) x1 = t(n)
      nn = n + k
      CALL INTRV(t,nn,x1,ilo,ileft,mflag)
      DO j = 1, k
        CALL BSPVD(t,k,j,xx,ileft,ldcc,cc,w)
        er = 0._SP
        DO jj = 1, k
          er = er + bc(ileft-k+jj)*cc(jj,j)
        END DO
        er = ABS(er-sv(j))
        x2 = ABS(sv(j))
        IF( x2>1._SP ) er = er/x2
        IF( er>atol ) THEN
          Ipass = 0
          IF( Kprint>=2 ) WRITE (Lun,99009)
          99009 FORMAT (' COMPARISONS FROM BSPEV AND BSPVD DO NOT AGREE')
        END IF
        atol = 10._SP*atol
      END DO
    END DO
    DO k = 2, 4
      n = ndata
      nmk = n - k
      DO i = 1, k
        t(i) = x(1)
        t(n+i) = x(n)
      END DO
      xl = x(n) - x(1)
      dn = n - k + 1
      DO i = 1, nmk
        t(k+i) = x(1) + i*xl/dn
      END DO
      CALL BINTK(x,y,t,n,k,bc,qq,w)
      !
      !       Error test on BINTK.
      !
      inbv = 1
      DO i = 1, n
        xx = x(i)
        bv = BVALU(t,bc,n,k,0,xx)
        er = ABS(y(i)-bv)
        IF( er>tol ) THEN
          Ipass = 0
          IF( Kprint>=2 ) WRITE (Lun,99010)
          99010 FORMAT (' ERROR TEST FOR INTERPOLATION BY BINTK NOT SATISFIED')
        END IF
      END DO
    END DO
    !
    !     Trigger error conditions.
    !
!    kontrl = control_xer
!    IF( Kprint<=2 ) THEN
!      control_xer = 0
!    ELSE
!      control_xer = 1
!    END IF
!    fatal = .FALSE.
!    num_xer = 0
!    !
!    IF( Kprint>=3 ) WRITE (Lun,99011)
!    99011 FORMAT (/' TRIGGER 52 ERROR CONDITIONS',/)
!    !
!    w(1) = 11._SP
!    w(2) = 4._SP
!    w(3) = 2._SP
!    w(4) = 0.5_SP
!    w(5) = 4._SP
!    ilo = 1
!    inev = 1
!    inbv = 1
!    CALL INTRV(t,n+1,w(4),ilo,ileft,mflag)
!    DO i = 1, 5
!      w(i) = -w(i)
!      n = INT( w(1) )
!      k = INT( w(2) )
!      id = INT( w(3) )
!      xx = w(4)
!      ldc = INT( w(5) )
!      IF( i<=4 ) THEN
!        bv = BVALU(t,bc,n,k,id,xx,inbv,qq)
!        IF( num_xer/=2 ) THEN
!          Ipass = 0
!          fatal = .TRUE.
!        END IF
!        num_xer = 0
!        !
!        CALL BSPEV(t,adif,n,k,id,xx,inev,sv,qq)
!        IF( num_xer/=2 ) THEN
!          Ipass = 0
!          fatal = .TRUE.
!        END IF
!        num_xer = 0
!        !
!        jhigh = n - 10
!        CALL BSPVN(t,jhigh,k,id,xx,ileft,sv,qq,iwork)
!        IF( num_xer/=2 ) THEN
!          Ipass = 0
!          fatal = .TRUE.
!        END IF
!        num_xer = 0
!        !
!        CALL BFQAD(FB,t,bc,n,k,id,xx,x2,tol,quad,ierr,qq)
!        IF( num_xer/=2 ) THEN
!          Ipass = 0
!          fatal = .TRUE.
!        END IF
!        num_xer = 0
!      END IF
!      !
!      IF( i/=3 .AND. i/=4 ) THEN
!        CALL BSPPP(t,bc,n,k,ldc,c,xi,lxi,qq)
!        IF( num_xer/=2 ) THEN
!          Ipass = 0
!          fatal = .TRUE.
!        END IF
!        num_xer = 0
!      END IF
!      !
!      IF( i<=3 ) THEN
!        CALL BSPDR(t,bc,n,k,id,adif)
!        IF( num_xer/=2 ) THEN
!          Ipass = 0
!          fatal = .TRUE.
!        END IF
!        num_xer = 0
!      END IF
!      !
!      IF( i/=3 .AND. i/=5 ) THEN
!        CALL BSQAD(t,bc,n,k,xx,x2,bquad,qq)
!        IF( num_xer/=2 ) THEN
!          Ipass = 0
!          fatal = .TRUE.
!        END IF
!        num_xer = 0
!      END IF
!      !
!      IF( i>1 ) THEN
!        CALL BSPVD(t,k,id,xx,ileft,ldc,c,qq)
!        IF( num_xer/=2 ) THEN
!          Ipass = 0
!          fatal = .TRUE.
!        END IF
!        num_xer = 0
!      END IF
!      !
!      IF( i<=2 ) THEN
!        CALL BINTK(x,y,t,n,k,bc,qq,adif)
!        IF( num_xer/=2 ) THEN
!          Ipass = 0
!          fatal = .TRUE.
!        END IF
!        num_xer = 0
!      END IF
!      !
!      IF( i/=4 ) THEN
!        kntopt = ldc - 2
!        ibcl = k - 2
!        CALL BINT4(x,y,n,ibcl,id,fbcl,fbcr,kntopt,t,bc,nn,kk,qq)
!        IF( num_xer/=2 ) THEN
!          Ipass = 0
!          fatal = .TRUE.
!        END IF
!        num_xer = 0
!      END IF
!      w(i) = -w(i)
!    END DO
!    kntopt = 1
!    x(1) = 1._SP
!    CALL BINT4(x,y,n,ibcl,ibcr,fbcl,fbcr,kntopt,t,bc,n,k,qq)
!    IF( num_xer/=2 ) THEN
!      Ipass = 0
!      fatal = .TRUE.
!    END IF
!    num_xer = 0
!    !
!    CALL BINTK(x,y,t,n,k,bc,qq,adif)
!    IF( num_xer/=2 ) THEN
!      Ipass = 0
!      fatal = .TRUE.
!    END IF
!    num_xer = 0
!    !
!    x(1) = 0._SP
!    atol = 1._SP
!    kntopt = 3
!    DO i = 1, 3
!      qq(i) = -0.30_SP + 0.10_SP*(i-1)
!      qq(i+3) = 1.1_SP + 0.10_SP*(i-1)
!    END DO
!    qq(1) = 1._SP
!    CALL BINT4(x,y,ndata,1,1,fbcl,fbcr,3,t,bc,n,k,qq)
!    IF( num_xer/=2 ) THEN
!      Ipass = 0
!      fatal = .TRUE.
!    END IF
!    num_xer = 0
!    !
!    CALL BFQAD(FB,t,bc,n,k,id,x1,x2,atol,quad,ierr,qq)
!    IF( num_xer/=2 ) THEN
!      Ipass = 0
!      fatal = .TRUE.
!    END IF
!    num_xer = 0
!    !
!    inppv = 1
!    DO i = 1, 5
!      w(i) = -w(i)
!      lxi = INT( w(1) )
!      k = INT( w(2) )
!      id = INT( w(3) )
!      xx = w(4)
!      ldc = INT( w(5) )
!      spv = PPVAL(ldc,c,xi,lxi,k,id,xx)
!      IF( (i/=4 .AND. num_xer/=2) .OR. (i==4 .AND. num_xer/=0) ) THEN
!        Ipass = 0
!        fatal = .TRUE.
!      END IF
!      num_xer = 0
!      !
!      CALL PFQAD(FB,ldc,c,xi,lxi,k,id,xx,x2,tol,quad,ierr)
!      IF( (i/=4 .AND. num_xer/=2) .OR. (i==4 .AND. num_xer/=0) ) THEN
!        Ipass = 0
!        fatal = .TRUE.
!      END IF
!      num_xer = 0
!      !
!      IF( i/=3 ) THEN
!        CALL PPQAD(ldc,c,xi,lxi,k,xx,x2,pquad)
!        IF( (i/=4 .AND. num_xer/=2) .OR. (i==4 .AND. num_xer/=0) ) THEN
!          Ipass = 0
!          fatal = .TRUE.
!        END IF
!        num_xer = 0
!      END IF
!      !
!      w(i) = -w(i)
!    END DO
!    ldc = INT( w(5) )
!    CALL PFQAD(FB,ldc,c,xi,lxi,k,id,x1,x2,atol,quad,ierr)
!    IF( num_xer/=2 ) THEN
!      Ipass = 0
!      fatal = .TRUE.
!    END IF
!    num_xer = 0
    !
    !     Restore KONTRL and check to see if the tests of error detection
    !     passed.
    !
!    control_xer = kontrl
!    IF( fatal ) THEN
!      IF( Kprint>=2 ) THEN
!        WRITE (Lun,99012)
!        99012 FORMAT (/' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
!      END IF
!    ELSEIF( Kprint>=3 ) THEN
!      WRITE (Lun,99013)
!      99013 FORMAT (/' ALL INCORRECT ARGUMENT TESTS PASSED')
!    END IF
    !
    !     Print PASS/FAIL message.
    !
    IF( Ipass==1 .AND. Kprint>=2 ) WRITE (Lun,99014)
    99014 FORMAT (/' **********B-SPLINE PACKAGE PASSED ALL TESTS**********')
    IF( Ipass==0 .AND. Kprint>=1 ) WRITE (Lun,99015)
    99015 FORMAT (/' *********B-SPLINE PACKAGE FAILED SOME TESTS**********')
    RETURN
  END SUBROUTINE BSPCK
  !** FB
  PURE REAL(SP) FUNCTION FB(X)
    !> Subsidiary to BSPCK.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      SINGLE PRECISION (FB-S, DFB-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   930214  Added TYPE statement.  (WRB)

    REAL(SP), INTENT(IN) :: X
    !* FIRST EXECUTABLE STATEMENT  FB
    FB = 1._SP
  END FUNCTION FB
END MODULE TEST30_MOD
!** TEST30
PROGRAM TEST30
  USE TEST30_MOD, ONLY : BSPCK
  USE slatec, ONLY : I1MACH, control_xer, max_xer
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !> Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  E, E1A, E3
  !***
  ! **Type:**      SINGLE PRECISION (TEST30-S, TEST31-D)
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
  !        BFQAD    BINT4    BINTK    BSPDR    BSPEV    BSPPP
  !        BSPVD    BSPVN    BSQAD    BVALU    INTRV    PFQAD
  !        PPQAD    PPVAL
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  BSPCK, I1MACH, XERMAX, XSETF

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  INTEGER :: ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST30
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
  !     Test single precision B-Spline package
  !
  CALL BSPCK(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST30 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST30 *************')
  END IF
  STOP
END PROGRAM TEST30
