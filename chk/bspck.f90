!DECK BSPCK
SUBROUTINE BSPCK(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  BSPCK
  !***PURPOSE  Quick check for the B-Spline package.
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (BSPCK-S, DBSPCK-D)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !   BSPCK is a quick check routine for the B-Spline package which
  !   tests consistency between results from higher level routines.
  !   Those routines not explicitly called are exercised at some lower
  !   level.  The routines exercised are BFQAD, BINT4, BINTK, BNFAC,
  !   BNSLV, BSGQ8, BSPDR, BSPEV, BSPPP, BSPVD, BSPVN, BSQAD, BVALU,
  !   INTRV, PFQAD, PPGQ8, PPQAD and PPVAL.
  !
  !***ROUTINES CALLED  BFQAD, BINT4, BINTK, BSPDR, BSPEV, BSPPP, BSPVD,
  !                    BSPVN, BSQAD, BVALU, FB, INTRV, PFQAD, PPQAD,
  !                    PPVAL, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891004  Removed unreachable code.  (WRB)
  !   891009  Removed unreferenced variables.  (WRB)
  !   891009  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   930214  Declarations sections added, code revised to test error
  !           returns for all values of KPRINT and code polished.  (WRB)
  !***END PROLOGUE  BSPCK
  !     .. Scalar Arguments ..
  INTEGER Ipass, Kprint, Lun
  !     .. Local Scalars ..
  REAL atol, bquad, bv, den, dn, er, fbcl, fbcr, pi, pquad, quad, &
    spv, tol, x1, x2, xl, xx
  INTEGER i, ibcl, ibcr, id, ierr, iknt, ileft, ilo, inbv, inev, &
    inppv, iwork, j, jhigh, jj, k, kk, knt, kntopt, kontrl, &
    ldc, ldcc, lxi, mflag, n, ndata, nerr, nmk, nn
  LOGICAL fatal
  !     .. Local Arrays ..
  REAL adif(52), bc(13), c(4,10), cc(4,4), q(3), qq(77), qsave(2), &
    sv(4), t(17), w(65), x(11), xi(11), y(11)
  !     .. External Functions ..
  REAL BVALU, FB, PPVAL, R1MACH
  INTEGER NUMXER
  EXTERNAL BVALU, FB, NUMXER, PPVAL, R1MACH
  !     .. External Subroutines ..
  EXTERNAL BFQAD, BINT4, BINTK, BSPDR, BSPEV, BSPPP, BSPVD, BSPVN, &
    BSQAD, INTRV, PFQAD, PPQAD, XGETF, XSETF
  !     .. Intrinsic Functions ..
  INTRINSIC ABS, SIN
  !***FIRST EXECUTABLE STATEMENT  BSPCK
  IF ( Kprint>=2 ) WRITE (Lun,99001)
  !
  99001 FORMAT ('1 QUICK CHECK FOR SPLINE ROUTINES',//)
  !
  Ipass = 1
  pi = 3.14159265358979324E0
  tol = 1000.0E0*R1MACH(4)
  !
  !     Generate data.
  !
  ndata = 11
  den = ndata - 1
  DO i = 1, ndata
    x(i) = (i-1)/den
    y(i) = SIN(pi*x(i))
  ENDDO
  x(3) = 2.0E0/den
  y(3) = SIN(pi*x(3))
  !
  !     Compute splines for two knot arrays.
  !
  DO iknt = 1, 2
    knt = 3 - iknt
    ibcl = 1
    ibcr = 2
    fbcl = pi
    fbcr = 0.0E0
    CALL BINT4(x,y,ndata,ibcl,ibcr,fbcl,fbcr,knt,t,bc,n,k,w)
    !
    !       Error test on BINT4.
    !
    inbv = 1
    DO i = 1, ndata
      xx = x(i)
      bv = BVALU(t,bc,n,k,0,xx,inbv,w)
      er = ABS(y(i)-bv)
      IF ( er>tol ) THEN
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Lun,99002)
        99002       FORMAT (' ERROR TEST FOR INTERPOLATION BY BINT4 NOT SATISFIED')
      ENDIF
    ENDDO
    inbv = 1
    bv = BVALU(t,bc,n,k,1,x(1),inbv,w)
    er = ABS(pi-bv)
    IF ( er>tol ) THEN
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (Lun,99003)
      99003     FORMAT (' ERROR TEST FOR INTERPOLATION BY BINT4 NOT SATISFIED ',&
        'BY FIRST DERIVATIVE')
    ENDIF
    bv = BVALU(t,bc,n,k,2,x(ndata),inbv,w)
    er = ABS(bv)
    IF ( er>tol ) THEN
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (Lun,99004)
      99004     FORMAT (' ERROR TEST FOR INTERPOLATION BY BINT4 NOT SATISFIED ',&
        'BY SECOND DERIVATIVE')
    ENDIF
    !
    !       Test for equality of area from 4 routines.
    !
    x1 = x(1)
    x2 = x(ndata)
    CALL BSQAD(t,bc,n,k,x1,x2,bquad,w)
    ldc = 4
    CALL BSPPP(t,bc,n,k,ldc,c,xi,lxi,w)
    CALL PPQAD(ldc,c,xi,lxi,k,x1,x2,q(1))
    CALL BFQAD(FB,t,bc,n,k,0,x1,x2,tol,q(2),ierr,w)
    CALL PFQAD(FB,ldc,c,xi,lxi,k,0,x1,x2,tol,q(3),ierr)
    !
    !       Error test for quadratures.
    !
    DO i = 1, 3
      er = ABS(bquad-q(i))
      IF ( er>tol ) THEN
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Lun,99005)
        99005       FORMAT (' ERROR IN QUADRATURE CHECKS')
      ENDIF
    ENDDO
    qsave(knt) = bquad
  ENDDO
  er = ABS(qsave(1)-qsave(2))
  IF ( er>tol ) THEN
    Ipass = 0
    IF ( Kprint>=2 ) WRITE (Lun,99006)
    99006   FORMAT (' ERROR IN QUADRATURE CHECK USING TWO SETS OF KNOTS')
  ENDIF
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
      spv = BVALU(t,bc,n,k,j-1,xx,inbv,w)
      er = ABS(spv-sv(j))
      x2 = ABS(sv(j))
      IF ( x2>1.0E0 ) er = er/x2
      IF ( er>atol ) THEN
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Lun,99007)
        99007       FORMAT (' COMPARISONS FROM BSPEV AND BVALU DO NOT AGREE')
      ENDIF
      atol = 10.0E0*atol
    ENDDO
    atol = tol
    DO j = 1, k
      spv = PPVAL(ldc,c,xi,lxi,k,j-1,xx,inppv)
      er = ABS(spv-sv(j))
      x2 = ABS(sv(j))
      IF ( x2>1.0E0 ) er = er/x2
      IF ( er>atol ) THEN
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Lun,99008)
        99008       FORMAT (' COMPARISONS FROM BSPEV AND PPVAL DO NOT AGREE')
      ENDIF
      atol = 10.0E0*atol
    ENDDO
    atol = tol
    ldcc = 4
    x1 = xx
    IF ( i+i-1==ndata ) x1 = t(n)
    nn = n + k
    CALL INTRV(t,nn,x1,ilo,ileft,mflag)
    DO j = 1, k
      CALL BSPVD(t,k,j,xx,ileft,ldcc,cc,w)
      er = 0.0E0
      DO jj = 1, k
        er = er + bc(ileft-k+jj)*cc(jj,j)
      ENDDO
      er = ABS(er-sv(j))
      x2 = ABS(sv(j))
      IF ( x2>1.0E0 ) er = er/x2
      IF ( er>atol ) THEN
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Lun,99009)
        99009       FORMAT (' COMPARISONS FROM BSPEV AND BSPVD DO NOT AGREE')
      ENDIF
      atol = 10.0E0*atol
    ENDDO
  ENDDO
  DO k = 2, 4
    n = ndata
    nmk = n - k
    DO i = 1, k
      t(i) = x(1)
      t(n+i) = x(n)
    ENDDO
    xl = x(n) - x(1)
    dn = n - k + 1
    DO i = 1, nmk
      t(k+i) = x(1) + i*xl/dn
    ENDDO
    CALL BINTK(x,y,t,n,k,bc,qq,w)
    !
    !       Error test on BINTK.
    !
    inbv = 1
    DO i = 1, n
      xx = x(i)
      bv = BVALU(t,bc,n,k,0,xx,inbv,w)
      er = ABS(y(i)-bv)
      IF ( er>tol ) THEN
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Lun,99010)
        99010       FORMAT (' ERROR TEST FOR INTERPOLATION BY BINTK NOT SATISFIED')
      ENDIF
    ENDDO
  ENDDO
  !
  !     Trigger error conditions.
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
  IF ( Kprint>=3 ) WRITE (Lun,99011)
  99011 FORMAT (/' TRIGGER 52 ERROR CONDITIONS',/)
  !
  w(1) = 11.0E0
  w(2) = 4.0E0
  w(3) = 2.0E0
  w(4) = 0.5E0
  w(5) = 4.0E0
  ilo = 1
  inev = 1
  inbv = 1
  CALL INTRV(t,n+1,w(4),ilo,ileft,mflag)
  DO i = 1, 5
    w(i) = -w(i)
    n = w(1)
    k = w(2)
    id = w(3)
    xx = w(4)
    ldc = w(5)
    IF ( i<=4 ) THEN
      bv = BVALU(t,bc,n,k,id,xx,inbv,qq)
      IF ( NUMXER(nerr)/=2 ) THEN
        Ipass = 0
        fatal = .TRUE.
      ENDIF
      CALL XERCLR
      !
      CALL BSPEV(t,adif,n,k,id,xx,inev,sv,qq)
      IF ( NUMXER(nerr)/=2 ) THEN
        Ipass = 0
        fatal = .TRUE.
      ENDIF
      CALL XERCLR
      !
      jhigh = n - 10
      CALL BSPVN(t,jhigh,k,id,xx,ileft,sv,qq,iwork)
      IF ( NUMXER(nerr)/=2 ) THEN
        Ipass = 0
        fatal = .TRUE.
      ENDIF
      CALL XERCLR
      !
      CALL BFQAD(FB,t,bc,n,k,id,xx,x2,tol,quad,ierr,qq)
      IF ( NUMXER(nerr)/=2 ) THEN
        Ipass = 0
        fatal = .TRUE.
      ENDIF
      CALL XERCLR
    ENDIF
    !
    IF ( i/=3.AND.i/=4 ) THEN
      CALL BSPPP(t,bc,n,k,ldc,c,xi,lxi,qq)
      IF ( NUMXER(nerr)/=2 ) THEN
        Ipass = 0
        fatal = .TRUE.
      ENDIF
      CALL XERCLR
    ENDIF
    !
    IF ( i<=3 ) THEN
      CALL BSPDR(t,bc,n,k,id,adif)
      IF ( NUMXER(nerr)/=2 ) THEN
        Ipass = 0
        fatal = .TRUE.
      ENDIF
      CALL XERCLR
    ENDIF
    !
    IF ( i/=3.AND.i/=5 ) THEN
      CALL BSQAD(t,bc,n,k,xx,x2,bquad,qq)
      IF ( NUMXER(nerr)/=2 ) THEN
        Ipass = 0
        fatal = .TRUE.
      ENDIF
      CALL XERCLR
    ENDIF
    !
    IF ( i>1 ) THEN
      CALL BSPVD(t,k,id,xx,ileft,ldc,c,qq)
      IF ( NUMXER(nerr)/=2 ) THEN
        Ipass = 0
        fatal = .TRUE.
      ENDIF
      CALL XERCLR
    ENDIF
    !
    IF ( i<=2 ) THEN
      CALL BINTK(x,y,t,n,k,bc,qq,adif)
      IF ( NUMXER(nerr)/=2 ) THEN
        Ipass = 0
        fatal = .TRUE.
      ENDIF
      CALL XERCLR
    ENDIF
    !
    IF ( i/=4 ) THEN
      kntopt = ldc - 2
      ibcl = k - 2
      CALL BINT4(x,y,n,ibcl,id,fbcl,fbcr,kntopt,t,bc,nn,kk,qq)
      IF ( NUMXER(nerr)/=2 ) THEN
        Ipass = 0
        fatal = .TRUE.
      ENDIF
      CALL XERCLR
    ENDIF
    w(i) = -w(i)
  ENDDO
  kntopt = 1
  x(1) = 1.0E0
  CALL BINT4(x,y,n,ibcl,ibcr,fbcl,fbcr,kntopt,t,bc,n,k,qq)
  IF ( NUMXER(nerr)/=2 ) THEN
    Ipass = 0
    fatal = .TRUE.
  ENDIF
  CALL XERCLR
  !
  CALL BINTK(x,y,t,n,k,bc,qq,adif)
  IF ( NUMXER(nerr)/=2 ) THEN
    Ipass = 0
    fatal = .TRUE.
  ENDIF
  CALL XERCLR
  !
  x(1) = 0.0E0
  atol = 1.0E0
  kntopt = 3
  DO i = 1, 3
    qq(i) = -0.30E0 + 0.10E0*(i-1)
    qq(i+3) = 1.1E0 + 0.10E0*(i-1)
  ENDDO
  qq(1) = 1.0E0
  CALL BINT4(x,y,ndata,1,1,fbcl,fbcr,3,t,bc,n,k,qq)
  IF ( NUMXER(nerr)/=2 ) THEN
    Ipass = 0
    fatal = .TRUE.
  ENDIF
  CALL XERCLR
  !
  CALL BFQAD(FB,t,bc,n,k,id,x1,x2,atol,quad,ierr,qq)
  IF ( NUMXER(nerr)/=2 ) THEN
    Ipass = 0
    fatal = .TRUE.
  ENDIF
  CALL XERCLR
  !
  inppv = 1
  DO i = 1, 5
    w(i) = -w(i)
    lxi = w(1)
    k = w(2)
    id = w(3)
    xx = w(4)
    ldc = w(5)
    spv = PPVAL(ldc,c,xi,lxi,k,id,xx,inppv)
    IF ( (i/=4.AND.NUMXER(nerr)/=2).OR.(i==4.AND.NUMXER(nerr)/=0) ) THEN
      Ipass = 0
      fatal = .TRUE.
    ENDIF
    CALL XERCLR
    !
    CALL PFQAD(FB,ldc,c,xi,lxi,k,id,xx,x2,tol,quad,ierr)
    IF ( (i/=4.AND.NUMXER(nerr)/=2).OR.(i==4.AND.NUMXER(nerr)/=0) ) THEN
      Ipass = 0
      fatal = .TRUE.
    ENDIF
    CALL XERCLR
    !
    IF ( i/=3 ) THEN
      CALL PPQAD(ldc,c,xi,lxi,k,xx,x2,pquad)
      IF ( (i/=4.AND.NUMXER(nerr)/=2).OR.(i==4.AND.NUMXER(nerr)/=0) ) THEN
        Ipass = 0
        fatal = .TRUE.
      ENDIF
      CALL XERCLR
    ENDIF
    !
    w(i) = -w(i)
  ENDDO
  ldc = w(5)
  CALL PFQAD(FB,ldc,c,xi,lxi,k,id,x1,x2,atol,quad,ierr)
  IF ( NUMXER(nerr)/=2 ) THEN
    Ipass = 0
    fatal = .TRUE.
  ENDIF
  CALL XERCLR
  !
  !     Restore KONTRL and check to see if the tests of error detection
  !     passed.
  !
  CALL XSETF(kontrl)
  IF ( fatal ) THEN
    IF ( Kprint>=2 ) THEN
      WRITE (Lun,99012)
      99012     FORMAT (/' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
    ENDIF
  ELSEIF ( Kprint>=3 ) THEN
    WRITE (Lun,99013)
    99013   FORMAT (/' ALL INCORRECT ARGUMENT TESTS PASSED')
  ENDIF
  !
  !     Print PASS/FAIL message.
  !
  IF ( Ipass==1.AND.Kprint>=2 ) WRITE (Lun,99014)
  99014 FORMAT (/' **********B-SPLINE PACKAGE PASSED ALL TESTS**********')
  IF ( Ipass==0.AND.Kprint>=1 ) WRITE (Lun,99015)
  99015 FORMAT (/' *********B-SPLINE PACKAGE FAILED SOME TESTS**********')
  RETURN
END SUBROUTINE BSPCK
