MODULE TEST43_MOD
  IMPLICIT NONE

CONTAINS
  !DECK JAC
  SUBROUTINE JAC(T,U,Pd,Nrowpd,Rpar,Ipar)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  JAC
    !***SUBSIDIARY
    !***PURPOSE  Evaluate Jacobian for DEBDF quick check.
    !***LIBRARY   SLATEC
    !***TYPE      SINGLE PRECISION (JAC-S, DJAC-D)
    !***AUTHOR  Chow, Jeff (LANL)
    !***ROUTINES CALLED  (NONE)
    !***REVISION HISTORY  (YYMMDD)
    !   810801  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   900415  Minor clean-up of prologue and code.  (WRB)
    !***END PROLOGUE  JAC
    INTEGER Ipar, Nrowpd
    REAL Pd, r, r5, Rpar, rsq, T, U, u1sq, u2sq, u1u2
    DIMENSION U(*), Pd(Nrowpd,*), Rpar(*), Ipar(*)
    !***FIRST EXECUTABLE STATEMENT  JAC
    u1sq = U(1)*U(1)
    u2sq = U(2)*U(2)
    u1u2 = U(1)*U(2)
    rsq = u1sq + u2sq
    r = SQRT(rsq)
    r5 = rsq*rsq*r
    Pd(3,1) = (3.E0*u1sq-rsq)/r5
    Pd(4,1) = 3.E0*u1u2/r5
    Pd(3,2) = Pd(4,1)
    Pd(4,2) = (3.E0*u2sq-rsq)/r5
    Pd(1,3) = 1.E0
    Pd(2,4) = 1.E0
  END SUBROUTINE JAC
  !DECK FDEQC
  SUBROUTINE FDEQC(T,U,Uprime,Rpar,Ipar)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  FDEQC
    !***SUBSIDIARY
    !***PURPOSE  Derivative evaluator for DEPAC quick checks.
    !***LIBRARY   SLATEC
    !***TYPE      SINGLE PRECISION (FDEQC-S, DFDEQC-D)
    !***AUTHOR  Chow, Jeff, (LANL)
    !***ROUTINES CALLED  (NONE)
    !***REVISION HISTORY  (YYMMDD)
    !   810801  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   900415  Name changed from F to FDEQC.  (WRB)
    !***END PROLOGUE  FDEQC
    !
    !     Declare arguments.
    !
    INTEGER Ipar(*)
    REAL Rpar(*), T, U(*), Uprime(*)
    !
    !     Declare local variables.
    !
    REAL r, rsq, r3
    !***FIRST EXECUTABLE STATEMENT  FDEQC
    rsq = U(1)*U(1) + U(2)*U(2)
    r = SQRT(rsq)
    r3 = rsq*r
    Uprime(1) = U(3)
    Uprime(2) = U(4)
    Uprime(3) = -(U(1)/r3)
    Uprime(4) = -(U(2)/r3)
  END SUBROUTINE FDEQC
  !DECK QXABM
  SUBROUTINE QXABM(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  QXABM
    !***SUBSIDIARY
    !***PURPOSE  Test the DEPAC routine DEABM.
    !***LIBRARY   SLATEC
    !***TYPE      SINGLE PRECISION (QXABM-S, QXDABM-D)
    !***KEYWORDS  QUICK CHECK
    !***AUTHOR  Chow, Jeff, (LANL)
    !***DESCRIPTION
    !
    ! *Usage:
    !
    !        INTEGER  LUN, KPRINT, IPASS
    !
    !        CALL QXABM (LUN, KPRINT, IPASS)
    !
    ! *Arguments:
    !
    !     LUN   :IN  is the unit number to which output is to be written.
    !
    !     KPRINT:IN  controls the amount of output, as specified in the
    !                SLATEC Guidelines.
    !
    !     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
    !                IPASS=0 indicates one or more tests failed.
    !
    ! *Description:
    !
    !   DEABM is tested by solving the equations of motion of a body
    !   moving in a plane about a spherical earth, namely
    !           (D/DT)(D/DT)X = -G*X/R**3
    !           (D/DT)(D/DT)Y = -G*Y/R**3
    !   where G = 1, R = SQRT(X**2 + Y**2) and
    !           X(0) = 1
    !           (D/DT)X(0) = 0
    !           Y(0) = 0
    !           (D/DT)Y(0) = 1.
    !
    !***ROUTINES CALLED  DEABM, FDEQC, R1MACH
    !***REVISION HISTORY  (YYMMDD)
    !   810801  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   900415  Code extensively revised.  (WRB)
    !***END PROLOGUE  QXABM
    !
    !     Declare arguments.
    !
    INTEGER Lun, Kprint, Ipass
    !
    !     Declare local variables.
    !
    INTEGER idid, info(15), ipar, iwork(51), n, liw, lrw
    REAL abserr, r, R1MACH, relerr, reltol, rpar, rwork(214), t, tout, u(4)
    !***FIRST EXECUTABLE STATEMENT  QXABM
    IF ( Kprint>=2 ) WRITE (Lun,99001)
    !
    ! FORMATs.
    !
    99001 FORMAT ('1'/' ------------  DEABM QUICK CHECK OUTPUT',' ------------')
    !
    !     Initialize problem.
    !
    n = 4
    lrw = 214
    liw = 51
    t = 0.0E0
    tout = 8.0E0*ATAN(1.0E0)
    u(1) = 1.0E0
    u(2) = 0.0E0
    u(3) = 0.0E0
    u(4) = 1.0E0
    Ipass = 1
    reltol = SQRT(R1MACH(4))
    relerr = 0.1E0*reltol
    abserr = relerr**1.5E0
    info(1) = 0
    info(2) = 0
    info(3) = 1
    info(4) = 0
    IF ( Kprint>2 ) WRITE (Lun,99002) relerr, abserr, t, (1.0E0)
    99002 FORMAT (/' RELERR = ',E16.8,'   ABSERR =',E16.8/12X,'T',19X,'R'/2E20.8)
    DO
      !
      CALL DEABM(FDEQC,n,t,u,tout,info,relerr,abserr,idid,rwork,lrw,iwork,liw,&
        rpar,ipar)
      r = SQRT(u(1)*u(1)+u(2)*u(2))
      IF ( ABS(r-1.0E0)>reltol ) Ipass = 0
      IF ( Kprint>2 ) WRITE (Lun,99003) t, r
      99003 FORMAT (2E20.8)
      info(1) = 1
      IF ( idid/=1 ) THEN
        !
        !     Finish up.
        !
        IF ( idid<1 ) Ipass = 0
        IF ( Kprint>1.AND.idid<1 ) WRITE (Lun,99004) idid
        99004 FORMAT (1X,'ERROR RETURN FROM DEABM.  IDID = ',I3)
        IF ( Kprint>1.AND.Ipass==1 ) WRITE (Lun,99005)
        99005 FORMAT (/' ------------  DEABM PASSED TESTS  ------------')
        IF ( Kprint>=1.AND.Ipass==0 ) WRITE (Lun,99006)
        99006 FORMAT (/' ************  DEABM FAILED TESTS  ************')
        RETURN
      ENDIF
    ENDDO
  END SUBROUTINE QXABM
  !DECK QXBDF
  SUBROUTINE QXBDF(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  QXBDF
    !***PURPOSE  Test the DEPAC routine DEBDF.
    !***LIBRARY   SLATEC
    !***TYPE      SINGLE PRECISION (QXBDF-S, QXDBDF-D)
    !***KEYWORDS  QUICK CHECK
    !***AUTHOR  Chow, Jeff, (LANL)
    !***DESCRIPTION
    !
    ! *Usage:
    !
    !        INTEGER  LUN, KPRINT, IPASS
    !
    !        CALL QXBDF (LUN, KPRINT, IPASS)
    !
    ! *Arguments:
    !
    !     LUN   :IN  is the unit number to which output is to be written.
    !
    !     KPRINT:IN  controls the amount of output, as specified in the
    !                SLATEC Guidelines.
    !
    !     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
    !                IPASS=0 indicates one or more tests failed.
    !
    ! *Description:
    !
    !   DEBDF is tested by solving the equations of motion of a body
    !   moving in a plane about a spherical earth, namely
    !           (D/DT)(D/DT)X = -G*X/R**3
    !           (D/DT)(D/DT)Y = -G*Y/R**3
    !   where G = 1, R = SQRT(X**2 + Y**2) and
    !           X(0) = 1
    !           (D/DT)X(0) = 0
    !           Y(0) = 0
    !           (D/DT)Y(0) = 1.
    !
    !***ROUTINES CALLED  DEBDF, FDEQC, JAC, R1MACH
    !***REVISION HISTORY  (YYMMDD)
    !   810801  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   900415  Code extensively revised.  (WRB)
    !***END PROLOGUE  QXBDF
    !
    !     Declare arguments.
    !
    INTEGER Lun, Kprint, Ipass
    !
    !     Declare local variables.
    !
    INTEGER idid, info(15), ipar, iwork(60), n, liw, lrw
    REAL abserr, r, R1MACH, relerr, reltol, rpar, rwork(306), t, tout, u(4)
    !***FIRST EXECUTABLE STATEMENT  QXBDF
    IF ( Kprint>=2 ) WRITE (Lun,99001)
    !
    ! FORMATs.
    !
    99001 FORMAT ('1'/' ------------  DEBDF QUICK CHECK OUTPUT',' ------------')
    !
    !     Initialize problem.
    !
    n = 4
    lrw = 306
    liw = 60
    t = 0.0E0
    tout = 8.0E0*ATAN(1.0E0)
    u(1) = 1.0E0
    u(2) = 0.0E0
    u(3) = 0.0E0
    u(4) = 1.0E0
    Ipass = 1
    reltol = SQRT(R1MACH(4))
    relerr = 0.001E0*reltol
    abserr = relerr**1.5E0
    info(1) = 0
    info(2) = 0
    info(3) = 1
    info(4) = 0
    info(5) = 1
    info(6) = 0
    IF ( Kprint>2 ) WRITE (Lun,99002) relerr, abserr, t, (1.0E0)
    99002 FORMAT (/' RELERR = ',E16.8,'   ABSERR =',E16.8/12X,'T',19X,'R'/2E20.8)
    DO
      !
      CALL DEBDF(FDEQC,n,t,u,tout,info,relerr,abserr,idid,rwork,lrw,iwork,liw,&
        rpar,ipar,JAC)
      r = SQRT(u(1)*u(1)+u(2)*u(2))
      IF ( ABS(r-1.0E0)>reltol ) Ipass = 0
      IF ( Kprint>2 ) WRITE (Lun,99003) t, r
      99003 FORMAT (2E20.8)
      info(1) = 1
      IF ( idid/=1 ) THEN
        !
        !     Finish up.
        !
        IF ( idid<1 ) Ipass = 0
        IF ( Kprint>1.AND.idid<1 ) WRITE (Lun,99004) idid
        99004 FORMAT (1X,'ERROR RETURN FROM DEBDF.  IDID = ',I3)
        IF ( Kprint>1.AND.Ipass==1 ) WRITE (Lun,99005)
        99005 FORMAT (/' ------------  DEBDF PASSED TESTS  ------------')
        IF ( Kprint>=1.AND.Ipass==0 ) WRITE (Lun,99006)
        99006 FORMAT (/' ************  DEBDF FAILED TESTS  ************')
        RETURN
      ENDIF
    ENDDO
  END SUBROUTINE QXBDF
  !DECK QXRKF
  SUBROUTINE QXRKF(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  QXRKF
    !***PURPOSE  Test the DEPAC routine DERKF.
    !***LIBRARY   SLATEC
    !***TYPE      SINGLE PRECISION (QXRKF-S, QXDRKF-D)
    !***KEYWORDS  QUICK CHECK
    !***AUTHOR  Chow, Jeff, (LANL)
    !***DESCRIPTION
    !
    ! *Usage:
    !
    !        INTEGER  LUN, KPRINT, IPASS
    !
    !        CALL QXRKF (LUN, KPRINT, IPASS)
    !
    ! *Arguments:
    !
    !     LUN   :IN  is the unit number to which output is to be written.
    !
    !     KPRINT:IN  controls the amount of output, as specified in the
    !                SLATEC Guidelines.
    !
    !     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
    !                IPASS=0 indicates one or more tests failed.
    !
    ! *Description:
    !
    !   DERKF is tested by solving the equations of motion of a body
    !   moving in a plane about a spherical earth, namely
    !           (D/DT)(D/DT)X = -G*X/R**3
    !           (D/DT)(D/DT)Y = -G*Y/R**3
    !   where G = 1, R = SQRT(X**2 + Y**2) and
    !           X(0) = 1
    !           (D/DT)X(0) = 0
    !           Y(0) = 0
    !           (D/DT)Y(0) = 1.
    !
    !***ROUTINES CALLED  DERKF, FDEQC, R1MACH
    !***REVISION HISTORY  (YYMMDD)
    !   810801  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   900415  Code extensively revised.  (WRB)
    !***END PROLOGUE  QXRKF
    !
    !     Declare arguments.
    !
    INTEGER Lun, Kprint, Ipass
    !
    !     Declare local variables.
    !
    INTEGER idid, info(15), ipar, iwork(34), n, liw, lrw
    REAL abserr, r, R1MACH, relerr, reltol, rpar, rwork(61), t, tout, u(4)
    !***FIRST EXECUTABLE STATEMENT  QXRKF
    IF ( Kprint>=2 ) WRITE (Lun,99001)
    !
    ! FORMATs.
    !
    99001 FORMAT ('1'/' ------------  DERKF QUICK CHECK OUTPUT',' ------------')
    !
    !     Initialize problem.
    !
    n = 4
    lrw = 61
    liw = 34
    t = 0.0E0
    tout = 8.0E0*ATAN(1.0E0)
    u(1) = 1.0E0
    u(2) = 0.0E0
    u(3) = 0.0E0
    u(4) = 1.0E0
    Ipass = 1
    reltol = SQRT(R1MACH(4))
    relerr = 0.1E0*reltol
    abserr = relerr**1.5E0
    info(1) = 0
    info(2) = 0
    info(3) = 1
    info(4) = 0
    IF ( Kprint>2 ) WRITE (Lun,99002) relerr, abserr, t, (1.0E0)
    99002 FORMAT (/' RELERR = ',E16.8,'   ABSERR =',E16.8/12X,'T',19X,'R'/2E20.8)
    DO
      !
      CALL DERKF(FDEQC,n,t,u,tout,info,relerr,abserr,idid,rwork,lrw,iwork,liw,&
        rpar,ipar)
      r = SQRT(u(1)*u(1)+u(2)*u(2))
      IF ( ABS(r-1.0E0)>reltol ) Ipass = 0
      IF ( Kprint>2 ) WRITE (Lun,99003) t, r
      99003 FORMAT (2E20.8)
      info(1) = 1
      IF ( idid/=1 ) THEN
        !
        !     Finish up.
        !
        IF ( idid<1 ) Ipass = 0
        IF ( Kprint>1.AND.idid<1 ) WRITE (Lun,99004) idid
        99004 FORMAT (1X,'ERROR RETURN FROM DERKF.  IDID = ',I3)
        IF ( Kprint>1.AND.Ipass==1 ) WRITE (Lun,99005)
        99005 FORMAT (/' ------------  DERKF PASSED TESTS  ------------')
        IF ( Kprint>=1.AND.Ipass==0 ) WRITE (Lun,99006)
        99006 FORMAT (/' ************  DERKF FAILED TESTS  ************')
        RETURN
      ENDIF
    ENDDO
  END SUBROUTINE QXRKF
  !DECK QXBVSP
  SUBROUTINE QXBVSP(Lun,Kprint,Ipass)
    IMPLICIT NONE
    REAL a, abser, ae, alpha, b, beta, re, reler, sve, TERm, tol, &
      work, xpts, XSAve, y, yans
    INTEGER i, iflag, igofx, Ipass, ipss, j, kont, kount, Kprint, l, &
      Lun, ncomp, ndiw, ndw, neqivp, nfc, nic, nrowa, nrowb, &
      nrowy
    INTEGER numort, nxpts
    !***BEGIN PROLOGUE  QXBVSP
    !***PURPOSE  Quick check for BVSUP.
    !***LIBRARY   SLATEC
    !***TYPE      SINGLE PRECISION (QXBVSP-S, QXDBVS-D)
    !***AUTHOR  (UNKNOWN)
    !***ROUTINES CALLED  BVSUP, PASS
    !***COMMON BLOCKS    SAVEX
    !***REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901014  Made editorial changes and added correct result to
    !           output.  (RWC)
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    !***END PROLOGUE  QXBVSP
    INTEGER itmp(9), iwork(100)
    DIMENSION y(4,15), xpts(15), a(2,4), alpha(2), b(2,4), beta(2), &
      yans(2,15), work(1000)
    CHARACTER(4) :: msg
    COMMON /SAVEX / XSAve, TERm
    DATA yans(1,1), yans(2,1), yans(1,2), yans(2,2), yans(1,3), yans(2,3)&
      , yans(1,4), yans(2,4), yans(1,5), yans(2,5), yans(1,6), &
      yans(2,6), yans(1,7), yans(2,7), yans(1,8), yans(2,8), yans(1,9)&
      , yans(2,9), yans(1,10), yans(2,10), yans(1,11), yans(2,11), &
      yans(1,12), yans(2,12), yans(1,13), yans(2,13), yans(1,14), &
      yans(2,14), yans(1,15), yans(2,15)/5.000000000E+00, &
      -6.888880126E-01, 8.609248635E+00, -1.083092311E+00, &
      1.674923836E+01, -2.072210073E+00, 3.351098494E+01, &
      -4.479263780E+00, 6.601103894E+01, -8.909222513E+00, &
      8.579580988E+01, -1.098742758E+01, 1.106536877E+02, &
      -1.402469444E+01, 1.421228220E+02, -1.742236546E+01, &
      1.803383474E+02, -2.086465851E+01, 2.017054332E+02, &
      -1.990879843E+01, 2.051622475E+02, -1.324886978E+01, &
      2.059197452E+02, 1.051529813E+01, 1.972191446E+02, &
      9.320592785E+01, 1.556894846E+02, 3.801682434E+02, &
      1.818989404E-12, 1.379853993E+03/
    DATA xpts(1), xpts(2), xpts(3), xpts(4), xpts(5), xpts(6), xpts(7), &
      xpts(8), xpts(9), xpts(10), xpts(11), xpts(12), xpts(13), &
      xpts(14), xpts(15)/60., 55., 50., 45., 40., 38., 36., 34., &
      32., 31., 30.8, 30.6, 30.4, 30.2, 30./
    !***FIRST EXECUTABLE STATEMENT  QXBVSP
    IF ( Kprint>=2 ) THEN
      WRITE (Lun,99001)
      !
      99001 FORMAT ('1')
      WRITE (Lun,99002)
      99002 FORMAT (/' BVSUP QUICK CHECK')
    ENDIF
    !
    !-----INITIALIZE VARIABLES FOR TEST PROBLEM.
    !
    DO i = 1, 9
      itmp(i) = 0
    ENDDO
    !
    tol = 1.0E-03
    XSAve = 0.
    nrowy = 4
    ncomp = 2
    nxpts = 15
    a(1,1) = 1.0
    a(1,2) = 0.0
    nrowa = 2
    alpha(1) = 5.0
    nic = 1
    b(1,1) = 1.0
    b(1,2) = 0.0
    nrowb = 2
    beta(1) = 0.0
    nfc = 1
    igofx = 1
    re = 1.0E-05
    ae = 1.0E-05
    ndw = 1000
    ndiw = 100
    neqivp = 0
    Ipass = 1
    !
    DO i = 1, 15
      iwork(i) = 0
    ENDDO
    !
    CALL BVSUP(y,nrowy,ncomp,xpts,nxpts,a,nrowa,alpha,nic,b,nrowb,beta,nfc,&
      igofx,re,ae,iflag,work,ndw,iwork,ndiw,neqivp)
    !
    !-----IF IFLAG = 0, WE HAVE A SUCCESSFUL SOLUTION; OTHERWISE, SKIP
    !     THE ARGUMENT CHECKING AND GO TO THE END.
    !
    IF ( iflag/=0 ) THEN
      Ipass = 0
      IF ( Kprint>1 ) WRITE (Lun,99003) iflag
      99003 FORMAT (10X,'IFLAG =',I2)
      GOTO 300
    ENDIF
    !
    !-----CHECK THE ACCURACY OF THE SOLUTION.
    !
    numort = iwork(1)
    DO j = 1, nxpts
      DO l = 1, 2
        abser = ABS(yans(l,j)-y(l,j))
        reler = abser/ABS(yans(l,j))
        IF ( reler>tol.AND.abser>tol ) Ipass = 0
      ENDDO
    ENDDO
    !
    !-----CHECK FOR SUPPRESSION OF PRINTING.
    !
    IF ( Kprint==0.OR.(Kprint==1.AND.Ipass==1) ) GOTO 400
    !
    IF ( Kprint/=1.OR.Ipass/=0 ) THEN
      IF ( Kprint>=3.OR.Ipass==0 ) THEN
        WRITE (Lun,99004)
        99004 FORMAT (/' ACCURACY TEST')
        WRITE (Lun,99005) numort
        99005 FORMAT (/' NUMBER OF ORTHONORMALIZATIONS =',I3)
        WRITE (Lun,99006) (work(j),j=1,numort)
        99006 FORMAT (/' ORTHONORMALIZATION POINTS ARE'/(1X,4F10.2))
        WRITE (Lun,99007)
        99007 FORMAT (//20X,'CALCULATION',30X,'TRUE SOLUTION'/2X,'X',14X,'Y',17X,&
          'Y-PRIME',15X,'Y',17X,'Y-PRIME'/)
        DO j = 1, nxpts
          msg = 'PASS'
          abser = ABS(yans(1,j)-y(1,j))
          reler = abser/ABS(yans(1,j))
          IF ( reler>tol.AND.abser>tol ) msg = 'FAIL'
          abser = ABS(yans(2,j)-y(2,j))
          reler = abser/ABS(yans(2,j))
          IF ( reler>tol.AND.abser>tol ) msg = 'FAIL'
          WRITE (Lun,99008) xpts(j), y(1,j), y(2,j), yans(1,j), yans(2,j), msg
          99008 FORMAT (F5.1,4E20.7,5X,A)
        ENDDO
      ENDIF
    ENDIF
    !
    !-----SEND MESSAGE INDICATING PASSAGE OR FAILURE OF TESTS.
    !
    CALL PASS(Lun,1,Ipass)
    !
    !-----ERROR MESSAGE TESTS.
    !
    IF ( Kprint==1 ) GOTO 400
    kont = 1
    WRITE (Lun,99009)
    99009 FORMAT (/' (7) TESTS OF IFLAG VALUES')
    !
    !-----NROWY LESS THAN NCOMP
    !
    kount = 1
    nrowy = 1
    100 CONTINUE
    DO
      DO i = 1, 15
        iwork(i) = 0
      ENDDO
      CALL BVSUP(y,nrowy,ncomp,xpts,nxpts,a,nrowa,alpha,nic,b,nrowb,beta,nfc,&
        igofx,re,ae,iflag,work,ndw,iwork,ndiw,neqivp)
      SELECT CASE (kount)
        CASE (2)
          !
          WRITE (Lun,99013) iflag
          IF ( iflag==-2 ) itmp(kont) = 1
          kont = kont + 1
          !
          !-----RE OR AE NEGATIVE
          !
          kount = 3
          igofx = 1
          re = -1.
          ae = -2.
        CASE (3)
          !
          WRITE (Lun,99013) iflag
          IF ( iflag==-2 ) itmp(kont) = 1
          kont = kont + 1
          !
          !-----NROWA LESS THAN NIC
          !
          kount = 4
          re = 1.0E-05
          ae = 1.0E-05
          nrowa = 0
          EXIT
        CASE (4)
          EXIT
        CASE (5)
          GOTO 200
        CASE (6)
          !
          WRITE (Lun,99010) iflag
          99010 FORMAT (/' IFLAG SHOULD BE -1, IFLAG =',I3)
          IF ( iflag==-1 ) itmp(kont) = 1
          kont = kont + 1
          !-----INCORRECT ORDERING OF XPTS
          kount = 7
          ndiw = 100
          sve = xpts(1)
          xpts(1) = xpts(4)
          xpts(4) = sve
        CASE (7)
          !
          WRITE (Lun,99013) iflag
          IF ( iflag==-2 ) itmp(kont) = 1
          GOTO 300
        CASE DEFAULT
          !
          WRITE (Lun,99013) iflag
          IF ( iflag==-2 ) itmp(kont) = 1
          kont = kont + 1
          !
          !-----IGOFX NOT EQUAL TO 0 OR 1
          !
          kount = 2
          nrowy = 2
          igofx = 3
      END SELECT
    ENDDO
    !
    WRITE (Lun,99013) iflag
    IF ( iflag==-2 ) itmp(kont) = 1
    kont = kont + 1
    !-----NROWB LESS THAN NFC
    kount = 5
    nrowa = 2
    nrowb = 0
    !
    200  WRITE (Lun,99013) iflag
    IF ( iflag==-2 ) itmp(kont) = 1
    kont = kont + 1
    !-----STORAGE ALLOCATION IS INSUFFICIENT
    kount = 6
    nrowb = 2
    ndiw = 17
    GOTO 100
    !
    !-----SEE IF IFLAG TESTS PASSED
    !
    300  ipss = 1
    DO i = 1, kont
      ipss = ipss*itmp(i)
    ENDDO
    !
    CALL PASS(Lun,2,ipss)
    !
    !     SEE IF ALL TESTS PASSED.
    !
    Ipass = Ipass*ipss
    !
    400 CONTINUE
    IF ( Ipass==1.AND.Kprint>1 ) WRITE (Lun,99011)
    99011 FORMAT (/' ****************BVSUP PASSED ALL TESTS***************')
    IF ( Ipass==0.AND.Kprint/=0 ) WRITE (Lun,99012)
    99012 FORMAT (/' ****************BVSUP FAILED SOME TESTS**************')
    RETURN
    99013 FORMAT (/' IFLAG SHOULD BE -2, IFLAG =',I3)
  END SUBROUTINE QXBVSP
END MODULE TEST43_MOD
!DECK TEST43
PROGRAM TEST43
  USE TEST43_MOD
  IMPLICIT NONE
  !***BEGIN PROLOGUE  TEST43
  !***PURPOSE  Driver for testing SLATEC subprograms
  !***LIBRARY   SLATEC
  !***CATEGORY  I1
  !***TYPE      SINGLE PRECISION (TEST43-S, TEST44-D)
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
  !        DEABM    DEBDF    DERKF    BVSUP
  !
  !***REFERENCES  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***ROUTINES CALLED  I1MACH, QXABM, QXBDF, QXBVSP, QXRKF, XERMAX, XSETF,
  !                    XSETUN
  !***REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  !***END PROLOGUE  TEST43
  INTEGER I1MACH
  INTEGER ipass, kprint, lin, lun, nfail
  !***FIRST EXECUTABLE STATEMENT  TEST43
  lun = I1MACH(2)
  lin = I1MACH(1)
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  CALL XERMAX(1000)
  CALL XSETUN(lun)
  IF ( kprint<=1 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  ENDIF
  !
  !     Test DEABM
  !
  CALL QXABM(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test DEBDF
  !
  CALL QXBDF(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test DERKF
  !
  CALL QXRKF(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test BVSUP
  !
  CALL QXBVSP(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST43 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST43 *************')
  ENDIF
  STOP
END PROGRAM TEST43
