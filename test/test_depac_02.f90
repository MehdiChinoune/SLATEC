MODULE TEST44_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  !** DJAC
  SUBROUTINE DJAC(T,U,Pd,Nrowpd)
    !> Evaluate Jacobian for DDEBDF quick check.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (JAC-S, DJAC-D)
    !***
    ! **Author:**  Chow, Jeff (LANL)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   810801  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   900415  Minor clean-up of prologue and code and name changed from
    !           DDJAC to DJAC.  (WRB)

    INTEGER :: Nrowpd
    REAL(DP) :: T, U(:), Pd(:,:)
    REAL(DP) :: r, r5, rsq, u1sq, u2sq, u1u2
    !* FIRST EXECUTABLE STATEMENT  DJAC
    u1sq = U(1)*U(1)
    u2sq = U(2)*U(2)
    u1u2 = U(1)*U(2)
    rsq = u1sq + u2sq
    r = SQRT(rsq)
    r5 = rsq*rsq*r
    Pd(3,1) = (3.D0*u1sq-rsq)/r5
    Pd(4,1) = 3.D0*u1u2/r5
    Pd(3,2) = Pd(4,1)
    Pd(4,2) = (3.D0*u2sq-rsq)/r5
    Pd(1,3) = 1.D0
    Pd(2,4) = 1.D0
  END SUBROUTINE DJAC
  !** DFDEQC
  SUBROUTINE DFDEQC(T,U,Uprime)
    !> Derivative evaluator for DDEPAC quick checks.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (FDEQC-S, DFDEQC-D)
    !***
    ! **Author:**  Chow, Jeff, (LANL)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   810801  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   900415  Name changed from DDF to DFDEQC.  (WRB)

    !
    !     Declare arguments.
    !
    REAL(DP) :: T, U(:), Uprime(:)
    !
    !     Declare local variables.
    !
    REAL(DP) :: r, rsq, r3
    !* FIRST EXECUTABLE STATEMENT  DFDEQC
    rsq = U(1)*U(1) + U(2)*U(2)
    r = SQRT(rsq)
    r3 = rsq*r
    Uprime(1) = U(3)
    Uprime(2) = U(4)
    Uprime(3) = -(U(1)/r3)
    Uprime(4) = -(U(2)/r3)
  END SUBROUTINE DFDEQC
  !** QXDABM
  SUBROUTINE QXDABM(Lun,Kprint,Ipass)
    !> Test the DEPAC routine DDEABM.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (QXABM-S, QXDABM-D)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Chow, Jeff, (LANL)
    !***
    ! **Description:**
    !
    !- Usage:
    !
    !        INTEGER  LUN, KPRINT, IPASS
    !
    !        CALL QXDABM (LUN, KPRINT, IPASS)
    !
    !- Arguments:
    !
    !     LUN   :IN  is the unit number to which output is to be written.
    !
    !     KPRINT:IN  controls the amount of output, as specified in the
    !                SLATEC Guidelines.
    !
    !     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
    !                IPASS=0 indicates one or more tests failed.
    !
    !- Description:
    !
    !   DDEABM is tested by solving the equations of motion of a body
    !   moving in a plane about a spherical earth, namely
    !           (D/DT)(D/DT)X = -G*X/R**3
    !           (D/DT)(D/DT)Y = -G*Y/R**3
    !   where G = 1, R = SQRT(X**2 + Y**2) and
    !           X(0) = 1
    !           (D/DT)X(0) = 0
    !           Y(0) = 0
    !           (D/DT)Y(0) = 1.
    !
    !***
    ! **Routines called:**  D1MACH, DDEABM, DFDEQC

    !* REVISION HISTORY  (YYMMDD)
    !   810801  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   900415  Code extensively revised.  (WRB)
    USE slatec, ONLY : D1MACH, DDEABM
    !
    !     Declare arguments.
    !
    INTEGER :: Lun, Kprint, Ipass
    !
    !     Declare local variables.
    !
    INTEGER :: idid, info(15), iwork(51), n, liw, lrw, nstep
    REAL(DP) :: abserr(1), r, relerr(1), reltol, rwork(214), t, tout, u(4)
    !* FIRST EXECUTABLE STATEMENT  QXDABM
    IF( Kprint>=2 ) WRITE (Lun,99001)
    !
    ! FORMATs.
    !
    99001 FORMAT ('1'/' ------------  DDEABM QUICK CHECK OUTPUT',' ------------')
    !
    !     Initialize problem.
    !
    n = 4
    lrw = 214
    liw = 51
    t = 0.0D0
    tout = 8.0D0*ATAN(1.0D0)
    u(1) = 1.0D0
    u(2) = 0.0D0
    u(3) = 0.0D0
    u(4) = 1.0D0
    Ipass = 1
    nstep = 0
    reltol = SQRT(D1MACH(4))
    relerr = 0.1D0*reltol
    abserr = relerr**1.5D0
    info(1) = 0
    info(2) = 0
    info(3) = 1
    info(4) = 0
    IF( Kprint>2 ) WRITE (Lun,99002) relerr, abserr, t, (1.0D0)
    99002 FORMAT (/' RELERR = ',D16.8,'   ABSERR =',D16.8/12X,'T',19X,'R'/2D20.8)
    DO
      !
      CALL DDEABM(DFDEQC,n,t,u,tout,info,relerr,abserr,idid,rwork,lrw,iwork,liw)
      r = SQRT(u(1)*u(1)+u(2)*u(2))
      IF( ABS(r-1.0D0)>reltol ) Ipass = 0
      IF( Kprint>2 ) WRITE (Lun,99003) t, r
      99003 FORMAT (2D20.8)
      info(1) = 1
      IF( idid/=1 ) THEN
        !
        !     For the double precision version, we allow the integrator to take
        !     up to 2000 steps before we declare failure.
        !
        IF( idid==-1 ) THEN
          nstep = nstep + 500
          IF( nstep<2000 ) CYCLE
        END IF
        !
        !     Finish up.
        !
        IF( idid<1 ) Ipass = 0
        IF( Kprint>1 .AND. idid<1 ) WRITE (Lun,99004) idid
        99004 FORMAT (1X,'ERROR RETURN FROM DDEABM.  IDID = ',I3)
        IF( Kprint>1 .AND. Ipass==1 ) WRITE (Lun,99005)
        99005 FORMAT (/' ------------  DDEABM PASSED TESTS  ------------')
        IF( Kprint>=1 .AND. Ipass==0 ) WRITE (Lun,99006)
        99006 FORMAT (/' ************  DDEABM FAILED TESTS  ************')
        RETURN
      END IF
    END DO
  END SUBROUTINE QXDABM
  !** QXDBDF
  SUBROUTINE QXDBDF(Lun,Kprint,Ipass)
    !> Test the DEPAC routine DDEBDF.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (QXBDF-S, QXDBDF-D)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Chow, Jeff, (LANL)
    !***
    ! **Description:**
    !
    !- Usage:
    !
    !        INTEGER  LUN, KPRINT, IPASS
    !
    !        CALL QXDBDF (LUN, KPRINT, IPASS)
    !
    !- Arguments:
    !
    !     LUN   :IN  is the unit number to which output is to be written.
    !
    !     KPRINT:IN  controls the amount of output, as specified in the
    !                SLATEC Guidelines.
    !
    !     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
    !                IPASS=0 indicates one or more tests failed.
    !
    !- Description:
    !
    !   DDEBDF is tested by solving the equations of motion of a body
    !   moving in a plane about a spherical earth, namely
    !           (D/DT)(D/DT)X = -G*X/R**3
    !           (D/DT)(D/DT)Y = -G*Y/R**3
    !   where G = 1, R = SQRT(X**2 + Y**2) and
    !           X(0) = 1
    !           (D/DT)X(0) = 0
    !           Y(0) = 0
    !           (D/DT)Y(0) = 1.
    !
    !***
    ! **Routines called:**  D1MACH, DDEBDF, DFDEQC, DJAC

    !* REVISION HISTORY  (YYMMDD)
    !   810801  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   900415  Code extensively revised.  (WRB)
    USE slatec, ONLY : D1MACH, DDEBDF
    !
    !     Declare arguments.
    !
    INTEGER :: Lun, Kprint, Ipass
    !
    !     Declare local variables.
    !
    INTEGER :: idid, info(15), iwork(60), n, liw, lrw, nstep
    REAL(DP) :: abserr(1), r, reltol, relerr(1), rwork(306), t, tout, u(4)
    !* FIRST EXECUTABLE STATEMENT  QXDBDF
    IF( Kprint>=2 ) WRITE (Lun,99001)
    !
    ! FORMATs.
    !
    99001 FORMAT ('1'/' ------------  DDEBDF QUICK CHECK OUTPUT',' ------------')
    !
    !     Initialize problem.
    !
    n = 4
    lrw = 306
    liw = 60
    t = 0.0D0
    tout = 8.0D0*ATAN(1.0D0)
    u(1) = 1.0D0
    u(2) = 0.0D0
    u(3) = 0.0D0
    u(4) = 1.0D0
    Ipass = 1
    nstep = 0
    reltol = MAX(SQRT(D1MACH(4)),1.D-9)
    relerr = MAX(0.0001D0*reltol,1.D-12)
    abserr = relerr**1.5D0
    info(1) = 0
    info(2) = 0
    info(3) = 1
    info(4) = 0
    info(5) = 1
    info(6) = 0
    IF( Kprint>2 ) WRITE (Lun,99002) relerr, abserr, t, (1.0D0)
    99002 FORMAT (/' RELERR = ',D16.8,'   ABSERR =',D16.8/12X,'T',19X,'R'/2D20.8)
    DO
      !
      CALL DDEBDF(DFDEQC,n,t,u,tout,info,relerr,abserr,idid,rwork,lrw,iwork,&
        liw,DJAC)
      r = SQRT(u(1)*u(1)+u(2)*u(2))
      IF( ABS(r-1.0D0)>reltol ) Ipass = 0
      IF( Kprint>2 ) WRITE (Lun,99003) t, r
      99003 FORMAT (2D20.8)
      info(1) = 1
      IF( idid/=1 ) THEN
        !
        !     For the double precision version, we allow the integrator to take
        !     up to 2000 steps before we declare failure.
        !
        IF( idid==-1 ) THEN
          nstep = nstep + 500
          IF( nstep<2000 ) CYCLE
        END IF
        !
        !     Finish up.
        !
        IF( idid<1 ) Ipass = 0
        IF( Kprint>1 .AND. idid<1 ) WRITE (Lun,99004) idid
        99004 FORMAT (1X,'ERROR RETURN FROM DDEBDF.  IDID = ',I3)
        IF( Kprint>1 .AND. Ipass==1 ) WRITE (Lun,99005)
        99005 FORMAT (/' ------------  DDEBDF PASSED TESTS  ------------')
        IF( Kprint>=1 .AND. Ipass==0 ) WRITE (Lun,99006)
        99006 FORMAT (/' ************  DDEBDF FAILED TESTS  ************')
        RETURN
      END IF
    END DO
  END SUBROUTINE QXDBDF
  !** QXDBVS
  SUBROUTINE QXDBVS(Lun,Kprint,Ipass)
    !> Quick check for DBVSUP.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (QXBVSP-S, QXDBVS-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  DBVSUP, PASS
    !***
    ! COMMON BLOCKS    DSAVEX

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901014  Made editorial changes and added correct result to
    !           output. (RWC)
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    USE DSAVEX, ONLY : xsave_com
    USE slatec, ONLY : DBVSUP
    USE common_mod, ONLY : PASS
    INTEGER :: i, iflag, igofx, Ipass, ipss, j, kont, kount, Kprint, l, &
      Lun, ncomp, ndiw, ndw, neqivp, nfc, nic, nrowa, nrowb, nrowy
    INTEGER :: numort, nxpts
    INTEGER :: itmp(9), iwork(100)
    REAL(DP) :: work(1000), ae, re, sve, tol
    REAL(DP) :: y(4,15), a(2,4), alpha(2), b(2,4), beta(2), reler, abser
    CHARACTER(4) :: msg
    REAL(DP), PARAMETER :: yans(2,15) = RESHAPE( [ &
      5.000000000D+00, -6.888880126D-01,   8.609248635D+00, -1.083092311D+00, &
      1.674923836D+01, -2.072210073D+00,   3.351098494D+01, -4.479263780D+00, &
      6.601103894D+01, -8.909222513D+00,   8.579580988D+01, -1.098742758D+01, &
      1.106536877D+02, -1.402469444D+01,   1.421228220D+02, -1.742236546D+01, &
      1.803383474D+02, -2.086465851D+01,   2.017054332D+02, -1.990879843D+01, &
      2.051622475D+02, -1.324886978D+01,   2.059197452D+02, 1.051529813D+01, &
      1.972191446D+02, 9.320592785D+01,    1.556894846D+02, 3.801682434D+02, &
      1.818989404D-12, 1.379853993D+03 ], [2,15] )
    REAL(DP) :: xpts(15) = [ 60.0D+00, 55.0D+00, 50.0D+00, 45.0D+00, &
      40.0D+00, 38.0D+00, 36.0D+00, 34.0D+00, 32.0D+00, 31.0D+00, &
      30.8D+00, 30.6D+00, 30.4D+00, 30.2D+00, 30.0D+00 ]
    !* FIRST EXECUTABLE STATEMENT  QXDBVS
    IF( Kprint>=2 ) THEN
      WRITE (Lun,99001)
      !
      99001 FORMAT ('1')
      WRITE (Lun,99002)
      99002 FORMAT (/' DBVSUP QUICK CHECK')
    END IF
    !
    !-----INITIALIZE VARIABLES FOR TEST PROBLEM.
    !
    DO i = 1, 9
      itmp(i) = 0
    END DO
    !
    tol = 1.0D-03
    xsave_com = 0.0D+00
    nrowy = 4
    ncomp = 2
    nxpts = 15
    a(1,1) = 1.0D+00
    a(1,2) = 0.0D+00
    nrowa = 2
    alpha(1) = 5.0D+00
    nic = 1
    b(1,1) = 1.0D+00
    b(1,2) = 0.0D+00
    nrowb = 2
    beta(1) = 0.0D+00
    nfc = 1
    igofx = 1
    re = 1.0D-05
    ae = 1.0D-05
    ndw = 1000
    ndiw = 100
    neqivp = 0
    Ipass = 1
    !
    DO i = 1, 15
      iwork(i) = 0
    END DO
    !
    CALL DBVSUP(y,nrowy,ncomp,xpts,nxpts,a,nrowa,alpha,nic,b,nrowb,beta,nfc,&
      igofx,re,ae,iflag,work,ndw,iwork,ndiw,neqivp)
    !
    !-----IF IFLAG = 0, WE HAVE A SUCCESSFUL SOLUTION; OTHERWISE, SKIP
    !     THE ARGUMENT CHECKING AND GO TO THE END.
    !
    IF( iflag/=0 ) THEN
      Ipass = 0
      IF( Kprint>1 ) WRITE (Lun,99003) iflag
      99003 FORMAT (10X,'IFLAG =',I2)
      GOTO 300
    END IF
    !
    !-----CHECK THE ACCURACY OF THE SOLUTION.
    !
    numort = iwork(1)
    DO j = 1, nxpts
      DO l = 1, 2
        abser = ABS(yans(l,j)-y(l,j))
        reler = abser/ABS(yans(l,j))
        IF( reler>tol .AND. abser>tol ) Ipass = 0
      END DO
    END DO
    !
    !-----CHECK FOR SUPPRESSION OF PRINTING.
    !
    IF( Kprint==0 .OR. (Kprint==1 .AND. Ipass==1) ) GOTO 400
    !
    IF( Kprint/=1 .OR. Ipass/=0 ) THEN
      IF( Kprint>=3 .OR. Ipass==0 ) THEN
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
          IF( reler>tol .AND. abser>tol ) msg = 'FAIL'
          abser = ABS(yans(2,j)-y(2,j))
          reler = abser/ABS(yans(2,j))
          IF( reler>tol .AND. abser>tol ) msg = 'FAIL'
          WRITE (Lun,99008) xpts(j), y(1,j), y(2,j), yans(1,j), yans(2,j), msg
          99008 FORMAT (F5.1,4E20.7,5X,A)
        END DO
      END IF
    END IF
    !
    !-----SEND MESSAGE INDICATING PASSAGE OR FAILURE OF TESTS.
    !
    CALL PASS(Lun,1,Ipass)
    !
    !-----ERROR MESSAGE TESTS.
    !
    IF( Kprint==1 ) GOTO 400
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
      END DO
      CALL DBVSUP(y,nrowy,ncomp,xpts,nxpts,a,nrowa,alpha,nic,b,nrowb,beta,nfc,&
        igofx,re,ae,iflag,work,ndw,iwork,ndiw,neqivp)
      SELECT CASE (kount)
        CASE (2)
          !
          WRITE (Lun,99013) iflag
          IF( iflag==-2 ) itmp(kont) = 1
          kont = kont + 1
          !
          !-----RE OR AE NEGATIVE
          !
          kount = 3
          igofx = 1
          re = -1.0D+00
          ae = -2.0D+00
        CASE (3)
          !
          WRITE (Lun,99013) iflag
          IF( iflag==-2 ) itmp(kont) = 1
          kont = kont + 1
          !
          !-----NROWA LESS THAN NIC
          !
          kount = 4
          re = 1.0D-05
          ae = 1.0D-05
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
          IF( iflag==-1 ) itmp(kont) = 1
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
          IF( iflag==-2 ) itmp(kont) = 1
          GOTO 300
        CASE DEFAULT
          !
          WRITE (Lun,99013) iflag
          IF( iflag==-2 ) itmp(kont) = 1
          kont = kont + 1
          !
          !-----IGOFX NOT EQUAL TO 0 OR 1
          !
          kount = 2
          nrowy = 2
          igofx = 3
      END SELECT
    END DO
    !
    WRITE (Lun,99013) iflag
    IF( iflag==-2 ) itmp(kont) = 1
    kont = kont + 1
    !-----NROWB LESS THAN NFC
    kount = 5
    nrowa = 2
    nrowb = 0
    !
    200  WRITE (Lun,99013) iflag
    IF( iflag==-2 ) itmp(kont) = 1
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
    END DO
    !
    CALL PASS(Lun,2,ipss)
    !
    !-----SEE IF ALL TESTS PASSED.
    !
    Ipass = Ipass*ipss
    !
    400 CONTINUE
    IF( Ipass==1 .AND. Kprint>1 ) WRITE (Lun,99011)
    99011 FORMAT (/' ***************DBVSUP PASSED ALL TESTS***************')
    IF( Ipass==0 .AND. Kprint/=0 ) WRITE (Lun,99012)
    99012 FORMAT (/' ***************DBVSUP FAILED SOME TESTS**************')
    RETURN
    99013 FORMAT (/' IFLAG SHOULD BE -2, IFLAG =',I3)
  END SUBROUTINE QXDBVS
  !** QXDRKF
  SUBROUTINE QXDRKF(Lun,Kprint,Ipass)
    !> Test the DEPAC routine DDERKF.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (QXRKF-S, QXDRKF-D)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Chow, Jeff, (LANL)
    !***
    ! **Description:**
    !
    !- Usage:
    !
    !        INTEGER  LUN, KPRINT, IPASS
    !
    !        CALL QXDRKF (LUN, KPRINT, IPASS)
    !
    !- Arguments:
    !
    !     LUN   :IN  is the unit number to which output is to be written.
    !
    !     KPRINT:IN  controls the amount of output, as specified in the
    !                SLATEC Guidelines.
    !
    !     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
    !                IPASS=0 indicates one or more tests failed.
    !
    !- Description:
    !
    !   DDERKF is tested by solving the equations of motion of a body
    !   moving in a plane about a spherical earth, namely
    !           (D/DT)(D/DT)X = -G*X/R**3
    !           (D/DT)(D/DT)Y = -G*Y/R**3
    !   where G = 1, R = SQRT(X**2 + Y**2) and
    !           X(0) = 1
    !           (D/DT)X(0) = 0
    !           Y(0) = 0
    !           (D/DT)Y(0) = 1.
    !
    !***
    ! **Routines called:**  D1MACH, DDERKF, DFDEQC

    !* REVISION HISTORY  (YYMMDD)
    !   810801  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   900415  Code extensively revised.  (WRB)
    USE slatec, ONLY : D1MACH, DDERKF
    !
    !     Declare arguments.
    !
    INTEGER :: Lun, Kprint, Ipass
    !
    !     Declare local variables.
    !
    INTEGER :: idid, info(15), iwork(34), n, liw, lrw, nstep
    REAL(DP) :: abserr(1), r, relerr(1), reltol, rwork(61), t, tout, u(4)
    !* FIRST EXECUTABLE STATEMENT  QXDRKF
    IF( Kprint>=2 ) WRITE (Lun,99001)
    !
    ! FORMATs.
    !
    99001 FORMAT ('1'/' ------------  DDERKF QUICK CHECK OUTPUT',' ------------')
    !
    !     Initialize problem.
    !
    n = 4
    lrw = 61
    liw = 34
    t = 0.0D0
    tout = 8.0D0*ATAN(1.0D0)
    u(1) = 1.0D0
    u(2) = 0.0D0
    u(3) = 0.0D0
    u(4) = 1.0D0
    Ipass = 1
    nstep = 0
    reltol = MAX(SQRT(D1MACH(4)),1.D-10)
    relerr = MAX(.1D0*reltol,1.D-12)
    abserr = relerr**1.5D0
    info(1) = 0
    info(2) = 0
    info(3) = 1
    info(4) = 0
    IF( Kprint>2 ) WRITE (Lun,99002) relerr, abserr, t, (1.0D0)
    99002 FORMAT (/' RELERR = ',D16.8,'   ABSERR =',D16.8/12X,'T',19X,'R'/2D20.8)
    DO
      !
      CALL DDERKF(DFDEQC,n,t,u,tout,info,relerr,abserr,idid,rwork,lrw,iwork,liw)
      r = SQRT(u(1)*u(1)+u(2)*u(2))
      IF( ABS(r-1.0D0)>reltol ) Ipass = 0
      IF( Kprint>2 ) WRITE (Lun,99003) t, r
      99003 FORMAT (2D20.8)
      info(1) = 1
      IF( idid/=1 ) THEN
        !
        !     For the double precision version, we allow the integrator to take
        !     up to 2000 steps before we declare failure.
        !
        IF( idid==-1 ) THEN
          nstep = nstep + 500
          IF( nstep<2000 ) CYCLE
        END IF
        !
        !     Finish up.
        !
        IF( idid<1 ) Ipass = 0
        IF( Kprint>1 .AND. idid<1 ) WRITE (Lun,99004) idid
        99004 FORMAT (1X,'ERROR RETURN FROM DDERKF.  IDID = ',I3)
        IF( Kprint>1 .AND. Ipass==1 ) WRITE (Lun,99005)
        99005 FORMAT (/' ------------  DDERKF PASSED TESTS  ------------')
        IF( Kprint>=1 .AND. Ipass==0 ) WRITE (Lun,99006)
        99006 FORMAT (/' ************  DDERKF FAILED TESTS  ************')
        RETURN
      END IF
    END DO
  END SUBROUTINE QXDRKF
END MODULE TEST44_MOD
!** TEST44
PROGRAM TEST44
  USE TEST44_MOD, ONLY : QXDABM, QXDBDF, QXDBVS, QXDRKF
  USE slatec, ONLY : I1MACH, XSETF, XSETUN, XERMAX
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !> Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  I1
  !***
  ! **Type:**      DOUBLE PRECISION (TEST43-S, TEST44-D)
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
  !        DDEABM   DDEBDF   DDERKF   DBVSUP
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  I1MACH, QXDABM, QXDBDF, QXDBVS, QXDRKF, XERMAX,
  !                    XSETF, XSETUN

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  INTEGER :: ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST44
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
  !     Test DDEABM
  !
  CALL QXDABM(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Test DDEBDF
  !
  CALL QXDBDF(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Test DDERKF
  !
  CALL QXDRKF(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Test DBVSUP
  !
  CALL QXDBVS(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST44 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST44 *************')
  END IF
  STOP
END PROGRAM TEST44
