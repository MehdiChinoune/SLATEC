!*==QXDABM.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK QXDABM
SUBROUTINE QXDABM(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--QXDABM5
  !*** Start of declarations inserted by SPAG
  REAL DFDEQC
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  QXDABM
  !***PURPOSE  Test the DEPAC routine DDEABM.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (QXABM-S, QXDABM-D)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Chow, Jeff, (LANL)
  !***DESCRIPTION
  !
  ! *Usage:
  !
  !        INTEGER  LUN, KPRINT, IPASS
  !
  !        CALL QXDABM (LUN, KPRINT, IPASS)
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
  !***ROUTINES CALLED  D1MACH, DDEABM, DFDEQC
  !***REVISION HISTORY  (YYMMDD)
  !   810801  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900415  Code extensively revised.  (WRB)
  !***END PROLOGUE  QXDABM
  !
  !     Declare arguments.
  !
  INTEGER Lun, Kprint, Ipass
  !
  !     Declare local variables.
  !
  INTEGER idid, info(15), ipar, iwork(51), n, liw, lrw, nstep
  REAL(8) :: abserr, D1MACH, r, relerr, reltol, rpar, rwork(214)&
    , t, tout, u(4)
  EXTERNAL DFDEQC
  !***FIRST EXECUTABLE STATEMENT  QXDABM
  IF ( Kprint>=2 ) WRITE (Lun,99001)
  !
  !     FORMATs.
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
  IF ( Kprint>2 ) WRITE (Lun,99002) relerr, abserr, t, (1.0D0)
  99002 FORMAT (/' RELERR = ',D16.8,'   ABSERR =',D16.8/12X,'T',19X,'R'/2D20.8)
  DO
    !
    CALL DDEABM(DFDEQC,n,t,u,tout,info,relerr,abserr,idid,rwork,lrw,iwork,&
      liw,rpar,ipar)
    r = SQRT(u(1)*u(1)+u(2)*u(2))
    IF ( ABS(r-1.0D0)>reltol ) Ipass = 0
    IF ( Kprint>2 ) WRITE (Lun,99003) t, r
    99003   FORMAT (2D20.8)
    info(1) = 1
    IF ( idid/=1 ) THEN
      !
      !     For the double precision version, we allow the integrator to take
      !     up to 2000 steps before we declare failure.
      !
      IF ( idid==-1 ) THEN
        nstep = nstep + 500
        IF ( nstep<2000 ) CYCLE
      ENDIF
      !
      !     Finish up.
      !
      IF ( idid<1 ) Ipass = 0
      IF ( Kprint>1.AND.idid<1 ) WRITE (Lun,99004) idid
      99004     FORMAT (1X,'ERROR RETURN FROM DDEABM.  IDID = ',I3)
      IF ( Kprint>1.AND.Ipass==1 ) WRITE (Lun,99005)
      99005     FORMAT (/' ------------  DDEABM PASSED TESTS  ------------')
      IF ( Kprint>=1.AND.Ipass==0 ) WRITE (Lun,99006)
      99006     FORMAT (/' ************  DDEABM FAILED TESTS  ************')
      RETURN
    ENDIF
  ENDDO
END SUBROUTINE QXDABM
