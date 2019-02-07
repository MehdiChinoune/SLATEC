!*==QXBDF.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK QXBDF
SUBROUTINE QXBDF(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--QXBDF5
  !*** Start of declarations inserted by SPAG
  REAL FDEQC
  INTEGER JAC
  !*** End of declarations inserted by SPAG
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
  REAL abserr, r, R1MACH, relerr, reltol, rpar, rwork(306), t, &
    tout, u(4)
  EXTERNAL FDEQC, JAC
  !***FIRST EXECUTABLE STATEMENT  QXBDF
  IF ( Kprint>=2 ) WRITE (Lun,99001)
  !
  !     FORMATs.
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
    99003   FORMAT (2E20.8)
    info(1) = 1
    IF ( idid/=1 ) THEN
      !
      !     Finish up.
      !
      IF ( idid<1 ) Ipass = 0
      IF ( Kprint>1.AND.idid<1 ) WRITE (Lun,99004) idid
      99004     FORMAT (1X,'ERROR RETURN FROM DEBDF.  IDID = ',I3)
      IF ( Kprint>1.AND.Ipass==1 ) WRITE (Lun,99005)
      99005     FORMAT (/' ------------  DEBDF PASSED TESTS  ------------')
      IF ( Kprint>=1.AND.Ipass==0 ) WRITE (Lun,99006)
      99006     FORMAT (/' ************  DEBDF FAILED TESTS  ************')
      RETURN
    ENDIF
  ENDDO
END SUBROUTINE QXBDF
