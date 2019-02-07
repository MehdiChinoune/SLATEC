!*==PCHQK2.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK PCHQK2
SUBROUTINE PCHQK2(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--PCHQK25
  !***BEGIN PROLOGUE  PCHQK2
  !***PURPOSE  Test the PCHIP integrators PCHIA and PCHID.
  !***LIBRARY   SLATEC (PCHIP)
  !***TYPE      SINGLE PRECISION (PCHQK2-S, DPCHQ2-D)
  !***KEYWORDS  PCHIP INTEGRATOR QUICK CHECK
  !***AUTHOR  Fritsch, F. N., (LLNL)
  !***DESCRIPTION
  !
  !              PCHIP QUICK CHECK NUMBER 2
  !
  !     TESTS THE INTEGRATORS:  PCHIA, PCHID.
  ! *Usage:
  !
  !        INTEGER  LUN, KPRINT, IPASS
  !
  !        CALL PCHQK2 (LUN, KPRINT, IPASS)
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
  !   This routine constructs data from a cubic, integrates it with PCHIA
  !   and compares the results with the correct answer.
  !   Since PCHIA calls PCHID, this tests both integrators.
  !
  !***ROUTINES CALLED  PCHIA, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   820601  DATE WRITTEN
  !   890306  Changed IPASS to the more accurate name IFAIL.  (FNF)
  !   890316  Added declarations as in DPCHQ2.  (FNF)
  !   890629  Appended E0 to real constants to reduce S.P./D.P.
  !           differences.
  !   890706  Cosmetic changes to prologue.  (WRB)
  !   891004  Cosmetic changes to prologue.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900314  Improved some output formats.  (FNF)
  !   900315  Revised prologue and improved some output formats.  (FNF)
  !   900316  Additional minor cosmetic changes.  (FNF)
  !   900321  Removed IFAIL from call sequence for SLATEC standards and
  !           made miscellaneous cosmetic changes.  (FNF)
  !   901130  Added 1P's to formats; changed to allow KPRINT.gt.3.  (FNF)
  !   910708  Minor modifications in use of KPRINT.  (WRB)
  !   930317  Improved output formats.  (FNF)
  !***END PROLOGUE  PCHQK2
  !
  !  Declare arguments.
  !
  INTEGER Lun, Kprint, Ipass
  !
  !  DECLARE VARIABLES.
  !
  INTEGER i, ierexp(17), ierr, ifail, n, npairs
  REAL a(17), b(17), calc, d(7), errmax, error, f(7), machep, one ,&
    three, thrqtr, tol, true, two, x(7)
  LOGICAL fail, skip
  !
  !  DECLARE EXTERNALS.
  !
  REAL PCHIA, R1MACH
  !
  !  INITIALIZE.
  !
  DATA thrqtr/0.75E0/, one/1.E0/, two/2.E0/, three/3.E0/
  DATA n/7/
  DATA x/ - 4.E0, -2.E0, -0.9E0, 0.E0, 0.9E0, 2.E0, 4.E0/
  DATA npairs/17/
  DATA a/ - 3.0E0, 3.0E0, -0.5E0, -0.5E0, -0.5E0, -4.0E0, -4.0E0 ,&
    3.0E0, -5.0E0, -5.0E0, -6.0E0, 6.0E0, -1.5E0, -1.5E0, -3.0E0 ,&
    3.0E0, 0.5E0/
  DATA b/3.0E0, -3.0E0, 1.0E0, 2.0E0, 5.0E0, -0.5E0, 4.0E0, 5.0E0 ,&
    -3.0E0, 5.0E0, -5.0E0, 5.0E0, -0.5E0, -1.0E0, -2.5E0, 3.5E0 ,&
    0.5E0/
  DATA ierexp/0, 0, 0, 0, 2, 0, 0, 2, 1, 3, 3, 3, 0, 0, 0 ,&
    0, 0/
  !
  !  SET PASS/FAIL TOLERANCE.
  !
  !***FIRST EXECUTABLE STATEMENT  PCHQK2
  machep = R1MACH(4)
  tol = 100.E0*machep
  !
  !  SET UP PCH FUNCTION DEFINITION.
  !
  DO i = 1, n
    f(i) = FCN(x(i))
    d(i) = DERIV(x(i))
  ENDDO
  !
  IF ( Kprint>=3 ) WRITE (Lun,99001)
  !
  !  FORMATS.
  !
  99001 FORMAT ('1'//10X,'TEST PCHIP INTEGRATORS')
  IF ( Kprint>=2 ) WRITE (Lun,99002)
  99002 FORMAT (//10X,'PCHQK2 RESULTS'/10X,'--------------')
  IF ( Kprint>=3 ) WRITE (Lun,99003) (x(i),f(i),d(i),i=1,n)
  99003 FORMAT (//5X,'DATA:'//11X,'X',9X,'F',9X,'D'/(5X,3F10.3))
  !
  !  LOOP OVER (A,B)-PAIRS.
  !
  IF ( Kprint>=3 ) WRITE (Lun,99004)
  99004 FORMAT (//5X,'TEST RESULTS:'//'    A     B    ERR     TRUE',16X,'CALC',&
    15X,'ERROR')
  !
  ifail = 0
  !
  skip = .FALSE.
  DO i = 1, npairs
    !               ---------------------------------------------
    calc = PCHIA(n,x,f,d,1,skip,a(i),b(i),ierr)
    !               ---------------------------------------------
    IF ( ierr>=0 ) THEN
      fail = ierr/=ierexp(i)
      true = ANTDER(b(i)) - ANTDER(a(i))
      error = calc - true
      IF ( Kprint>=3 ) THEN
        IF ( fail ) THEN
          WRITE (Lun,99005) a(i), b(i), ierr, true, calc, error ,&
            ierexp(i)
          99005         FORMAT (2F6.1,I5,1P,2E20.10,E15.5,'  (',I1,') *****')
        ELSE
          WRITE (Lun,99010) a(i), b(i), ierr, true, calc, error
        ENDIF
      ENDIF
      !
      error = ABS(error)/MAX(one,ABS(true))
      IF ( fail.OR.(error>tol) ) ifail = ifail + 1
      IF ( i==1 ) THEN
        errmax = error
      ELSE
        errmax = MAX(errmax,error)
      ENDIF
    ELSE
      IF ( Kprint>=3 ) WRITE (Lun,99010) a(i), b(i), ierr
      ifail = ifail + 1
    ENDIF
  ENDDO
  !
  !  PRINT SUMMARY.
  !
  IF ( Kprint>=2 ) THEN
    WRITE (Lun,99006) errmax, tol
    99006   FORMAT (/'  MAXIMUM RELATIVE ERROR IS:',1P,E15.5,',   TOLERANCE:',1P,&
      E15.5)
    IF ( ifail/=0 ) WRITE (Lun,99007) ifail
    99007   FORMAT (/' *** TROUBLE ***',I5,' INTEGRATION TESTS FAILED.')
  ENDIF
  !
  !  TERMINATE.
  !
  IF ( ifail==0 ) THEN
    Ipass = 1
    IF ( Kprint>=2 ) WRITE (Lun,99008)
    99008   FORMAT (/' ------------  PCHIP PASSED  ALL INTEGRATION TESTS',&
      ' ------------')
  ELSE
    Ipass = 0
    IF ( Kprint>=1 ) WRITE (Lun,99009)
    99009   FORMAT (/' ************  PCHIP FAILED SOME INTEGRATION TESTS',&
      ' ************')
  ENDIF
  !
  RETURN
  99010 FORMAT (2F6.1,I5,1P,2E20.10,E15.5)
  !------------- LAST LINE OF PCHQK2 FOLLOWS -----------------------------
CONTAINS

  !
  !  DEFINE TEST FUNCTIONS.
  !
  REAL FUNCTION FCN(ax)
    REAL, INTENT(IN) :: ax
    FCN = three*ax*ax*(ax-two)
  END FUNCTION FCN
  REAL FUNCTION DERIV(ax)
    REAL, INTENT(IN) :: ax
    DERIV = three*ax*(two*(ax-two)+ax)
  END FUNCTION DERIV
  REAL FUNCTION ANTDER(ax)
    REAL, INTENT(IN) :: ax
    ANTDER = ax**3*(thrqtr*ax-two)
  END FUNCTION ANTDER
END SUBROUTINE PCHQK2
