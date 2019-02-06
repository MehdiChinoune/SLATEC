!*==DPCHQ2.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DPCHQ2
      SUBROUTINE DPCHQ2(Lun,Kprint,Ipass)
      IMPLICIT NONE
!*--DPCHQ25
!***BEGIN PROLOGUE  DPCHQ2
!***PURPOSE  Test the PCHIP integrators DPCHIA and DPCHID.
!***LIBRARY   SLATEC (PCHIP)
!***TYPE      DOUBLE PRECISION (PCHQK2-S, DPCHQ2-D)
!***KEYWORDS  PCHIP INTEGRATOR QUICK CHECK
!***AUTHOR  Fritsch, F. N., (LLNL)
!***DESCRIPTION
!
!             DPCHIP QUICK CHECK NUMBER 2
!
!     TESTS THE INTEGRATORS:  DPCHIA, DPCHID.
! *Usage:
!
!        INTEGER  LUN, KPRINT, IPASS
!
!        CALL DPCHQ2 (LUN, KPRINT, IPASS)
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
!   This routine constructs data from a cubic, integrates it with DPCHIA
!   and compares the results with the correct answer.
!   Since DPCHIA calls DPCHID, this tests both integrators.
!
!***ROUTINES CALLED  D1MACH, DPCHIA
!***REVISION HISTORY  (YYMMDD)
!   820601  DATE WRITTEN
!   890306  Changed IPASS to the more accurate name IFAIL.  (FNF)
!   890316  1. Removed IMPLICIT statement.                  (FNF)
!           2. Eliminated unnecessary variable N1.          (FNF)
!           3. Miscellaneous cosmetic changes.              (FNF)
!   891004  Cosmetic changes to prologue.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900314  Improved some output formats.  (FNF)
!   900315  Revised prologue and improved some output formats.  (FNF)
!   900316  Additional minor cosmetic changes.  (FNF)
!   900321  Removed IFAIL from call sequence for SLATEC standards and
!           made miscellaneous cosmetic changes.  (FNF)
!   900323  Corrected list of routines called.  (FNF)
!   901130  Added 1P's to formats; changed to allow KPRINT.gt.3.  (FNF)
!   910708  Minor modifications in use of KPRINT.  (WRB)
!   930317  Improved output formats.  (FNF)
!***END PROLOGUE  DPCHQ2
!
!  Declare arguments.
!
      INTEGER Lun , Kprint , Ipass
!
!  DECLARE VARIABLES.
!
      INTEGER i , ierexp(17) , ierr , ifail , n , npairs
      DOUBLE PRECISION a(17) , b(17) , calc , d(7) , errmax , error , f(7) ,
     &                 machep , one , three , thrqtr , tol , true , two , x(7)
      LOGICAL fail , skip
!
!  DECLARE EXTERNALS.
!
      DOUBLE PRECISION DPCHIA , D1MACH
!
!  INITIALIZE.
!
      DATA thrqtr/0.75D0/ , one/1.D0/ , two/2.D0/ , three/3.D0/
      DATA n/7/
      DATA x/ - 4.D0 , -2.D0 , -0.9D0 , 0.D0 , 0.9D0 , 2.D0 , 4.D0/
      DATA npairs/17/
      DATA a/ - 3.0D0 , 3.0D0 , -0.5D0 , -0.5D0 , -0.5D0 , -4.0D0 , -4.0D0 ,
     &     3.0D0 , -5.0D0 , -5.0D0 , -6.0D0 , 6.0D0 , -1.5D0 , -1.5D0 , -3.0D0 ,
     &     3.0D0 , 0.5D0/
      DATA b/3.0D0 , -3.0D0 , 1.0D0 , 2.0D0 , 5.0D0 , -0.5D0 , 4.0D0 , 5.0D0 ,
     &     -3.0D0 , 5.0D0 , -5.0D0 , 5.0D0 , -0.5D0 , -1.0D0 , -2.5D0 , 3.5D0 ,
     &     0.5D0/
      DATA ierexp/0 , 0 , 0 , 0 , 2 , 0 , 0 , 2 , 1 , 3 , 3 , 3 , 0 , 0 , 0 ,
     &     0 , 0/
!
!  SET PASS/FAIL TOLERANCE.
!
!***FIRST EXECUTABLE STATEMENT  DPCHQ2
      machep = D1MACH(4)
      tol = 100.D0*machep
!
!  SET UP PCH FUNCTION DEFINITION.
!
      DO i = 1 , n
        f(i) = FCN(x(i))
        d(i) = DERIV(x(i))
      ENDDO
!
      IF ( Kprint>=3 ) WRITE (Lun,99001)
!
!  FORMATS.
!
99001 FORMAT ('1'//10X,'TEST DPCHIP INTEGRATORS')
      IF ( Kprint>=2 ) WRITE (Lun,99002)
99002 FORMAT (//10X,'DPCHQ2 RESULTS'/10X,'--------------')
      IF ( Kprint>=3 ) WRITE (Lun,99003) (x(i),f(i),d(i),i=1,n)
99003 FORMAT (//5X,'DATA:'//11X,'X',9X,'F',9X,'D'/(5X,3F10.3))
!
!  LOOP OVER (A,B)-PAIRS.
!
      IF ( Kprint>=3 ) WRITE (Lun,99004)
99004 FORMAT (//5X,'TEST RESULTS:'//'    A     B    ERR     TRUE',16X,'CALC',
     &        15X,'ERROR')
!
      ifail = 0
!
      skip = .FALSE.
      DO i = 1 , npairs
!               ---------------------------------------------
        calc = DPCHIA(n,x,f,d,1,skip,a(i),b(i),ierr)
!               ---------------------------------------------
        IF ( ierr>=0 ) THEN
          fail = ierr/=ierexp(i)
          true = ANTDER(b(i)) - ANTDER(a(i))
          error = calc - true
          IF ( Kprint>=3 ) THEN
            IF ( fail ) THEN
              WRITE (Lun,99005) a(i) , b(i) , ierr , true , calc , error ,
     &                          ierexp(i)
99005         FORMAT (2F6.1,I5,1P,2D20.10,D15.5,'  (',I1,') *****')
            ELSE
              WRITE (Lun,99010) a(i) , b(i) , ierr , true , calc , error
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
          IF ( Kprint>=3 ) WRITE (Lun,99010) a(i) , b(i) , ierr
          ifail = ifail + 1
        ENDIF
      ENDDO
!
!  PRINT SUMMARY.
!
      IF ( Kprint>=2 ) THEN
        WRITE (Lun,99006) errmax , tol
99006   FORMAT (/'  MAXIMUM RELATIVE ERROR IS:',1P,D15.5,',   TOLERANCE:',1P,
     &          D15.5)
        IF ( ifail/=0 ) WRITE (Lun,99007) ifail
99007   FORMAT (/' *** TROUBLE ***',I5,' INTEGRATION TESTS FAILED.')
      ENDIF
!
!  TERMINATE.
!
      IF ( ifail==0 ) THEN
        Ipass = 1
        IF ( Kprint>=2 ) WRITE (Lun,99008)
99008   FORMAT (/' ------------ DPCHIP PASSED  ALL INTEGRATION TESTS',
     &          ' ------------')
      ELSE
        Ipass = 0
        IF ( Kprint>=1 ) WRITE (Lun,99009)
99009   FORMAT (/' ************ DPCHIP FAILED SOME INTEGRATION TESTS',
     &          ' ************')
      ENDIF
!
      RETURN
99010 FORMAT (2F6.1,I5,1P,2D20.10,D15.5)
!------------- LAST LINE OF DPCHQ2 FOLLOWS -----------------------------
      CONTAINS

!
!  DEFINE TEST FUNCTIONS.
!
        REAL(8) FUNCTION FCN(ax)
          REAL(8), INTENT(IN) :: ax
          FCN = three*ax*ax*(ax-two)
        END FUNCTION FCN
        REAL(8) FUNCTION DERIV(ax)
          REAL(8), INTENT(IN) :: ax
          DERIV = three*ax*(two*(ax-two)+ax)
        END FUNCTION DERIV
        REAL(8) FUNCTION ANTDER(ax)
          REAL(8), INTENT(IN) :: ax
          ANTDER = ax**3*(thrqtr*ax-two)
        END FUNCTION ANTDER
      END SUBROUTINE DPCHQ2
