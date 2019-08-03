MODULE TEST29_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  !** PNTCHK
  SUBROUTINE PNTCHK(Lun,Kprint,Ipass)
    !> Quick check for POLINT, POLCOF and POLYVL
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      SINGLE PRECISION (PNTCHK-S, DPNTCK-D)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  NUMXER, POLCOF, POLINT, POLYVL, R1MACH, XERCLR,
    !                    XGETF, XSETF

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Changed usage of eps_2_sp to eps_sp.  (RWC)
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    !   920212  Code completely restructured to test errors for all values of KPRINT.  (WRB)
    USE service, ONLY : eps_sp
    USE interpolation, ONLY : POLCOF, POLINT, POLYVL
    !     .. Scalar Arguments ..
    INTEGER :: Ipass, Kprint, Lun
    !     .. Local Scalars ..
    REAL(SP) :: tol, yf
    INTEGER :: i, ierr, n
    LOGICAL :: fatal
    !     .. Local Arrays ..
    REAL(SP) :: c(6), d(6), w(12)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, SQRT
    !     .. Data statements ..
    REAL(SP) :: x(6) = [ 1._SP, 2._SP, 3._SP, -1._SP, -2._SP, -3._SP ]
    REAL(SP), PARAMETER :: y(6) = [ 0._SP, 9._SP, 64._SP, 0._SP, 9._SP, 64._SP ]
    REAL(SP), PARAMETER :: xchk(6) = [ 1._SP, 0._SP, -2._SP, 0._SP, 1._SP, 0._SP ]
    REAL(SP), PARAMETER :: dchk(6) = [ 1._SP, 0._SP, -4._SP, 0._SP, 24._SP, 0._SP ]
    !* FIRST EXECUTABLE STATEMENT  PNTCHK
    IF( Kprint>=2 ) WRITE (Lun,99001)
    !
    99001 FORMAT ('1'/' Test POLINT, POLCOF and POLYVL')
    !
    !     Initialize variables for tests.
    !
    tol = SQRT(eps_sp)
    Ipass = 1
    n = 6
    !
    !     Set up polynomial test.
    !
    CALL POLINT(n,x,y,c)
    CALL POLCOF(0._SP,n,x,c,d,w)
    !
    !     Check to see if POLCOF test passed.
    !
    fatal = .FALSE.
    DO i = 1, n
      IF( ABS(d(i)-xchk(i))>tol ) THEN
        Ipass = 0
        fatal = .TRUE.
      END IF
    END DO
    IF( fatal ) THEN
      IF( Kprint>=2 ) WRITE (Lun,99007) 'FAILED', (d(i),i=1,n)
    ELSE
      IF( Kprint>=3 ) WRITE (Lun,99007) 'PASSED', (d(i),i=1,n)
    END IF
    !
    !     Test POLYVL.
    !
    CALL POLYVL(5,0._SP,yf,d,n,x,c,w,ierr)
    IF( ABS(dchk(1)-yf)<=tol ) THEN
      IF( Kprint>=3 ) WRITE (Lun,99008) 'PASSED', yf, (d(i),i=1,5)
    ELSE
      Ipass = 0
      IF( Kprint>=2 ) WRITE (Lun,99008) 'FAILED', yf, (d(i),i=1,5)
    END IF
    !
    fatal = .FALSE.
    DO i = 1, 5
      IF( ABS(dchk(i+1)-d(i))>tol ) THEN
        Ipass = 0
        fatal = .TRUE.
      END IF
    END DO
    !
    !     Trigger 2 error conditions
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
!    IF( Kprint>=3 ) WRITE (Lun,99002)
!    99002 FORMAT (/' 2 Error messages expected')
!    CALL POLINT(0,x,y,c)
!    IF( num_xer/=2 ) THEN
!      Ipass = 0
!      fatal = .TRUE.
!    END IF
!    num_xer = 0
!    !
!    x(1) = -1._SP
!    CALL POLINT(n,x,y,c)
!    IF( num_xer/=2 ) THEN
!      Ipass = 0
!      fatal = .TRUE.
!    END IF
!    num_xer = 0
!    !
!    control_xer = kontrl
!    IF( fatal ) THEN
!      IF( Kprint>=2 ) THEN
!        WRITE (Lun,99003)
!        99003 FORMAT (/' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
!      END IF
!    ELSEIF( Kprint>=3 ) THEN
!      WRITE (Lun,99004)
!      99004 FORMAT (/' ALL INCORRECT ARGUMENT TESTS PASSED')
!    END IF
    !
    IF( Ipass==1 .AND. Kprint>=2 ) WRITE (Lun,99005)
    99005 FORMAT (/' ****************POLINT PASSED ALL TESTS**************')
    IF( Ipass==0 .AND. Kprint>=1 ) WRITE (Lun,99006)
    99006 FORMAT (/' ***************POLINT FAILED SOME TESTS**************')
    RETURN
    99007 FORMAT (/'POLCOF ',A,&
      ' test'/' Taylor coefficients for the quintic should be'/6X,&
      '1.000',5X,'0.000',4X,'-2.000',5X,'0.000',5X,'1.000',5X,&
      '0.000'/' Taylor coefficients from POLCOF are'/1X,6F10.3/)
    99008 FORMAT (' Derivative test ',&
      A/' The derivatives of the polynomial at zero as ',&
      'computed by POLYVL are'/1X,6F10.3/)
  END SUBROUTINE PNTCHK
  !** DPNTCK
  SUBROUTINE DPNTCK(Lun,Kprint,Ipass)
    !> Quick check for DPLINT, DPOLCF and DPOLVL
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (PNTCHK-S, DPNTCK-D)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Boland, W. Robert, (LANL)
    !***
    ! **Routines called:**  D1MACH, DPLINT, DPOLCF, DPOLVL, NUMXER, XERCLR,
    !                    XGETF, XSETF

    !* REVISION HISTORY  (YYMMDD)
    !   920212  DATE WRITTEN
    USE service, ONLY : eps_dp
    USE interpolation, ONLY : DPLINT, DPOLCF, DPOLVL
    !     .. Scalar Arguments ..
    INTEGER :: Ipass, Kprint, Lun
    !     .. Local Scalars ..
    REAL(DP) :: tol, yf
    INTEGER :: i, ierr, n
    LOGICAL :: fatal
    !     .. Local Arrays ..
    REAL(DP) :: c(6), d(6), w(12)
    !     .. Data statements ..
    REAL(DP) :: x(6) = [ 1._DP, 2._DP, 3._DP, -1._DP, -2._DP, -3._DP ]
    REAL(DP), PARAMETER :: y(6) = [ 0._DP, 9._DP, 64._DP, 0._DP, 9._DP, 64._DP ]
    REAL(DP), PARAMETER :: xchk(6) = [ 1._DP, 0._DP, -2._DP, 0._DP, 1._DP, 0._DP ]
    REAL(DP), PARAMETER :: dchk(6) = [ 1._DP, 0._DP, -4._DP, 0._DP, 24._DP, 0._DP ]
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, SQRT
    !* FIRST EXECUTABLE STATEMENT  DPNTCK
    IF( Kprint>=2 ) WRITE (Lun,99001)
    !
    99001 FORMAT ('1'/' Test DPLINT, DPOLCF and DPOLVL')
    !
    !     Initialize variables for tests.
    !
    tol = SQRT(eps_dp)
    Ipass = 1
    n = 6
    !
    !     Set up polynomial test.
    !
    CALL DPLINT(n,x,y,c)
    CALL DPOLCF(0._DP,n,x,c,d,w)
    !
    !     Check to see if DPOLCF test passed.
    !
    fatal = .FALSE.
    DO i = 1, n
      IF( ABS(d(i)-xchk(i))>tol ) THEN
        Ipass = 0
        fatal = .TRUE.
      END IF
    END DO
    IF( fatal ) THEN
      IF( Kprint>=2 ) WRITE (Lun,99007) 'FAILED', (d(i),i=1,n)
    ELSE
      IF( Kprint>=3 ) WRITE (Lun,99007) 'PASSED', (d(i),i=1,n)
    END IF
    !
    !     Test DPOLVL.
    !
    CALL DPOLVL(5,0._DP,yf,d,n,x,c,w,ierr)
    IF( ABS(dchk(1)-yf)<=tol ) THEN
      IF( Kprint>=3 ) WRITE (Lun,99008) 'PASSED', yf, (d(i),i=1,5)
    ELSE
      Ipass = 0
      IF( Kprint>=2 ) WRITE (Lun,99008) 'FAILED', yf, (d(i),i=1,5)
    END IF
    !
    fatal = .FALSE.
    DO i = 1, 5
      IF( ABS(dchk(i+1)-d(i))>tol ) THEN
        Ipass = 0
        fatal = .TRUE.
      END IF
    END DO
    !
    !     Trigger 2 error conditions
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
!    IF( Kprint>=3 ) WRITE (Lun,99002)
!    99002 FORMAT (/' 2 Error messages expected')
!    CALL DPLINT(0,x,y,c)
!    IF( num_xer/=2 ) THEN
!      Ipass = 0
!      fatal = .TRUE.
!    END IF
!    num_xer = 0
!    !
!    x(1) = -1._DP
!    CALL DPLINT(n,x,y,c)
!    IF( num_xer/=2 ) THEN
!      Ipass = 0
!      fatal = .TRUE.
!    END IF
!    num_xer = 0
!    !
!    control_xer = kontrl
!    IF( fatal ) THEN
!      IF( Kprint>=2 ) THEN
!        WRITE (Lun,99003)
!        99003 FORMAT (/' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
!      END IF
!    ELSEIF( Kprint>=3 ) THEN
!      WRITE (Lun,99004)
!      99004 FORMAT (/' ALL INCORRECT ARGUMENT TESTS PASSED')
!    END IF
    !
    IF( Ipass==1 .AND. Kprint>=2 ) WRITE (Lun,99005)
    99005 FORMAT (/' ****************DPLINT PASSED ALL TESTS**************')
    IF( Ipass==0 .AND. Kprint>=1 ) WRITE (Lun,99006)
    99006 FORMAT (/' ***************DPLINT FAILED SOME TESTS**************')
    RETURN
    99007 FORMAT (/'DPOLCF ',A,&
      ' test'/' Taylor coefficients for the quintic should be'/6X,&
      '1.000',5X,'0.000',4X,'-2.000',5X,'0.000',5X,'1.000',5X,&
      '0.000'/' Taylor coefficients from DPOLCF are'/1X,6F10.3/)
    99008 FORMAT (' Derivative test ',&
      A/' The derivatives of the polynomial at zero as ',&
      'computed by DPOLVL are'/1X,6F10.3/)
  END SUBROUTINE DPNTCK
END MODULE TEST29_MOD
!** TEST29
PROGRAM TEST29
  USE TEST29_MOD, ONLY : DPNTCK, PNTCHK
  USE ISO_FORTRAN_ENV, ONLY : INPUT_UNIT, OUTPUT_UNIT
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !> Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  E1, E3
  !***
  ! **Type:**      ALL (TEST29-A)
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
  !        POLINT   POLCOF   POLYVL
  !        DPLINT   DPOLCF   DPOLVL
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  DPNTCK, I1MACH, PNTCHK, XERMAX, XSETF

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  !   920225  Added CALL to DPNTCK.  (WRB)
  INTEGER :: ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST29
  lun = OUTPUT_UNIT
  lin = INPUT_UNIT
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  !
  !     Test POLINT, POLCOF and POLYVL.
  !
  CALL PNTCHK(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Test DPLINT, DPOLCF and DPOLVL.
  !
  CALL DPNTCK(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST29 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST29 *************')
  END IF
  STOP
END PROGRAM TEST29
