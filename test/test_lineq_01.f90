MODULE TEST21_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  !** SGEQC
  SUBROUTINE SGEQC(Lun,Kprint,Nerr)
    !> Quick check for SGEFS and SGEIR.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      SINGLE PRECISION (SGEQC-S, DGEQC-D, CGEQC-C)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Jacobsen, Nancy, (LANL)
    !***
    ! **Description:**
    !
    !   Let A*X=B be a SINGLE PRECISION linear system where the
    !   matrix is of the proper type for the Linpack subroutines
    !   being called.  The values of A and B and the pre-computed
    !   values of BXEX (the solution vector) are given in DATA
    !   statements.  The computed test results for X are compared to
    !   the stored pre-computed values.  Failure of the test occurs
    !   when there is less than 80% agreement between the absolute
    !   values.  There are 2 tests - one for the normal case and one
    !   for the singular case.  A message is printed indicating
    !   whether each subroutine has passed or failed for each case.
    !
    !   On return, NERR (INTEGER type) contains the total count of
    !   all failures detected.
    !
    !***
    ! **Routines called:**  R1MACH, SGEFS, SGEIR

    !* REVISION HISTORY  (YYMMDD)
    !   801022  DATE WRITTEN
    !   891009  Removed unreferenced statement label.  (WRB)
    !   891009  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   920601  Code reworked and TYPE section added.  (RWC, WRB)
    USE slatec, ONLY : eps_sp, SGEFS, SGEIR
    !     .. Scalar Arguments ..
    INTEGER :: Kprint, Lun, Nerr
    !     .. Local Scalars ..
    REAL(SP) :: errcmp, errmax
    INTEGER :: i, ind, itask, j, kprog, lda, n
    !     .. Local Arrays ..
    REAL(SP) :: atemp(5,4), btemp(4), work(20)
    INTEGER :: iwork(4)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX
    !     .. Data statements ..
    REAL(SP), PARAMETER :: a(5,4) = RESHAPE( [ 5._SP, 1._SP, 0.3_SP, 2.1_SP, 0._SP, &
      -1._SP, -0.5_SP, 1._SP, 1._SP, 0._SP, 4.5_SP, -1._SP, -1.7_SP, 2._SP, 0._SP, &
      0.5_SP, 2._SP, 0.6_SP, 1.3_SP, 0._SP ], [5,4] )
    REAL(SP), PARAMETER :: b(4) = [ 0._SP, 3.5_SP, 3.6_SP, 2.4_SP ]
    REAL(SP), PARAMETER :: bxex(4) = [ 1._SP, 1._SP, -1._SP, 1._SP ]
    CHARACTER(4), PARAMETER :: list(2) = [ 'GEFS', 'GEIR' ]
    !* FIRST EXECUTABLE STATEMENT  SGEQC
    n = 4
    lda = 5
    Nerr = 0
    errcmp = eps_sp**0.8_SP
    IF( Kprint>=2 ) WRITE (Lun,99001)
    !
    99001 FORMAT (//,2X,'SGEFS and SGEIR Quick Check'/)
    !
    DO kprog = 1, 2
      !
      !     First test case - normal
      !
      itask = 1
      DO i = 1, n
        btemp(i) = b(i)
      END DO
      DO j = 1, n
        DO i = 1, n
          atemp(i,j) = a(i,j)
        END DO
      END DO
      IF( kprog==1 ) THEN
        CALL SGEFS(atemp,lda,n,btemp,itask,ind,work,iwork)
      ELSE
        CALL SGEIR(atemp,lda,n,btemp,itask,ind,work,iwork)
      END IF
      IF( ind<0 ) THEN
        IF( Kprint>=2 ) WRITE (Lun,FMT=99004) list(kprog), ind
        Nerr = Nerr + 1
      END IF
      !
      !       Calculate error for first test
      !
      errmax = 0._SP
      !
      DO i = 1, n
        errmax = MAX(errmax,ABS(btemp(i)-bxex(i)))
      END DO
      IF( errcmp>errmax ) THEN
        IF( Kprint>=3 ) WRITE (Lun,FMT=99004) list(kprog)
      ELSE
        IF( Kprint>=2 ) WRITE (Lun,FMT=99005) list(kprog), errmax
        Nerr = Nerr + 1
      END IF
      !
      !       Second test case - singular matrix
      !
!      itask = 1
!      DO i = 1, n
!        btemp(i) = b(i)
!      END DO
!      DO j = 1, n
!        DO i = 1, n
!          atemp(i,j) = a(i,j)
!        END DO
!      END DO
!      DO j = 1, n
!        atemp(1,j) = 0._SP
!      END DO
!      IF( kprog==1 ) THEN
!        CALL SGEFS(atemp,lda,n,btemp,itask,ind,work,iwork)
!      ELSE
!        CALL SGEIR(atemp,lda,n,btemp,itask,ind,work,iwork)
!      END IF
!      IF( ind==-4 ) THEN
!        IF( Kprint>=3 ) WRITE (Lun,FMT=99006) list(kprog)
!      ELSE
!        IF( Kprint>=2 ) WRITE (Lun,FMT=99007) list(kprog), ind
!        Nerr = Nerr + 1
!      END IF
      !
    END DO
    !
    IF( Kprint>=3 .AND. Nerr==0 ) WRITE (Lun,99002)
    99002 FORMAT (/,2X,'SGEFS and SGEIR Quick Check PASSED'/)
    IF( Kprint>=2 .AND. Nerr/=0 ) WRITE (Lun,99003)
    99003 FORMAT (/,2X,'SGEFS and SGEIR Quick Check FAILED'/)
    RETURN
    99004 FORMAT (/,5X,'S',A,' Normal test PASSED')
    99005 FORMAT (/,5X,'S',A,' Test FAILED, MAX ABS(ERROR) is',E13.5)
!    99006 FORMAT (/,5X,'S',A,' Singular test PASSED')
!    99007 FORMAT (/,5X,'S',A,' Singular test FAILED, IND=',I3)
  END SUBROUTINE SGEQC
  !** DGEQC
  SUBROUTINE DGEQC(Lun,Kprint,Nerr)
    !> Quick check for DGEFS.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (SGEQC-S, DGEQC-D, CGEQC-C)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Jacobsen, Nancy, (LANL)
    !***
    ! **Description:**
    !
    !   Let A*X=B be a DOUBLE PRECISION linear system where the
    !   matrix is of the proper type for the Linpack subroutines
    !   being called.  The values of A and B and the pre-computed
    !   values of BXEX (the solution vector) are given in DATA
    !   statements.  The computed test results for X are compared to
    !   the stored pre-computed values.  Failure of the test occurs
    !   when there is less than 80% agreement between the absolute
    !   values.  There are 2 tests - one for the normal case and one
    !   for the singular case.  A message is printed indicating
    !   whether each subroutine has passed or failed for each case.
    !
    !   On return, NERR (INTEGER type) contains the total count of
    !   all failures detected.
    !
    !***
    ! **Routines called:**  D1MACH, DGEFS

    !* REVISION HISTORY  (YYMMDD)
    !   801022  DATE WRITTEN
    !   891009  Removed unreferenced statement label.  (WRB)
    !   891009  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   920601  Code reworked and TYPE section added.  (RWC, WRB)
    USE slatec, ONLY : eps_dp, DGEFS
    !     .. Scalar Arguments ..
    INTEGER :: Kprint, Lun, Nerr
    !     .. Local Scalars ..
    REAL(DP) :: errcmp, errmax
    INTEGER :: i, ind, itask, j, kprog, lda, n
    !     .. Local Arrays ..
    REAL(DP) :: atemp(5,4), btemp(4), work(20)
    INTEGER :: iwork(4)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX
    !     .. Data statements ..
    REAL(DP), PARAMETER :: a(5,4) = RESHAPE( [ 5._DP, 1._DP, 0.3_DP, 2.1_DP, 0._DP, &
      -1._DP, -0.5_DP, 1._DP, 1._DP, 0._DP, &
      4.5_DP,  -1._DP, -1.7_DP, 2._DP, 0._DP, &
      0.5_DP, 2._DP, 0.6_DP, 1.3_DP, 0._DP ], [5,4] )
    REAL(DP), PARAMETER :: b(4) = [ 0._DP, 3.5_DP, 3.6_DP, 2.4_DP ]
    REAL(DP), PARAMETER :: bxex(4) = [ 1._DP, 1._DP, -1._DP, 1._DP ]
    CHARACTER(4), PARAMETER :: list(2) = [ 'GEFS', 'GEIR' ]
    !* FIRST EXECUTABLE STATEMENT  DGEQC
    n = 4
    lda = 5
    Nerr = 0
    errcmp = eps_dp**0.8_DP
    IF( Kprint>=2 ) WRITE (Lun,99001)
    !
    99001 FORMAT (//,2X,'DGEFS Quick Check'/)
    !
    kprog = 1
    !
    !     First test case - normal
    !
    itask = 1
    DO i = 1, n
      btemp(i) = b(i)
    END DO
    DO j = 1, n
      DO i = 1, n
        atemp(i,j) = a(i,j)
      END DO
    END DO
    CALL DGEFS(atemp,lda,n,btemp,itask,ind,work,iwork)
    IF( ind<0 ) THEN
      IF( Kprint>=2 ) WRITE (Lun,FMT=99005) list(kprog), ind
      Nerr = Nerr + 1
    END IF
    !
    !     Calculate error for first test
    !
    errmax = 0._DP
    !
    DO i = 1, n
      errmax = MAX(errmax,ABS(btemp(i)-bxex(i)))
    END DO
    IF( errcmp>errmax ) THEN
      IF( Kprint>=3 ) WRITE (Lun,FMT=99004) list(kprog)
    ELSE
      IF( Kprint>=2 ) WRITE (Lun,FMT=99005) list(kprog), errmax
      Nerr = Nerr + 1
    END IF
    !
    !     Second test case - singular matrix
    !
!    itask = 1
!    DO i = 1, n
!      btemp(i) = b(i)
!    END DO
!    DO j = 1, n
!      DO i = 1, n
!        atemp(i,j) = a(i,j)
!      END DO
!    END DO
!    DO j = 1, n
!      atemp(1,j) = 0._DP
!    END DO
!    CALL DGEFS(atemp,lda,n,btemp,itask,ind,work,iwork)
!    IF( ind==-4 ) THEN
!      IF( Kprint>=3 ) WRITE (Lun,FMT=99006) list(kprog)
!    ELSE
!      IF( Kprint>=2 ) WRITE (Lun,FMT=99007) list(kprog), ind
!      Nerr = Nerr + 1
!    END IF
    !
    IF( Kprint>=3 .AND. Nerr==0 ) WRITE (Lun,99002)
    99002 FORMAT (/,2X,'DGEFS Quick Check PASSED'/)
    IF( Kprint>=2 .AND. Nerr/=0 ) WRITE (Lun,99003)
    99003 FORMAT (/,2X,'DGEFS Quick Check FAILED'/)
    RETURN
    99004 FORMAT (/,5X,'D',A,' Normal test PASSED')
    99005 FORMAT (/,5X,'D',A,' Test FAILED, MAX ABS(ERROR) is',E13.5)
!    99006 FORMAT (/,5X,'D',A,' Singular test PASSED')
!    99007 FORMAT (/,5X,'D',A,' Singular test FAILED, IND=',I3)
  END SUBROUTINE DGEQC
  !** CGEQC
  SUBROUTINE CGEQC(Lun,Kprint,Nerr)
    !> Quick check for CGEFS and CGEIR.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      COMPLEX (SGEQC-S, DGEQC-D, CGEQC-C)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Jacobsen, Nancy, (LANL)
    !***
    ! **Description:**
    !
    !   Let A*X=B be a COMPLEX linear system where the
    !   matrix is of the proper type for the Linpack subroutines
    !   being called.  The values of A and B and the pre-computed
    !   values of BXEX (the solution vector) are given in DATA
    !   statements.  The computed test results for X are compared to
    !   the stored pre-computed values.  Failure of the test occurs
    !   when there is less than 80% agreement between the absolute
    !   values.  There are 2 tests - one for the normal case and one
    !   for the singular case.  A message is printed indicating
    !   whether each subroutine has passed or failed for each case.
    !
    !   On return, NERR (INTEGER type) contains the total count of
    !   all failures detected.
    !
    !***
    ! **Routines called:**  CGEFS, CGEIR

    !* REVISION HISTORY  (YYMMDD)
    !   801029  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   920601  Code reworked and TYPE section added.  (RWC, WRB)
    USE slatec, ONLY : CGEFS, CGEIR
    USE blas, ONLY : SCABS1
    !     .. Scalar Arguments ..
    INTEGER :: Kprint, Lun, Nerr
    !     .. Local Scalars ..
    INTEGER :: i, ind, indx, itask, j, kprog, lda, n
    !     .. Local Arrays ..
    COMPLEX(SP) :: atemp(5,3), btemp(3), work(12)
    INTEGER :: iwork(3)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, AIMAG, REAL
    !     .. Data statements ..
    COMPLEX(SP), PARAMETER :: a(3,3) = RESHAPE( [ &
      (2._SP,3._SP), (1._SP,1._SP), (1._SP,2._SP), &
      (2._SP,0._SP), (1._SP,-1._SP), (0._SP,0._SP), &
      (0._SP,0._SP), (2._SP,5._SP), (3._SP,2._SP) ], [3,3] )
    COMPLEX(SP), PARAMETER :: b(3) = [ (-1._SP,1._SP), (-5._SP,4._SP), (-4._SP,7._SP) ]
    COMPLEX(SP), PARAMETER :: bxex(3) = [ (.21459E-01_SP,.209012E+01_SP), &
      (.261373E+01_SP,-.162231E+01_SP), (.785407E+00_SP,.109871E+01_SP) ]
    CHARACTER(4), PARAMETER :: list(2) = [ 'GEFS', 'GEIR' ]
    !* FIRST EXECUTABLE STATEMENT  CGEQC
    n = 3
    lda = 5
    Nerr = 0
    IF( Kprint>=2 ) WRITE (Lun,99001)
    !
    99001 FORMAT (//,2X,'CGEFS and CGEIR Quick Check'/)
    !
    DO kprog = 1, 2
      !
      !     First test case - normal
      !
      itask = 1
      DO i = 1, n
        btemp(i) = b(i)
      END DO
      DO j = 1, n
        DO i = 1, n
          atemp(i,j) = a(i,j)
        END DO
      END DO
      IF( kprog==1 ) THEN
        CALL CGEFS(atemp,lda,n,btemp,itask,ind,work,iwork)
      ELSE
        CALL CGEIR(atemp,lda,n,btemp,itask,ind,work,iwork)
      END IF
      IF( ind<0 ) THEN
        IF( Kprint>=2 ) WRITE (Lun,FMT=99005) list(kprog), ind
        Nerr = Nerr + 1
      END IF
      !
      !       Calculate error for first test
      !
      indx = 0
      DO i = 1, n
        IF( SCABS1(bxex(i)-btemp(i))>.0001 ) indx = indx + 1
      END DO
      IF( indx==0 ) THEN
        IF( Kprint>=3 ) WRITE (Lun,FMT=99004) list(kprog)
      ELSE
        IF( Kprint>=2 ) WRITE (Lun,FMT=99005) list(kprog)
        Nerr = Nerr + 1
      END IF
      !
      !       Second test case - singular matrix
      !
!      itask = 1
!      DO i = 1, n
!        btemp(i) = b(i)
!      END DO
!      DO j = 1, n
!        DO i = 1, n
!          atemp(i,j) = a(i,j)
!        END DO
!      END DO
!      DO j = 1, n
!        atemp(1,j) = (0._SP,0._SP)
!      END DO
!      IF( kprog==1 ) THEN
!        CALL CGEFS(atemp,lda,n,btemp,itask,ind,work,iwork)
!      ELSE
!        CALL CGEIR(atemp,lda,n,btemp,itask,ind,work,iwork)
!      END IF
!      IF( ind==-4 ) THEN
!        IF( Kprint>=3 ) WRITE (Lun,FMT=99006) list(kprog)
!      ELSE
!        IF( Kprint>=2 ) WRITE (Lun,FMT=99007) list(kprog), ind
!        Nerr = Nerr + 1
!      END IF
    END DO
    !
    IF( Kprint>=3 .AND. Nerr==0 ) WRITE (Lun,99002)
    99002 FORMAT (/,2X,'CGEFS and CGEIR Quick Check PASSED'/)
    IF( Kprint>=2 .AND. Nerr/=0 ) WRITE (Lun,99003)
    99003 FORMAT (/,2X,'CGEFS and CGEIR Quick Check FAILED'/)
    RETURN
    99004 FORMAT (/,5X,'C',A,' Normal test PASSED')
    99005 FORMAT (/,5X,'C',A,' Test FAILED')
!    99006 FORMAT (/,5X,'C',A,' Singular test PASSED')
!    99007 FORMAT (/,5X,'C',A,' Singular test FAILED, IND=',I3)
  END SUBROUTINE CGEQC
END MODULE TEST21_MOD
!** TEST21
PROGRAM TEST21
  USE TEST21_MOD, ONLY : CGEQC, DGEQC, SGEQC
  USE ISO_FORTRAN_ENV, ONLY : INPUT_UNIT, OUTPUT_UNIT
  USE slatec, ONLY : control_xer, max_xer
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !> Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  D2
  !***
  ! **Type:**      ALL (TEST21-A)
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
  !        SGEFS    SGEIR
  !        DGEFS
  !        CGEFS    CGEIR
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  CGEQC, DGEQC, I1MACH, SGEQC, XERMAX, XSETF

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  INTEGER :: kprint, lin, lun, nerr, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST21
  lun = OUTPUT_UNIT
  lin = INPUT_UNIT
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
  !     Test LINPACK routines
  !
  CALL SGEQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL DGEQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CGEQC(lun,kprint,nerr)
  nfail = nfail + nerr
  !
  !     Write PASS or FAIL message
  !
  IF( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST21 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST21 *************')
  END IF
  STOP
END PROGRAM TEST21
