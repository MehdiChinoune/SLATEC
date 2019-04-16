MODULE TEST21_MOD
  use slatec
  IMPLICIT NONE

CONTAINS
  !** SGEQC
  SUBROUTINE SGEQC(Lun,Kprint,Nerr)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for SGEFS and SGEIR.
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

    !     .. Scalar Arguments ..
    INTEGER Kprint, Lun, Nerr
    !     .. Local Scalars ..
    REAL errcmp, errmax
    INTEGER i, ind, itask, j, kprog, lda, n
    !     .. Local Arrays ..
    REAL atemp(5,4), btemp(4), work(20)
    INTEGER iwork(4)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX
    !     .. Data statements ..
    REAL, PARAMETER :: a(5,4) = RESHAPE( [ 5.0E0, 1.0E0, 0.3E0, 2.1E0, 0.0E0, &
      -1.0E0, -0.5E0, 1.0E0, 1.0E0, 0.0E0, 4.5E0, -1.0E0, -1.7E0, 2.0E0, 0.0E0, &
      0.5E0, 2.0E0, 0.6E0, 1.3E0, 0.0E0 ], [5,4] )
    REAL, PARAMETER :: b(4) = [ 0.0E0, 3.5E0, 3.6E0, 2.4E0 ]
    REAL, PARAMETER :: bxex(4) = [ 0.10E+01, 0.10E+01, -0.10E+01, 0.10E+01 ]
    CHARACTER(4), PARAMETER :: list(2) = [ 'GEFS', 'GEIR' ]
    !* FIRST EXECUTABLE STATEMENT  SGEQC
    n = 4
    lda = 5
    Nerr = 0
    errcmp = R1MACH(4)**0.8E0
    IF ( Kprint>=2 ) WRITE (Lun,99001)
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
      IF ( kprog==1 ) THEN
        CALL SGEFS(atemp,lda,n,btemp,itask,ind,work,iwork)
      ELSE
        CALL SGEIR(atemp,lda,n,btemp,itask,ind,work,iwork)
      END IF
      IF ( ind<0 ) THEN
        IF ( Kprint>=2 ) WRITE (Lun,FMT=99004) list(kprog), ind
        Nerr = Nerr + 1
      END IF
      !
      !       Calculate error for first test
      !
      errmax = 0.0E0
      !
      DO i = 1, n
        errmax = MAX(errmax,ABS(btemp(i)-bxex(i)))
      END DO
      IF ( errcmp>errmax ) THEN
        IF ( Kprint>=3 ) WRITE (Lun,FMT=99004) list(kprog)
      ELSE
        IF ( Kprint>=2 ) WRITE (Lun,FMT=99005) list(kprog), errmax
        Nerr = Nerr + 1
      END IF
      !
      !       Second test case - singular matrix
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
      DO j = 1, n
        atemp(1,j) = 0.0E0
      END DO
      IF ( kprog==1 ) THEN
        CALL SGEFS(atemp,lda,n,btemp,itask,ind,work,iwork)
      ELSE
        CALL SGEIR(atemp,lda,n,btemp,itask,ind,work,iwork)
      END IF
      IF ( ind==-4 ) THEN
        IF ( Kprint>=3 ) WRITE (Lun,FMT=99006) list(kprog)
      ELSE
        IF ( Kprint>=2 ) WRITE (Lun,FMT=99007) list(kprog), ind
        Nerr = Nerr + 1
      END IF
      !
    END DO
    !
    IF ( Kprint>=3.AND.Nerr==0 ) WRITE (Lun,99002)
    99002 FORMAT (/,2X,'SGEFS and SGEIR Quick Check PASSED'/)
    IF ( Kprint>=2.AND.Nerr/=0 ) WRITE (Lun,99003)
    99003 FORMAT (/,2X,'SGEFS and SGEIR Quick Check FAILED'/)
    RETURN
    99004 FORMAT (/,5X,'S',A,' Normal test PASSED')
    99005 FORMAT (/,5X,'S',A,' Test FAILED, MAX ABS(ERROR) is',E13.5)
    99006 FORMAT (/,5X,'S',A,' Singular test PASSED')
    99007 FORMAT (/,5X,'S',A,' Singular test FAILED, IND=',I3)
  END SUBROUTINE SGEQC
  !** DGEQC
  SUBROUTINE DGEQC(Lun,Kprint,Nerr)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for DGEFS.
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

    !     .. Scalar Arguments ..
    INTEGER Kprint, Lun, Nerr
    !     .. Local Scalars ..
    REAL(8) :: errcmp, errmax
    INTEGER i, ind, itask, j, kprog, lda, n
    !     .. Local Arrays ..
    REAL(8) :: atemp(5,4), btemp(4), work(20)
    INTEGER iwork(4)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX
    !     .. Data statements ..
    REAL(8), PARAMETER :: a(5,4) = RESHAPE( [ 5.0D0, 1.0D0, 0.3D0, 2.1D0, 0.0D0, &
      -1.0D0, -0.5D0, 1.0D0, 1.0D0, 0.0D0, &
      4.5D0,  -1.0D0, -1.7D0, 2.0D0, 0.0D0, &
      0.5D0, 2.0D0, 0.6D0, 1.3D0, 0.0D0 ], [5,4] )
    REAL(8), PARAMETER :: b(4) = [ 0.0D0, 3.5D0, 3.6D0, 2.4D0 ]
    REAL(8), PARAMETER :: bxex(4) = [ 0.10D+01, 0.10D+01, -0.10D+01, 0.10D+01 ]
    CHARACTER(4), PARAMETER :: list(2) = [ 'GEFS', 'GEIR' ]
    !* FIRST EXECUTABLE STATEMENT  DGEQC
    n = 4
    lda = 5
    Nerr = 0
    errcmp = D1MACH(4)**0.8D0
    IF ( Kprint>=2 ) WRITE (Lun,99001)
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
    IF ( ind<0 ) THEN
      IF ( Kprint>=2 ) WRITE (Lun,FMT=99005) list(kprog), ind
      Nerr = Nerr + 1
    END IF
    !
    !     Calculate error for first test
    !
    errmax = 0.0D0
    !
    DO i = 1, n
      errmax = MAX(errmax,ABS(btemp(i)-bxex(i)))
    END DO
    IF ( errcmp>errmax ) THEN
      IF ( Kprint>=3 ) WRITE (Lun,FMT=99004) list(kprog)
    ELSE
      IF ( Kprint>=2 ) WRITE (Lun,FMT=99005) list(kprog), errmax
      Nerr = Nerr + 1
    END IF
    !
    !     Second test case - singular matrix
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
    DO j = 1, n
      atemp(1,j) = 0.0D0
    END DO
    CALL DGEFS(atemp,lda,n,btemp,itask,ind,work,iwork)
    IF ( ind==-4 ) THEN
      IF ( Kprint>=3 ) WRITE (Lun,FMT=99006) list(kprog)
    ELSE
      IF ( Kprint>=2 ) WRITE (Lun,FMT=99007) list(kprog), ind
      Nerr = Nerr + 1
    END IF
    !
    IF ( Kprint>=3.AND.Nerr==0 ) WRITE (Lun,99002)
    99002 FORMAT (/,2X,'DGEFS Quick Check PASSED'/)
    IF ( Kprint>=2.AND.Nerr/=0 ) WRITE (Lun,99003)
    99003 FORMAT (/,2X,'SGEFS and SGEIR Quick Check FAILED'/)
    RETURN
    99004 FORMAT (/,5X,'D',A,' Normal test PASSED')
    99005 FORMAT (/,5X,'D',A,' Test FAILED, MAX ABS(ERROR) is',E13.5)
    99006 FORMAT (/,5X,'D',A,' Singular test PASSED')
    99007 FORMAT (/,5X,'D',A,' Singular test FAILED, IND=',I3)
  END SUBROUTINE DGEQC
  !** CGEQC
  SUBROUTINE CGEQC(Lun,Kprint,Nerr)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for CGEFS and CGEIR.
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

    !     .. Scalar Arguments ..
    INTEGER Kprint, Lun, Nerr
    !     .. Local Scalars ..
    INTEGER i, ind, indx, itask, j, kprog, lda, n
    !     .. Local Arrays ..
    COMPLEX atemp(5,3), btemp(3), work(12)
    INTEGER iwork(3)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, AIMAG, REAL
    !     .. Data statements ..
    COMPLEX, PARAMETER :: a(3,3) = RESHAPE( [ (2.,3.), (1.,1.), (1.,2.), &
      (2.,0.), (1.,-1.), (0.,0.),   (0.,0.), (2.,5.), (3.,2.) ], [3,3] )
    COMPLEX, PARAMETER :: b(3) = [ (-1.,1.), (-5.,4.), (-4.,7.) ]
    COMPLEX, PARAMETER :: bxex(3) = [ (.21459E-01,.209012E+01), &
      (.261373E+01,-.162231E+01), (.785407E+00,.109871E+01) ]
    CHARACTER(4), PARAMETER :: list(2) = [ 'GEFS', 'GEIR' ]
    !* FIRST EXECUTABLE STATEMENT  CGEQC
    n = 3
    lda = 5
    Nerr = 0
    IF ( Kprint>=2 ) WRITE (Lun,99001)
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
      IF ( kprog==1 ) THEN
        CALL CGEFS(atemp,lda,n,btemp,itask,ind,work,iwork)
      ELSE
        CALL CGEIR(atemp,lda,n,btemp,itask,ind,work,iwork)
      END IF
      IF ( ind<0 ) THEN
        IF ( Kprint>=2 ) WRITE (Lun,FMT=99005) list(kprog), ind
        Nerr = Nerr + 1
      END IF
      !
      !       Calculate error for first test
      !
      indx = 0
      DO i = 1, n
        IF ( CABS1(bxex(i)-btemp(i))>.0001 ) indx = indx + 1
      END DO
      IF ( indx==0 ) THEN
        IF ( Kprint>=3 ) WRITE (Lun,FMT=99004) list(kprog)
      ELSE
        IF ( Kprint>=2 ) WRITE (Lun,FMT=99005) list(kprog)
        Nerr = Nerr + 1
      END IF
      !
      !       Second test case - singular matrix
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
      DO j = 1, n
        atemp(1,j) = (0.E0,0.E0)
      END DO
      IF ( kprog==1 ) THEN
        CALL CGEFS(atemp,lda,n,btemp,itask,ind,work,iwork)
      ELSE
        CALL CGEIR(atemp,lda,n,btemp,itask,ind,work,iwork)
      END IF
      IF ( ind==-4 ) THEN
        IF ( Kprint>=3 ) WRITE (Lun,FMT=99006) list(kprog)
      ELSE
        IF ( Kprint>=2 ) WRITE (Lun,FMT=99007) list(kprog), ind
        Nerr = Nerr + 1
      END IF
    END DO
    !
    IF ( Kprint>=3.AND.Nerr==0 ) WRITE (Lun,99002)
    99002 FORMAT (/,2X,'CGEFS and CGEIR Quick Check PASSED'/)
    IF ( Kprint>=2.AND.Nerr/=0 ) WRITE (Lun,99003)
    99003 FORMAT (/,2X,'CGEFS and CGEIR Quick Check FAILED'/)
    RETURN
    99004 FORMAT (/,5X,'C',A,' Normal test PASSED')
    99005 FORMAT (/,5X,'C',A,' Test FAILED')
    99006 FORMAT (/,5X,'C',A,' Singular test PASSED')
    99007 FORMAT (/,5X,'C',A,' Singular test FAILED, IND=',I3)
  END SUBROUTINE CGEQC
END MODULE TEST21_MOD
!** TEST21
PROGRAM TEST21
  USE TEST21_MOD
  use slatec
  IMPLICIT NONE
  !>
  !***
  !  Driver for testing SLATEC subprograms
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
  ! **Routines called:**  CGEQC, DGEQC, I1MACH, SGEQC, XERMAX, XSETF, XSETUN

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  INTEGER kprint, lin, lun, nerr, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST21
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
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST21 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST21 *************')
  END IF
  STOP
END PROGRAM TEST21
