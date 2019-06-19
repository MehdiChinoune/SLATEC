MODULE TEST23_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  !** CGBQC
  SUBROUTINE CGBQC(Lun,Kprint,Nerr)
    !> Quick check for CGBFA, CGBCO, CGBSL and CGBDI.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Voorhees, E. A., (LANL)
    !***
    ! **Description:**
    !
    !    LET  A*X=B  BE A COMPLEX LINEAR SYSTEM WHERE THE MATRIX  A  IS
    !    OF THE PROPER TYPE FOR THE LINPACK SUBROUTINES BEING TESTED.
    !    THE VALUES OF  A  AND  B  AND THE PRE-COMPUTED VALUES OF  C
    !    (THE SOLUTION VECTOR),  DC  (DETERMINANT OF  A ), AND
    !    RCND  (RCOND) ARE ENTERED WITH DATA STATEMENTS.
    !
    !    THE COMPUTED TEST RESULTS FOR  X,  RCOND  AND THE DETER-
    !    MINANT ARE COMPARED TO THE STORED PRE-COMPUTED VALUES.
    !    FAILURE OF THE TEST OCCURS WHEN AGREEMENT TO 3 SIGNIFICANT
    !    DIGITS IS NOT ACHIEVED AND AN ERROR MESSAGE INDICATING WHICH
    !    LINPACK SUBROUTINE FAILED AND WHICH QUANTITY WAS INVOLVED IS
    !    PRINTED.  A SUMMARY LINE IS ALWAYS PRINTED.
    !
    !    NO INPUT ARGUMENTS ARE REQUIRED.
    !    ON RETURN,  NERR  (INTEGER TYPE) CONTAINS THE TOTAL COUNT OF
    !    ALL FAILURES DETECTED BY CGBQC.
    !
    !***
    ! **Routines called:**  CGBCO, CGBDI, CGBFA, CGBSL

    !* REVISION HISTORY  (YYMMDD)
    !   801015  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Restructured using IF-THEN-ELSE-ENDIF, moved an ARITHMETIC
    !           STATEMENT FUNCTION ahead of the FIRST EXECUTABLE STATEMENT
    !           record and cleaned up FORMATs.  (RWC)
    USE linpack, ONLY : CGBFA
    USE blas, ONLY : SCABS1
    INTEGER :: Kprint, Lun
    COMPLEX(SP) :: at(7,4), bt(4)
    INTEGER :: lda, n, ipvt(4), info, i, j, Nerr
    INTEGER :: ml, mu
    COMPLEX(SP), PARAMETER :: abd(6,4) = RESHAPE( [ &
      (0.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), (2.E0,0.E0), (0.E0,1.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), (0.E0,-1.E0), (2.E0,0.E0), (0.E0,0.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), (3.E0,0.E0), (0.E0,1.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), (0.E0,-1.E0), (4.E0,0.E0), (0.E0,0.E0) ], &
      [6,4] )
    COMPLEX(SP), PARAMETER :: b(4) = [ (3.E0,2.E0), (-1.E0,3.E0), (0.E0,-4.E0), &
      (5.E0,0.E0) ]
    CHARACTER(19), PARAMETER :: kprog = 'GBFA GBCO GBSL GBDI'
    CHARACTER(39), PARAMETER :: kfail = 'INFO RCOND SOLUTION DETERMINANT INVERSE'
    !* FIRST EXECUTABLE STATEMENT  CGBQC
    lda = 7
    n = 4
    ml = 1
    mu = 3
    Nerr = 0
    !
    !     FORM AT FOR CGBFA AND BT FOR CGBSL, TEST CGBFA
    !
    DO j = 1, n
      bt(j) = b(j)
      DO i = 1, 6
        at(i,j) = abd(i,j)
      END DO
    END DO
    !
    CALL CGBFA(at,lda,n,ml,mu,ipvt,info)
    IF( info/=0 ) THEN
      WRITE (Lun,99002) kprog(1:4), kfail(1:4)
      Nerr = Nerr + 1
    END IF
    !
    IF( Kprint>=2 .OR. Nerr/=0 ) WRITE (Lun,99001) Nerr
    !
    99001 FORMAT (/' * CGBQC - TEST FOR CGBFA FOUND ',I1,' ERRORS.'/)
    RETURN
    99002 FORMAT (/' *** C',A,' FAILURE - ERROR IN ',A)
  END SUBROUTINE CGBQC
  !** CGECK
  SUBROUTINE CGECK(Lun,Kprint,Nerr)
    !> Quick check for CGEFA, CGECO, CGESL and CGEDI.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Voorhees, E. A., (LANL)
    !***
    ! **Description:**
    !
    !    LET  A*X=B  BE A COMPLEX LINEAR SYSTEM WHERE THE MATRIX  A  IS
    !    OF THE PROPER TYPE FOR THE LINPACK SUBROUTINES BEING TESTED.
    !    THE VALUES OF  A  AND  B  AND THE PRE-COMPUTED VALUES OF  C
    !    (THE SOLUTION VECTOR),  AINV  (INVERSE OF MATRIX  A ),  DC
    !    (DETERMINANT OF  A ), AND  RCND  ( RCOND ) ARE ENTERED
    !    WITH DATA STATEMENTS.
    !
    !    THE COMPUTED TEST RESULTS FOR  X, RCOND, THE DETERMINANT, AND
    !    THE INVERSE ARE COMPARED TO THE STORED PRE-COMPUTED VALUES.
    !    FAILURE OF THE TEST OCCURS WHEN AGREEMENT TO 3 SIGNIFICANT
    !    DIGITS IS NOT ACHIEVED AND AN ERROR MESSAGE INDICATING WHICH
    !    LINPACK SUBROUTINE FAILED AND WHICH QUANTITY WAS INVOLVED IS
    !    PRINTED.  A SUMMARY LINE IS ALWAYS PRINTED.
    !
    !    NO INPUT ARGUMENTS ARE REQUIRED.
    !    ON RETURN,  NERR  (INTEGER TYPE) CONTAINS THE TOTAL COUNT OF
    !    ALL FAILURES DETECTED BY CGECK.
    !
    !***
    ! **Routines called:**  CGECO, CGEDI, CGEFA, CGESL

    !* REVISION HISTORY  (YYMMDD)
    !   801014  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Restructured using IF-THEN-ELSE-ENDIF and cleaned up
    !           FORMATs.  (RWC)
    USE linpack, ONLY : CGECO, CGEFA, CGESL
    USE blas, ONLY : SCABS1
    INTEGER :: Kprint, Lun
    COMPLEX(SP) :: at(5,4), bt(4), z(4)
    REAL(SP) :: r, rcond
    INTEGER :: lda, n, ipvt(4), info, i, j, indx, Nerr
    COMPLEX(SP), PARAMETER :: a(4,4) = RESHAPE( [ &
      (2.E0,0.E0), (0.E0,1.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,-1.E0), (2.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (3.E0,0.E0), (0.E0,1.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (0.E0,-1.E0), (4.E0,0.E0) ], [4,4] )
    COMPLEX(SP), PARAMETER :: b(4) = [ (3.E0,2.E0), (-1.E0,3.E0), (0.E0,-4.E0), &
      (5.E0,0.E0) ]
    COMPLEX(SP), PARAMETER :: c(4) = [ (1.E0,1.E0), (0.E0,1.E0), (0.E0,-1.E0), &
      (1.E0,0.E0) ]
    CHARACTER(19), PARAMETER :: kprog = 'GEFA GECO GESL GEDI'
    CHARACTER(39), PARAMETER :: kfail = 'INFO RCOND SOLUTION DETERMINANT INVERSE'
    REAL, PARAMETER :: rcnd = .24099E0
    !* FIRST EXECUTABLE STATEMENT  CGECK
    lda = 5
    n = 4
    Nerr = 0
    !
    !     FORM AT FOR CGEFA AND BT FOR CGESL, TEST CGEFA
    !
    DO j = 1, n
      bt(j) = b(j)
      DO i = 1, n
        at(i,j) = a(i,j)
      END DO
    END DO
    !
    CALL CGEFA(at,lda,n,ipvt,info)
    IF( info/=0 ) THEN
      WRITE (Lun,99002) kprog(1:4), kfail(1:4)
      Nerr = Nerr + 1
    END IF
    !
    !     TEST CGESL FOR JOB=0
    !
    CALL CGESL(at,lda,n,ipvt,bt,0)
    indx = 0
    DO i = 1, n
      IF( SCABS1(c(i)-bt(i))>.0001 ) indx = indx + 1
    END DO
    !
    IF( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(11:14), kfail(12:19)
      Nerr = Nerr + 1
    END IF
    !
    !     FORM AT FOR CGECO AND BT FOR CGESL, TEST CGECO
    !
    DO j = 1, n
      bt(j) = b(j)
      DO i = 1, n
        at(i,j) = a(i,j)
      END DO
    END DO
    !
    CALL CGECO(at,lda,n,ipvt,rcond,z)
    r = ABS(rcnd-rcond)
    IF( r>=.0001 ) THEN
      WRITE (Lun,99002) kprog(6:9), kfail(6:10)
      Nerr = Nerr + 1
    END IF
    !
    !     TEST CGESL FOR JOB NOT EQUAL TO 0
    !
    CALL CGESL(at,lda,n,ipvt,bt,1)
    indx = 0
    DO i = 1, n
      IF( SCABS1(c(i)-bt(i))>.0001 ) indx = indx + 1
    END DO
    !
    IF( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(11:14), kfail(12:19)
      Nerr = Nerr + 1
    END IF
    !
    IF( Kprint>=2 .OR. Nerr/=0 ) WRITE (Lun,99001) Nerr
    !
    99001 FORMAT (/' * CGECK - TEST FOR CGEFA, CGECO AND CGESL FOUND ',I1,&
      ' ERRORS.'/)
    RETURN
    99002 FORMAT (/' *** C',A,' FAILURE - ERROR IN ',A)
  END SUBROUTINE CGECK
  !** CPOQC
  SUBROUTINE CPOQC(Lun,Kprint,Nerr)
    !> Quick check for CPOFA, CPOCO, CPOSL and CPODI.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Voorhees, E. A., (LANL)
    !***
    ! **Description:**
    !
    !    LET  A*X=B  BE A COMPLEX LINEAR SYSTEM WHERE THE MATRIX  A  IS
    !    OF THE PROPER TYPE FOR THE LINPACK SUBROUTINES BEING TESTED.
    !    THE VALUES OF  A  AND  B  AND THE PRE-COMPUTED VALUES OF  C
    !    (THE SOLUTION VECTOR),  AINV  (INVERSE OF MATRIX  A ),  DC
    !    (DETERMINANT OF  A ), AND  RCND  ( RCOND ) ARE ENTERED
    !    WITH DATA STATEMENTS.
    !
    !    THE COMPUTED TEST RESULTS FOR  X, RCOND, THE DETERMINANT, AND
    !    THE INVERSE ARE COMPARED TO THE STORED PRE-COMPUTED VALUES.
    !    FAILURE OF THE TEST OCCURS WHEN AGREEMENT TO 3 SIGNIFICANT
    !    DIGITS IS NOT ACHIEVED AND AN ERROR MESSAGE INDICATING WHICH
    !    LINPACK SUBROUTINE FAILED AND WHICH QUANTITY WAS INVOLVED IS
    !    PRINTED.  A SUMMARY LINE IS ALWAYS PRINTED.
    !
    !    NO INPUT ARGUMENTS ARE REQUIRED.
    !    ON RETURN,  NERR  (INTEGER TYPE) CONTAINS THE TOTAL COUNT OF
    !    ALL FAILURES DETECTED BY CPOQC.
    !
    !***
    ! **Routines called:**  CPOCO, CPODI, CPOFA, CPOSL

    !* REVISION HISTORY  (YYMMDD)
    !   801016  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Restructured using IF-THEN-ELSE-ENDIF and cleaned up
    !           FORMATs.  (RWC)
    USE linpack, ONLY : CPOCO, CPOFA, CPOSL
    USE blas, ONLY : SCABS1
    INTEGER :: Kprint, Lun
    COMPLEX(SP) :: at(5,4), bt(4), z(4)
    REAL(SP) :: r, rcond
    INTEGER :: lda, n, info, i, j, indx, Nerr
    COMPLEX(SP), PARAMETER :: a(4,4) = RESHAPE( [ &
      (2.E0,0.E0), (0.E0,1.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,-1.E0), (2.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (3.E0,0.E0), (0.E0,1.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (0.E0,-1.E0), (4.E0,0.E0) ], [4,4] )
    COMPLEX(SP), PARAMETER :: b(4) = [ (3.E0,2.E0), (-1.E0,3.E0), (0.E0,-4.E0), &
      (5.E0,0.E0) ]
    COMPLEX(SP), PARAMETER :: c(4) = [ (1.E0,1.E0), (0.E0,1.E0), (0.E0,-1.E0), &
      (1.E0,0.E0) ]
    CHARACTER(19), PARAMETER :: kprog = 'POFA POCO POSL PODI'
    CHARACTER(39), PARAMETER :: kfail = 'INFO RCOND SOLUTION DETERMINANT INVERSE'
    REAL, PARAMETER :: rcnd = .24099E0
    !* FIRST EXECUTABLE STATEMENT  CPOQC
    lda = 5
    n = 4
    Nerr = 0
    !
    !     FORM AT FOR CPOFA AND BT FOR CPOSL, TEST CPOFA
    !
    DO j = 1, n
      bt(j) = b(j)
      DO i = 1, n
        at(i,j) = a(i,j)
      END DO
    END DO
    !
    CALL CPOFA(at,lda,n,info)
    IF( info/=0 ) THEN
      WRITE (Lun,99002) kprog(1:4), kfail(1:4)
      Nerr = Nerr + 1
    END IF
    !
    !     TEST CPOSL
    !
    CALL CPOSL(at,lda,n,bt)
    indx = 0
    DO i = 1, n
      IF( SCABS1(c(i)-bt(i))>.0001 ) indx = indx + 1
    END DO
    !
    IF( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(11:14), kfail(12:19)
      Nerr = Nerr + 1
    END IF
    !
    !     FORM AT FOR CPOCO, TEST CPOCO
    !
    DO j = 1, n
      DO i = 1, n
        at(i,j) = a(i,j)
      END DO
    END DO
    !
    CALL CPOCO(at,lda,n,rcond,z,info)
    r = ABS(rcnd-rcond)
    IF( r>=.0001 ) THEN
      WRITE (Lun,99002) kprog(6:9), kfail(6:10)
      Nerr = Nerr + 1
    END IF
    !
    IF( info/=0 ) THEN
      WRITE (Lun,99002) kprog(6:9), kfail(1:4)
      Nerr = Nerr + 1
    END IF
    !
    IF( Kprint>=2 .OR. Nerr/=0 ) WRITE (Lun,99001) Nerr
    !
    99001 FORMAT (/' * CPOQC - TEST FOR CPOFA, CPOCO AND CPOSL FOUND ',I1,&
      ' ERRORS.'/)
    RETURN
    99002 FORMAT (/' *** C',A,' FAILURE - ERROR IN ',A)
  END SUBROUTINE CPOQC
END MODULE TEST23_MOD
!** TEST23
PROGRAM TEST23
  USE TEST23_MOD, ONLY : CGBQC, CGECK, CPOQC
  USE slatec, ONLY : I1MACH, XSETF, XSETUN, XERMAX
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !> Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  D2
  !***
  ! **Type:**      COMPLEX (TEST23-S)
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
  !        CGECO    CGEDI    CGEFA    CGESL
  !        CGBCO    CGBDI    CGBFA    CGBSL
  !        CPOCO    CPODI    CPOFA    CPOSL
  !        CPPCO    CPPDI    CPPFA    CPPSL
  !        CPBCO    CPBDI    CPBFA    CPBSL
  !        CSICO    CSIDI    CSIFA    CSISL
  !        CSPCO    CSPDI    CSPFA    CSPSL
  !        CHICO    CHIDI    CHIFA    CHISL
  !        CHPCO    CHPDI    CHPFA    CHPSL
  !        CTRCO    CTRDI      -      CTRSL
  !        CGTSL
  !        CPTSL
  !        CCHDC
  !        CQRDC    CQRSL
  !        CSVDC
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  CCHQC, CGBQC, CGECK, CGTQC, CHIQC, CHPQC, CPBQC,
  !                    CPOQC, CPPQC, CPTQC, CQRQC, CSIQC, CSPQC, CSVQC,
  !                    CTRQC, I1MACH, XERMAX, XSETF, XSETUN

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  INTEGER :: kprint, lin, lun, nerr, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST23
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
  !     Test LINPACK routines
  !
  CALL CGECK(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CGBQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CPOQC(lun,kprint,nerr)
  nfail = nfail + nerr
  !
  !
  !     Write PASS or FAIL message
  !
  IF( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST23 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST23 *************')
  END IF
  STOP
END PROGRAM TEST23
