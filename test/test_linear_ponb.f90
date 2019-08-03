MODULE TEST22_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  !** SQCK
  SUBROUTINE SQCK(Lun,Kprint,Nerr)
    !> Quick check for SPOFS, SPOIR, SNBFS and SNBIR.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Voorhees, E. A., (LANL)
    !***
    ! **Description:**
    !
    !    QUICK CHECK SUBROUTINE SQCK TESTS THE EXECUTION OF THE
    !    SLATEC SUBROUTINES SPOFS, SPOIR, SNBFS AND SNBIR.
    !    A TITLE LINE AND A SUMMARY LINE ARE ALWAYS OUTPUTTED.
    !
    !    THE SUMMARY LINE GIVES A COUNT OF THE NUMBER OF
    !    PROBLEMS ENCOUNTERED IN THE TEST IF ANY EXIST.  SQCK
    !    CHECKS COMPUTED VS. EXACT SOLUTIONS TO AGREE TO
    !    WITHIN 0.8 TIMES THE WORD LENGTH OF THE COMPUTER
    !    (1.6 IF DOUBLE PRECISION) FOR CASE 1.  SQCK ALSO
    !    TESTS ERROR HANDLING BY THE SUBROUTINE (CALLS TO
    !    XERMSG (SQCK SETS IFLAG/KONTRL TO 0))
    !    USING A SINGULAR MATRIX FOR CASE 2.  EACH EXECUTION
    !    PROBLEM DETECTED BY SQCK RESULTS IN AN ADDITIONAL
    !    EXPLANATORY LINE OF OUTPUT.
    !
    !    SQCK REQUIRES NO INPUT ARGUMENTS.
    !    ON RETURN, NERR (INTEGER TYPE) CONTAINS THE TOTAL COUNT
    !    OF ALL PROBLEMS DETECTED BY SQCK.
    !
    !***
    ! **Routines called:**  R1MACH, SNBFS, SNBIR, SPOFS, SPOIR

    !* REVISION HISTORY  (YYMMDD)
    !   800930  DATE WRITTEN
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   891009  Removed unreferenced statement label.  (WRB)
    !   891009  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901009  Routine writes illegal character to column 1, fixed.
    !           Editorial changes made, code fixed to test all four routines.  (RWC)
    !   901009  Restructured using IF-THEN-ELSE-ENDIF, cleaned up FORMATs,
    !           including removing an illegal character from column 1, and
    !           fixed code to test all four routines.  (RWC)
    USE service, ONLY : eps_sp
    USE linear, ONLY : SNBFS, SNBIR, SPOFS, SPOIR
    INTEGER :: Kprint, Lun
    REAL(SP) :: at(5,4), abe(5,7), abet(5,7), b(4), bt(4), c(4), &
      work(35), r, delx, delmax, signn
    INTEGER :: lda, n, ml, mu, ind, iwork(4), Nerr, i, j, j1, j2, jd, &
      mlp, k, kcase, kprog
    REAL(SP), PARAMETER :: a(4,4) = RESHAPE( [ 5._SP, 4._SP, 1._SP, 1._SP, &
      4._SP, 5._SP, 1._SP, 1._SP,    1._SP, 1._SP, 4._SP, 2._SP, &
      1._SP, 1._SP, 2._SP, 4._SP ], [4,4] )
    CHARACTER(4), PARAMETER :: list(4) = [ 'POFS', 'POIR', 'NBFS', 'NBIR' ]
    !* FIRST EXECUTABLE STATEMENT  SQCK
    IF( Kprint>=3 ) WRITE (Lun,99001)
    !
    99001 FORMAT (/' *    SQCK - QUICK CHECK FOR SPOFS, SPOIR, SNBFS AND ','SNBIR'/)
    lda = 5
    n = 4
    ml = 2
    mu = 1
    jd = 2*ml + mu + 1
    Nerr = 0
    r = eps_sp**0.8_SP
    !
    !     COMPUTE C VECTOR.
    !
    signn = 1._SP
    DO i = 1, n
      c(i) = signn/i
      signn = -signn
    END DO
    !
    !     CASE 1 FOR WELL-CONDITIONED MATRIX, CASE 2 FOR SINGULAR MATRIX.
    !
    DO kcase = 1, 1 !2
      DO kprog = 1, 4
        !           SET VECTOR B TO ZERO.
        DO i = 1, n
          b(i) = 0._SP
        END DO
        !
        !           FORM VECTOR B FOR NON-BANDED.
        !
        IF( kprog<=2 ) THEN
          DO i = 1, n
            DO j = 1, n
              b(i) = b(i) + a(i,j)*c(j)
            END DO
          END DO
        ELSE
          !
          !              FORM ABE(NB ARRAY) FROM MATRIX A
          !              AND FORM VECTOR B FOR BANDED.
          !
          DO j = 1, jd
            DO i = 1, n
              abe(i,j) = 0._SP
            END DO
          END DO
          !
          mlp = ml + 1
          DO i = 1, n
            j1 = MAX(1,i-ml)
            j2 = MIN(n,i+mu)
            DO j = j1, j2
              k = j - i + mlp
              abe(i,k) = a(i,j)
              b(i) = b(i) + (a(i,j)*c(j))
            END DO
          END DO
        END IF
        !
        !           FORM BT FROM B, AT FROM A, AND ABET FROM ABE.
        !
        DO i = 1, n
          bt(i) = b(i)
          DO j = 1, n
            at(i,j) = a(i,j)
          END DO
        END DO
        !
        DO j = 1, jd
          DO i = 1, n
            abet(i,j) = abe(i,j)
          END DO
        END DO
        !
        !           MAKE AT AND ABET SINGULAR FOR CASE  =  2
        !
        IF( kcase==2 ) THEN
          DO j = 1, n
            at(1,j) = 0._SP
          END DO
          !
          DO j = 1, jd
            abet(1,j) = 0._SP
          END DO
        END IF
        !
        !           SOLVE FOR X
        !
        IF( kprog==1 ) CALL SPOFS(at,lda,n,bt,1,ind,work)
        IF( kprog==2 ) CALL SPOIR(at,lda,n,bt,1,ind,work)
        IF( kprog==3 ) CALL SNBFS(abet,lda,n,ml,mu,bt,1,ind,work,iwork)
        IF( kprog==4 ) CALL SNBIR(abet,lda,n,ml,mu,bt,1,ind,work,iwork)
        !
        !           COMPARE EXACT AND COMPUTED SOLUTIONS FOR CASE 1
        !
        IF( kcase==1 ) THEN
          delmax = 0._SP
          DO i = 1, n
            delx = ABS(bt(i)-c(i))
            delmax = MAX(delmax,delx)
          END DO
          !
          IF( r<=delmax ) THEN
            Nerr = Nerr + 1
            WRITE (Lun,99002) list(kprog), kcase, delmax
            99002 FORMAT ('   PROBLEM WITH S',A,', CASE ',I1,'.  MAX ABS ERROR OF',&
              E11.4/)
          END IF
          !
          !              CHECK CONTROL FOR SINGULAR MATRIX FOR CASE 2
          !
        ELSEIF( ind/=-4 ) THEN
          Nerr = Nerr + 1
          WRITE (Lun,99003) list(kprog), kcase, ind
          99003 FORMAT ('   PROBLEM WITH S',A,', CASE ',I1,'.  IND = ',I2,&
            ' INSTEAD OF -4'/)
        END IF
      END DO
    END DO
    !
    !     SUMMARY PRINT
    !
    IF( Nerr/=0 ) WRITE (Lun,99004) Nerr
    99004 FORMAT (/' **** SQCK DETECTED A TOTAL OF ',I2,' PROBLEMS. ****'/)
    IF( Kprint>=2 .AND. Nerr==0 ) WRITE (Lun,99005)
    99005 FORMAT ('     SQCK DETECTED NO PROBLEMS.'/)
    RETURN
  END SUBROUTINE SQCK
  !** DQCK
  SUBROUTINE DQCK(Lun,Kprint,Nerr)
    !> Quick check for DPOFS AND DNBFS.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Voorhees, E. A., (LANL)
    !***
    ! **Description:**
    !
    !    QUICK CHECK SUBROUTINE DQCK TESTS THE EXECUTION OF THE
    !    SLATEC SUBROUTINES DPOFS AND DNBFS.
    !    A TITLE LINE AND A SUMMARY LINE ARE ALWAYS OUTPUTTED.
    !
    !    THE SUMMARY LINE GIVES A COUNT OF THE NUMBER OF
    !    PROBLEMS ENCOUNTERED IN THE TEST IF ANY EXIST.  DQCK
    !    CHECKS COMPUTED VS. EXACT SOLUTIONS TO AGREE TO
    !    WITHIN 0.8 TIMES THE WORD LENGTH OF THE COMPUTER
    !    (1.6 IF DOUBLE PRECISION) FOR CASE 1.  DQCK ALSO
    !    TESTS ERROR HANDLING BY THE SUBROUTINE (CALLS TO
    !    XERMSG (DQCK SETS IFLAG/KONTRL TO 0))
    !    USING A SINGULAR MATRIX FOR CASE 2.  EACH EXECUTION
    !    PROBLEM DETECTED BY DQCK RESULTS IN AN ADDITIONAL
    !    EXPLANATORY LINE OF OUTPUT.
    !
    !    DQCK REQUIRES NO INPUT ARGUMENTS.
    !    ON RETURN, NERR (INTEGER TYPE) CONTAINS THE TOTAL COUNT
    !    OF ALL PROBLEMS DETECTED BY DQCK.
    !
    !***
    ! **Routines called:**  D1MACH, DNBFS, DPOFS

    !* REVISION HISTORY  (YYMMDD)
    !   801002  DATE WRITTEN
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   890911  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901009  Restructured using IF-THEN-ELSE-ENDIF, cleaned up FORMATs,
    !           including removing an illegal character from column 1, and
    !           editorial changes.  (RWC)
    USE service, ONLY : eps_dp
    USE linear, ONLY : DNBFS, DPOFS
    !
    INTEGER :: Kprint, Lun
    REAL(DP) :: at(5,4), abe(5,7), abet(5,7), b(4), bt(4), c(4), work(35), signn
    REAL(SP) :: r, delx, delmax
    INTEGER :: lda, n, ml, mu, ind, iwork(4), Nerr, i, j, j1, j2, jd, &
      mlp, k, kcase, kprog
    REAL(DP), PARAMETER :: a(4,4) = RESHAPE( [ 5._DP, 4._DP, 1._DP, 1._DP, &
      4._DP, 5._DP, 1._DP, 1._DP,    1._DP, 1._DP, 4._DP, 2._DP, &
      1._DP, 1._DP, 2._DP, 4._DP ], [4,4] )
    CHARACTER(4), PARAMETER :: list(2) = [ 'POFS', 'NBFS' ]
    !* FIRST EXECUTABLE STATEMENT  DQCK
    IF( Kprint>=3 ) WRITE (Lun,99001)
    !
    99001 FORMAT (/' *    DQCK - QUICK CHECK FOR  DPOFS AND DNBFS'/)
    lda = 5
    n = 4
    ml = 2
    mu = 1
    jd = 2*ml + mu + 1
    Nerr = 0
    r = REAL( eps_dp**0.8_DP, SP )
    !
    !     COMPUTE C VECTOR.
    !
    signn = 1._DP
    DO i = 1, n
      c(i) = signn/i
      signn = -signn
    END DO
    !
    !     CASE 1 FOR WELL-CONDITIONED MATRIX, CASE 2 FOR SINGULAR MATRIX.
    !
    DO kcase = 1, 1 !2
      DO kprog = 1, 2
        !           SET VECTOR B TO ZERO.
        DO i = 1, n
          b(i) = 0._DP
        END DO
        !
        !           FORM VECTOR B FOR NON-BANDED.
        !
        IF( kprog==1 ) THEN
          DO i = 1, n
            DO j = 1, n
              b(i) = b(i) + a(i,j)*c(j)
            END DO
          END DO
        ELSE
          !
          !              FORM ABE(NB ARRAY) FROM MATRIX A
          !              AND FORM VECTOR B FOR BANDED.
          !
          DO j = 1, jd
            DO i = 1, n
              abe(i,j) = 0._DP
            END DO
          END DO
          !
          mlp = ml + 1
          DO i = 1, n
            j1 = MAX(1,i-ml)
            j2 = MIN(n,i+mu)
            DO j = j1, j2
              k = j - i + mlp
              abe(i,k) = a(i,j)
              b(i) = b(i) + (a(i,j)*c(j))
            END DO
          END DO
        END IF
        !
        !           FORM BT FROM B, AT FROM A, AND ABET FROM ABE.
        !
        DO i = 1, n
          bt(i) = b(i)
          DO j = 1, n
            at(i,j) = a(i,j)
          END DO
        END DO
        !
        DO j = 1, jd
          DO i = 1, n
            abet(i,j) = abe(i,j)
          END DO
        END DO
        !
        !           MAKE AT AND ABET SINGULAR FOR CASE  =  2
        !
        IF( kcase==2 ) THEN
          DO j = 1, n
            at(1,j) = 0._DP
          END DO
          !
          DO j = 1, jd
            abet(1,j) = 0._DP
          END DO
        END IF
        !
        !           SOLVE FOR X
        !
        IF( kprog==1 ) CALL DPOFS(at,lda,n,bt,1,ind,work)
        IF( kprog==2 ) CALL DNBFS(abet,lda,n,ml,mu,bt,1,ind,work,iwork)
        !
        !           COMPARE EXACT AND COMPUTED SOLUTIONS FOR CASE 1
        !
        IF( kcase==1 ) THEN
          delmax = 0._SP
          DO i = 1, n
            delx = REAL( ABS(bt(i)-c(i)), SP )
            delmax = MAX(delmax,delx)
          END DO
          !
          IF( r<=delmax ) THEN
            Nerr = Nerr + 1
            WRITE (Lun,99002) list(kprog), kcase, delmax
            99002 FORMAT ('   PROBLEM WITH D',A,', CASE ',I1,'.  MAX ABS ERROR OF',&
              E11.4/)
          END IF
          !
          !              CHECK CONTROL FOR SINGULAR MATRIX FOR CASE 2
          !
        ELSEIF( ind/=-4 ) THEN
          Nerr = Nerr + 1
          WRITE (Lun,99003) list(kprog), kcase, ind
          99003 FORMAT ('   PROBLEM WITH D',A,', CASE ',I1,'.  IND = ',I2,&
            ' INSTEAD OF -4'/)
        END IF
      END DO
    END DO
    !
    !     SUMMARY PRINT
    !
    IF( Nerr/=0 ) WRITE (Lun,99004) Nerr
    99004 FORMAT (/' **** DQCK DETECTED A TOTAL OF ',I2,' PROBLEMS. ****'/)
    IF( Kprint>=2 .AND. Nerr==0 ) WRITE (Lun,99005)
    99005 FORMAT ('     DQCK DETECTED NO PROBLEMS.'/)
    RETURN
  END SUBROUTINE DQCK
  !** CQCK
  SUBROUTINE CQCK(Lun,Kprint,Nerr)
    !> Quick check for CPOFS, CPOIR, CNBFS and CNBIR.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Voorhees, E. A., (LANL)
    !***
    ! **Description:**
    !
    !    QUICK CHECK SUBROUTINE CQCK TESTS THE EXECUTION OF THE
    !    SLATEC SUBROUTINES CPOFS, CPOIR, CNBFS AND CNBIR.
    !    A TITLE LINE AND A SUMMARY LINE ARE ALWAYS OUTPUTTED.
    !
    !    THE SUMMARY LINE GIVES A COUNT OF THE NUMBER OF
    !    PROBLEMS ENCOUNTERED IN THE TEST IF ANY EXIST.  CQCK
    !    CHECKS COMPUTED VS. EXACT SOLUTIONS TO AGREE TO
    !    WITHIN 0.8 TIMES THE WORD LENGTH OF THE COMPUTER
    !    (1.6 IF DOUBLE PRECISION) FOR CASE 1.  CQCK ALSO
    !    TESTS ERROR HANDLING BY THE SUBROUTINE (CALLS TO
    !    XERMSG (CQCK SETS IFLAG/KONTRL TO 0))
    !    USING A SINGULAR MATRIX FOR CASE 2.  EACH EXECUTION
    !    PROBLEM DETECTED BY CQCK RESULTS IN AN ADDITIONAL
    !    EXPLANATORY LINE OF OUTPUT.
    !
    !    CQCK REQUIRES NO INPUT ARGUMENTS.
    !    ON RETURN, NERR (INTEGER TYPE) CONTAINS THE TOTAL COUNT
    !    OF ALL PROBLEMS DETECTED BY CQCK.
    !
    !***
    ! **Routines called:**  CNBFS, CNBIR, CPOFS, CPOIR, R1MACH

    !* REVISION HISTORY  (YYMMDD)
    !   801002  DATE WRITTEN
    !   891009  Removed unreferenced statement labels.  (WRB)
    !   891009  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901009  Restructured using IF-THEN-ELSE-ENDIF, cleaned up FORMATs,
    !           including removing an illegal character from column 1, and
    !           editorial changes.  (RWC)
    USE service, ONLY : eps_sp
    USE linear, ONLY : CNBFS, CNBIR, CPOFS, CPOIR
    !
    INTEGER :: Kprint, Lun
    REAL(SP) :: r, delx, delmax
    COMPLEX(SP) :: at(5,4), abe(5,7), abet(5,7), bt(4), work(35)
    INTEGER :: lda, n, ml, mu, ind, iwork(4), Nerr, i, j, j1, j2, jd, &
      mlp, k, kcase, kprog
    COMPLEX(SP), PARAMETER :: a(4,4) = RESHAPE( [ &
      (2._SP,0._SP), (0._SP,1._SP), (0._SP,0._SP), (0._SP,0._SP), &
      (0._SP,-1._SP), (2._SP,0._SP), (0._SP,0._SP), (0._SP,0._SP), &
      (0._SP,0._SP), (0._SP,0._SP), (3._SP,0._SP), (0._SP,1._SP), &
      (0._SP,0._SP), (0._SP,0._SP), (0._SP,-1._SP), (4._SP,0._SP) ], [4,4] )
    COMPLEX(SP), PARAMETER :: c(4) = [ (1._SP,1._SP), (0._SP,1._SP), (0._SP,-1._SP), &
      (1._SP,0._SP) ]
    COMPLEX(SP), PARAMETER :: b(4) = [ (3._SP,2._SP), (-1._SP,3._SP), (0._SP,-4._SP), &
      (5._SP,0._SP) ]
    CHARACTER(4), PARAMETER :: list(4) = [ 'POFS', 'POIR', 'NBFS', 'NBIR' ]
    !* FIRST EXECUTABLE STATEMENT  CQCK
    IF( Kprint>=3 ) WRITE (Lun,99001)
    !
    99001 FORMAT (/' *    CQCK - QUICK CHECK FOR CPOFS, CPOIR, CNBFS AND ','CNBIR'/)
    lda = 5
    n = 4
    ml = 2
    mu = 1
    jd = 2*ml + mu + 1
    Nerr = 0
    r = eps_sp**0.8_SP
    !
    !     FORM ABE(NB ARRAY) FROM MATRIX A.
    !
    DO j = 1, jd
      DO i = 1, n
        abe(i,j) = (0._SP,0._SP)
      END DO
    END DO
    !
    mlp = ml + 1
    DO i = 1, n
      j1 = MAX(1,i-ml)
      j2 = MIN(n,i+mu)
      DO j = j1, j2
        k = j - i + mlp
        abe(i,k) = a(i,j)
      END DO
    END DO
    !
    !     CASE 1 FOR WELL-CONDITIONED MATRIX, CASE 2 FOR SINGULAR MATRIX
    !
    DO kcase = 1, 1 !2
      DO kprog = 1, 4
        !           FORM BT FROM B, AT FROM A, AND ABET FROM ABE.
        DO i = 1, n
          bt(i) = b(i)
          DO j = 1, n
            at(i,j) = a(i,j)
          END DO
        END DO
        !
        DO j = 1, jd
          DO i = 1, n
            abet(i,j) = abe(i,j)
          END DO
        END DO
        !
        !           MAKE AT AND ABET SINGULAR FOR CASE  =  2
        !
        IF( kcase==2 ) THEN
          DO j = 1, n
            at(1,j) = (0._SP,0._SP)
          END DO
          !
          DO j = 1, jd
            abet(1,j) = (0._SP,0._SP)
          END DO
        END IF
        !
        !           SOLVE FOR X
        !
        IF( kprog==1 ) CALL CPOFS(at,lda,n,bt,1,ind,work)
        IF( kprog==2 ) CALL CPOIR(at,lda,n,bt,1,ind,work)
        IF( kprog==3 ) CALL CNBFS(abet,lda,n,ml,mu,bt,1,ind,work,iwork)
        IF( kprog==4 ) CALL CNBIR(abet,lda,n,ml,mu,bt,1,ind,work,iwork)
        !
        !           COMPARE EXACT AND COMPUTED SOLUTIONS FOR CASE 1
        !
        IF( kcase==1 ) THEN
          delmax = 0._SP
          DO i = 1, n
            delx = ABS(REAL(bt(i))-REAL(c(i)))
            delmax = MAX(delmax,delx)
            delx = ABS(AIMAG(bt(i))-AIMAG(c(i)))
            delmax = MAX(delmax,delx)
          END DO
          !
          IF( r<=delmax ) THEN
            Nerr = Nerr + 1
            WRITE (Lun,99002) list(kprog), kcase, delmax
            99002 FORMAT ('   PROBLEM WITH C',A,', CASE ',I1,'.  MAX ABS ERROR OF',&
              E11.4/)
          END IF
          !              CHECK CONTROL FOR SINGULAR MATRIX FOR CASE 2
          !
        ELSEIF( ind/=-4 ) THEN
          Nerr = Nerr + 1
          WRITE (Lun,99003) list(kprog), kcase, ind
          99003 FORMAT ('   PROBLEM WITH C',A,', CASE ',I1,'.  IND = ',I2,&
            ' INSTEAD OF -4'/)
        END IF
      END DO
    END DO
    !
    !     SUMMARY PRINT
    !
    IF( Nerr/=0 ) WRITE (Lun,99004) Nerr
    99004 FORMAT (/' **** CQCK DETECTED A TOTAL OF ',I2,' PROBLEMS. ****'/)
    IF( Kprint>=2 .AND. Nerr==0 ) WRITE (Lun,99005)
    99005 FORMAT ('     CQCK DETECTED NO PROBLEMS.'/)
    RETURN
  END SUBROUTINE CQCK
END MODULE TEST22_MOD
!** TEST22
PROGRAM TEST22
  USE TEST22_MOD, ONLY : CQCK, DQCK, SQCK
  USE ISO_FORTRAN_ENV, ONLY : INPUT_UNIT, OUTPUT_UNIT
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !> Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  D2
  !***
  ! **Type:**      ALL (TEST22-A)
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
  !        SNBFS    SNBIR    SPOFS    SPOIR
  !        DNBFS             DPOFS
  !        CNBFS    CNBIR    CPOFS    CPOIR
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  CQCK, DQCK, I1MACH, SQCK, XERMAX, XSETF

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  INTEGER :: kprint, lin, lun, nerr, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST22
  lun = OUTPUT_UNIT
  lin = INPUT_UNIT
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  !
  !     Test LINPACK routines
  !
  CALL SQCK(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL DQCK(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CQCK(lun,kprint,nerr)
  nfail = nfail + nerr
  !
  !     Write PASS or FAIL message
  !
  IF( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST22 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST22 *************')
  END IF
  STOP
END PROGRAM TEST22
