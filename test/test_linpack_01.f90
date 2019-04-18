MODULE TEST23_MOD
  IMPLICIT NONE

CONTAINS
  !** CCHQC
  SUBROUTINE CCHQC(Lun,Kprint,Nerr)
    !>
    !  Quick check for CCHDC.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Voorhees, E. A., (LANL)
    !***
    ! **Description:**
    !
    !    QUICK CHECK FOR LINPACK SUBROUTINE CCHDC.
    !
    !    THE CHOLESKY FACTORIZATION OF MATRIX  A  IS COMPARED TO
    !    THE STORED PRE-COMPUTED FACTORIZATION OF  A  (ENTERED
    !    WITH A DATA STATEMENT).  FAILURE OF THE TEST OCCURS WHEN
    !    AGREEMENT TO 3 SIGNIFICANT DIGITS IS NOT ACHIEVED AND AN
    !    ERROR MESSAGE IS PRINTED.
    !
    !    THE INTEGER VALUES OF JPVT AND INFO ARE SIMILARLY TESTED.
    !    LACK OF AGREEMENT RESULTS IN AN ERROR MESSAGE.  A SUMMARY
    !    LINE IS ALWAYS PRINTED.
    !
    !    NO INPUT ARGUMENTS ARE REQUIRED.  ON RETURN, NERR (INTEGER
    !    TYPE) CONTAINS THE TOTAL COUNT OF ALL FAILURES DETECTED.
    !
    !***
    ! **Routines called:**  CCHDC

    !* REVISION HISTORY  (YYMMDD)
    !   801027  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Restructured using IF-THEN-ELSE-ENDIF and cleaned up
    !           FORMATs.  (RWC)
    USE slatec, ONLY : CCHDC
    INTEGER Kprint, Lun, Nerr
    COMPLEX work(4), at(5,4)
    INTEGER lda, p, job, info, jpvtt(4), i, j
    INTEGER indx
    REAL delx
    COMPLEX, PARAMETER :: a(4,4) = RESHAPE( [ &
      (2.E0,0.E0), (0.E0,1.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,-1.E0), (2.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (3.E0,0.E0), (0.E0,1.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (0.E0,-1.E0), (4.E0,0.E0) ], [4,4] )
    INTEGER, PARAMETER :: jpvt(4) = [ 0, -1, 1, 0 ]
    COMPLEX, PARAMETER :: af(4,4) = RESHAPE( [ &
      (1.73205E0,0.E0), (0.E0,1.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,-.57735E0), (1.91485E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (1.41421E0,0.E0), (0.E0,1.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (0.E0,-.70711E0), (1.22475E0,0.E0) ], [4,4] )
    INTEGER, PARAMETER :: infoc = 4
    INTEGER, PARAMETER :: jpvtc(4) = [ 3, 4, 1, 2 ]
    CHARACTER(20), PARAMETER :: kfail = 'FACTORING JPVT INFO '
    !* FIRST EXECUTABLE STATEMENT  CCHQC
    job = 1
    lda = 5
    p = 4
    Nerr = 0
    !
    !     FORM AT AND JPVTT.
    !
    DO j = 1, p
      jpvtt(j) = jpvt(j)
      DO i = 1, p
        at(i,j) = a(i,j)
      END DO
    END DO
    !
    !     TEST CCHDC.
    !
    CALL CCHDC(at,lda,p,work,jpvtt,job,info)
    indx = 0
    DO j = 1, p
      DO i = 1, p
        delx = ABS(REAL(at(i,j)-af(i,j))) + ABS(AIMAG(at(i,j)-af(i,j)))
        IF ( delx>.0001 ) indx = indx + 1
      END DO
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kfail(1:9)
      Nerr = Nerr + 1
    END IF
    !
    indx = 0
    DO i = 1, p
      IF ( jpvtt(i)/=jpvtc(i) ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kfail(11:14)
      Nerr = Nerr + 1
    END IF
    !
    IF ( info/=infoc ) THEN
      WRITE (Lun,99002) kfail(16:19)
      Nerr = Nerr + 1
    END IF
    !
    IF ( Kprint>=2.OR.Nerr/=0 ) WRITE (Lun,99001) Nerr
    !
    99001 FORMAT (/' * CCHQC - TEST FOR CCHDC FOUND ',I1,' ERRORS.'/6X,&
      '(NO TEST FOR CCHUD, CCHDD OR CCHEX)'/)
    RETURN
    99002 FORMAT (/' *** CCHDC FAILURE - ERROR IN ',A)
  END SUBROUTINE CCHQC
  !** CGBQC
  SUBROUTINE CGBQC(Lun,Kprint,Nerr)
    !>
    !  Quick check for CGBFA, CGBCO, CGBSL and CGBDI.
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
    USE slatec, ONLY : CGBCO, CGBDI, CGBFA, CGBSL, CABS1
    INTEGER Kprint, Lun
    COMPLEX at(7,4), bt(4), det(2), z(4)
    REAL r, rcond
    INTEGER lda, n, ipvt(4), info, i, j, indx, Nerr
    INTEGER ml, mu
    COMPLEX, PARAMETER :: abd(6,4) = RESHAPE( [ &
      (0.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), (2.E0,0.E0), (0.E0,1.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), (0.E0,-1.E0), (2.E0,0.E0), (0.E0,0.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), (3.E0,0.E0), (0.E0,1.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), (0.E0,-1.E0), (4.E0,0.E0), (0.E0,0.E0) ], &
      [6,4] )
    COMPLEX, PARAMETER :: b(4) = [ (3.E0,2.E0), (-1.E0,3.E0), (0.E0,-4.E0), (5.E0,0.E0) ]
    COMPLEX, PARAMETER :: c(4) = [ (1.E0,1.E0), (0.E0,1.E0), (0.E0,-1.E0), (1.E0,0.E0) ]
    COMPLEX, PARAMETER :: dc(2) = [ (3.3E0,0.E0), (1.0E0,0.E0) ]
    CHARACTER(19), PARAMETER :: kprog = 'GBFA GBCO GBSL GBDI'
    CHARACTER(39), PARAMETER :: kfail = 'INFO RCOND SOLUTION DETERMINANT INVERSE'
    REAL, PARAMETER :: rcnd = .24099E0
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
    IF ( info/=0 ) THEN
      WRITE (Lun,99002) kprog(1:4), kfail(1:4)
      Nerr = Nerr + 1
    END IF
    !
    !     TEST CGBSL FOR JOB=0
    !
    CALL CGBSL(at,lda,n,ml,mu,ipvt,bt,0)
    indx = 0
    DO i = 1, n
      IF ( CABS1(c(i)-bt(i))>.0001 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(11:14), kfail(12:19)
      Nerr = Nerr + 1
    END IF
    !
    !     FORM AT FOR CGBCO AND BT FOR CGBSL, TEST CGBCO
    !
    DO j = 1, n
      bt(j) = b(j)
      DO i = 1, 6
        at(i,j) = abd(i,j)
      END DO
    END DO
    !
    CALL CGBCO(at,lda,n,ml,mu,ipvt,rcond,z)
    r = ABS(rcnd-rcond)
    IF ( r>=.0001 ) THEN
      WRITE (Lun,99002) kprog(6:9), kfail(6:10)
      Nerr = Nerr + 1
    END IF
    !
    !     TEST CGBSL FOR JOB NOT EQUAL TO 0
    !
    CALL CGBSL(at,lda,n,ml,mu,ipvt,bt,1)
    indx = 0
    DO i = 1, n
      IF ( CABS1(c(i)-bt(i))>.0001 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(11:14), kfail(12:19)
      Nerr = Nerr + 1
    END IF
    !
    !     TEST CGBDI
    !
    CALL CGBDI(at,lda,n,ml,mu,ipvt,det)
    indx = 0
    DO i = 1, 2
      IF ( CABS1(dc(i)-det(i))>.0001 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(16:19), kfail(21:31)
      Nerr = Nerr + 1
    END IF
    !
    IF ( Kprint>=2.OR.Nerr/=0 ) WRITE (Lun,99001) Nerr
    !
    99001 FORMAT (/' * CGBQC - TEST FOR CGBFA, CGBCO, CGBSL AND CGBDI FOUND ',I1,&
      ' ERRORS.'/)
    RETURN
    99002 FORMAT (/' *** C',A,' FAILURE - ERROR IN ',A)
  END SUBROUTINE CGBQC
  !** CGECK
  SUBROUTINE CGECK(Lun,Kprint,Nerr)
    !>
    !  Quick check for CGEFA, CGECO, CGESL and CGEDI.
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
    USE slatec, ONLY : CGECO, CGEDI, CGEFA, CGESL, CABS1
    INTEGER Kprint, Lun
    COMPLEX at(5,4), bt(4), det(2), z(4)
    REAL r, rcond
    INTEGER lda, n, ipvt(4), info, i, j, indx, Nerr
    COMPLEX, PARAMETER :: a(4,4) = RESHAPE( [ &
      (2.E0,0.E0), (0.E0,1.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,-1.E0), (2.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (3.E0,0.E0), (0.E0,1.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (0.E0,-1.E0), (4.E0,0.E0) ], [4,4] )
    COMPLEX, PARAMETER :: b(4) = [ (3.E0,2.E0), (-1.E0,3.E0), (0.E0,-4.E0), (5.E0,0.E0) ]
    COMPLEX, PARAMETER :: c(4) = [ (1.E0,1.E0), (0.E0,1.E0), (0.E0,-1.E0), (1.E0,0.E0) ]
    COMPLEX, PARAMETER :: ainv(4,4) = RESHAPE( [ &
      (.66667E0,0.E0), (0.E0,-.33333E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,.33333E0), (.66667E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (.36364E0,0.E0), (0.E0,-.09091E0), &
      (0.E0,0.E0), (0.E0,0.E0), (0.E0,.09091E0), (.27273E0,0.E0) ], [4,4] )
    COMPLEX, PARAMETER :: dc(2) = [ (3.3E0,0.E0), (1.0E0,0.E0) ]
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
    IF ( info/=0 ) THEN
      WRITE (Lun,99002) kprog(1:4), kfail(1:4)
      Nerr = Nerr + 1
    END IF
    !
    !     TEST CGESL FOR JOB=0
    !
    CALL CGESL(at,lda,n,ipvt,bt,0)
    indx = 0
    DO i = 1, n
      IF ( CABS1(c(i)-bt(i))>.0001 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
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
    IF ( r>=.0001 ) THEN
      WRITE (Lun,99002) kprog(6:9), kfail(6:10)
      Nerr = Nerr + 1
    END IF
    !
    !     TEST CGESL FOR JOB NOT EQUAL TO 0
    !
    CALL CGESL(at,lda,n,ipvt,bt,1)
    indx = 0
    DO i = 1, n
      IF ( CABS1(c(i)-bt(i))>.0001 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(11:14), kfail(12:19)
      Nerr = Nerr + 1
    END IF
    !
    !     TEST CGEDI FOR JOB=11
    !
    CALL CGEDI(at,lda,n,ipvt,det,z,11)
    indx = 0
    DO i = 1, 2
      IF ( CABS1(dc(i)-det(i))>.0001 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(16:19), kfail(21:31)
      Nerr = Nerr + 1
    END IF
    !
    indx = 0
    DO i = 1, n
      DO j = 1, n
        IF ( CABS1(ainv(i,j)-at(i,j))>.0001 ) indx = indx + 1
      END DO
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(16:19), kfail(33:39)
      Nerr = Nerr + 1
    END IF
    !
    IF ( Kprint>=2.OR.Nerr/=0 ) WRITE (Lun,99001) Nerr
    !
    99001 FORMAT (/' * CGECK - TEST FOR CGEFA, CGECO, CGESL AND CGEDI FOUND ',I1,&
      ' ERRORS.'/)
    RETURN
    99002 FORMAT (/' *** C',A,' FAILURE - ERROR IN ',A)
  END SUBROUTINE CGECK
  !** CGTQC
  SUBROUTINE CGTQC(Lun,Kprint,Nerr)
    !>
    !  Quick check for CGTSL.
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
    !    OF THE PROPER TYPE FOR THE LINPACK SUBROUTINE BEING TESTED.
    !    THE VALUES OF  A  AND  B  AND THE PRE-COMPUTED VALUES OF  CX
    !    (THE SOLUTION VECTOR) ARE ENTERED WITH DATA STATEMENTS.
    !
    !    THE COMPUTED VALUES OF  X  ARE COMPARED TO THE STORED
    !    PRE-COMPUTED VALUES OF CX.  FAILURE OF THE TEST OCCURS WHEN
    !    AGREEMENT TO 3 SIGNIFICANT DIGITS IS NOT ACHIEVED AND AN
    !    ERROR MESSAGE IS PRINTED.  A SUMMARY LINE IS ALWAYS PRINTED.
    !
    !    NO INPUT ARGUMENTS ARE REQUIRED.
    !    ON RETURN,  NERR  (INTEGER TYPE) CONTAINS THE TOTAL COUNT
    !    OF ALL FAILURES DETECTED BY CGTQC.
    !
    !***
    ! **Routines called:**  CGTSL

    !* REVISION HISTORY  (YYMMDD)
    !   801024  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Restructured using IF-THEN-ELSE-ENDIF, moved an ARITHMETIC
    !           STATEMENT FUNCTION ahead of the FIRST EXECUTABLE STATEMENT
    !           record and cleaned up FORMATs.  (RWC)
    USE slatec, ONLY : CGTSL, CABS1
    INTEGER Kprint, Lun
    COMPLEX ct(4), dt(4), et(4), bt(4)
    INTEGER n, info, i, indx, Nerr
    REAL delx
    COMPLEX, PARAMETER :: c(4) = [ (0.E0,0.E0), (0.E0,1.E0), (0.E0,0.E0), (0.E0,1.E0) ]
    COMPLEX, PARAMETER :: d(4) = [ (2.E0,0.E0), (2.E0,0.E0), (3.E0,0.E0), (4.E0,0.E0) ]
    COMPLEX, PARAMETER :: e(4) = [ (0.E0,-1.E0), (0.E0,0.E0), (0.E0,-1.E0), (0.E0,0.E0) ]
    COMPLEX, PARAMETER :: b(4) = [ (3.E0,2.E0), (-1.E0,3.E0), (0.E0,-4.E0), (5.E0,0.E0) ]
    COMPLEX, PARAMETER :: cx(4) = [ (1.E0,1.E0), (0.E0,1.E0), (0.E0,-1.E0), (1.E0,0.E0) ]
    CHARACTER(13), PARAMETER :: kfail = 'INFO SOLUTION'
    !* FIRST EXECUTABLE STATEMENT  CGTQC
    n = 4
    Nerr = 0
    DO i = 1, n
      ct(i) = c(i)
      dt(i) = d(i)
      et(i) = e(i)
      bt(i) = b(i)
    END DO
    !
    CALL CGTSL(n,ct,dt,et,bt,info)
    IF ( info/=0 ) THEN
      WRITE (Lun,99002) kfail(1:4)
      Nerr = Nerr + 1
    END IF
    !
    indx = 0
    DO i = 1, n
      delx = CABS1(bt(i)-cx(i))
      IF ( delx>.0001 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kfail(6:13)
      Nerr = Nerr + 1
    END IF
    !
    IF ( Kprint>=2.OR.Nerr/=0 ) WRITE (Lun,99001) Nerr
    !
    99001 FORMAT (/' * CGTQC - TEST FOR CGTSL FOUND ',I1,' ERRORS.'/)
    RETURN
    99002 FORMAT (/' *** CGTSL FAILURE - ERROR IN ',A)
  END SUBROUTINE CGTQC
  !** CHIQC
  SUBROUTINE CHIQC(Lun,Kprint,Nerr)
    !>
    !  Quick check for CHIFA, CHICO, CHISL and CHIDI.
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
    !    ALL FAILURES DETECTED BY CHIQC.
    !
    !***
    ! **Routines called:**  CHICO, CHIDI, CHIFA, CHISL

    !* REVISION HISTORY  (YYMMDD)
    !   801022  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Restructured using IF-THEN-ELSE-ENDIF and cleaned up
    !           FORMATs.  (RWC)
    USE slatec, ONLY : CHICO, CHIDI, CHIFA, CHISL, CABS1
    INTEGER Kprint, Lun
    COMPLEX at(5,4), bt(4), z(4)
    REAL r, rcond, det(2)
    INTEGER lda, n, ipvt(4), info, i, j, indx, Nerr, inert(3)
    COMPLEX, PARAMETER :: a(4,4) = RESHAPE( [ &
      (2.E0,0.E0), (0.E0,1.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,-1.E0), (2.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (3.E0,0.E0), (0.E0,1.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (0.E0,-1.E0), (4.E0,0.E0) ], [4,4] )
    COMPLEX, PARAMETER :: b(4) = [ (3.E0,2.E0), (-1.E0,3.E0), (0.E0,-4.E0), (5.E0,0.E0) ]
    COMPLEX, PARAMETER :: c(4) = [ (1.E0,1.E0), (0.E0,1.E0), (0.E0,-1.E0), (1.E0,0.E0) ]
    COMPLEX, PARAMETER :: ainv(4,4) = RESHAPE( [ &
      (.66667E0,0.E0), (0.E0,1.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,.33333E0), (.66667E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (.36364E0,0.E0), (0.E0,1.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (0.E0,.09091E0), (.27273E0,0.E0) ], [4,4] )
    REAL, PARAMETER :: dc(2) = [ 3.3E0, 1.0E0 ]
    CHARACTER(19), PARAMETER :: kprog = 'HIFA HICO HISL HIDI'
    CHARACTER(47), PARAMETER :: kfail = 'INFO RCOND SOLUTION DETERMINANT INVERSE INERTIA'
    REAL, PARAMETER :: rcnd  = .24099E0
    INTEGER, PARAMETER :: irt(3) = [ 4, 0, 0 ]
    !* FIRST EXECUTABLE STATEMENT  CHIQC
    lda = 5
    n = 4
    Nerr = 0
    !
    !     FORM AT FOR CHIFA AND BT FOR CHISL, TEST CHIFA
    !
    DO j = 1, n
      bt(j) = b(j)
      DO i = 1, n
        at(i,j) = a(i,j)
      END DO
    END DO
    !
    CALL CHIFA(at,lda,n,ipvt,info)
    IF ( info/=0 ) THEN
      WRITE (Lun,99002) kprog(1:4), kfail(1:4)
      Nerr = Nerr + 1
    END IF
    !
    !     TEST CHISL
    !
    CALL CHISL(at,lda,n,ipvt,bt)
    indx = 0
    DO i = 1, n
      IF ( CABS1(c(i)-bt(i))>.0001 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(11:14), kfail(12:19)
      Nerr = Nerr + 1
    END IF
    !
    !     FORM AT FOR CHICO, TEST CHICO
    !
    DO j = 1, n
      DO i = 1, n
        at(i,j) = a(i,j)
      END DO
    END DO
    !
    CALL CHICO(at,lda,n,ipvt,rcond,z)
    r = ABS(rcnd-rcond)
    IF ( r>=.0001 ) THEN
      WRITE (Lun,99002) kprog(6:9), kfail(6:10)
      Nerr = Nerr + 1
    END IF
    !
    !     TEST CHIDI FOR JOB=111
    !
    CALL CHIDI(at,lda,n,ipvt,det,inert,z,111)
    indx = 0
    DO i = 1, 2
      IF ( ABS(dc(i)-det(i))>.0001 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(16:19), kfail(21:31)
      Nerr = Nerr + 1
    END IF
    !
    indx = 0
    DO i = 1, n
      DO j = 1, n
        IF ( CABS1(ainv(i,j)-at(i,j))>.0001 ) indx = indx + 1
      END DO
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(16:19), kfail(33:39)
      Nerr = Nerr + 1
    END IF
    !
    indx = 0
    DO i = 1, 3
      IF ( (inert(i)-irt(i))/=0 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(16:19), kfail(41:47)
      Nerr = Nerr + 1
    END IF
    !
    IF ( Kprint>=2.OR.Nerr/=0 ) WRITE (Lun,99001) Nerr
    !
    99001 FORMAT (/' * CHIQC - TEST FOR CHIFA, CHICO, CHISL AND CHIDI FOUND ',I1,&
      ' ERRORS.'/)
    RETURN
    99002 FORMAT (/' *** C',A,' FAILURE - ERROR IN ',A)
  END SUBROUTINE CHIQC
  !** CHPQC
  SUBROUTINE CHPQC(Lun,Kprint,Nerr)
    !>
    !  Quick check for CHPFA, CHPCO, CHPSL and CHPDI.
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
    !    ALL FAILURES DETECTED BY CHPQC.
    !
    !***
    ! **Routines called:**  CHPCO, CHPDI, CHPFA, CHPSL

    !* REVISION HISTORY  (YYMMDD)
    !   801022  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Restructured using IF-THEN-ELSE-ENDIF and cleaned up
    !           FORMATs.  (RWC)
    USE slatec, ONLY : CHPCO, CHPDI, CHPFA, CHPSL, CABS1
    INTEGER Kprint, Lun
    COMPLEX at(10), bt(4), z(4)
    REAL r, rcond, det(2)
    INTEGER n, ipvt(4), info, i, j, indx, Nerr, inert(3)
    COMPLEX, PARAMETER :: ap(10) = [ (2.E0,0.E0), (0.E0,-1.E0), (2.E0,0.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (3.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,-1.E0), (4.E0,0.E0) ]
    COMPLEX, PARAMETER :: b(4) = [ (3.E0,2.E0), (-1.E0,3.E0), (0.E0,-4.E0), (5.E0,0.E0) ]
    COMPLEX, PARAMETER :: c(4) = [ (1.E0,1.E0), (0.E0,1.E0), (0.E0,-1.E0), (1.E0,0.E0) ]
    COMPLEX, PARAMETER :: ainv(10) = [ (.66667E0,0.E0), (0.E0,.33333E0), &
      (.66667E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), (.36364E0,0.E0), (0.E0,0.E0), &
      (0.E0,0.E0), (0.E0,.09091E0), (.27273E0,0.E0) ]
    REAL, PARAMETER :: dc(2) = [ 3.3E0, 1.0E0 ]
    CHARACTER(19), PARAMETER :: kprog = 'HPFA HPCO HPSL HPDI'
    CHARACTER(47), PARAMETER :: kfail = 'INFO RCOND SOLUTION DETERMINANT INVERSE INERTIA'
    REAL, PARAMETER :: rcnd  = .24099E0
    INTEGER, PARAMETER :: irt(3) = [ 4, 0, 0 ]
    !* FIRST EXECUTABLE STATEMENT  CHPQC
    n = 4
    Nerr = 0
    !
    !     FORM AT FOR CHPFA AND BT FOR CHPSL, TEST CHPFA
    !
    DO j = 1, n
      bt(j) = b(j)
    END DO
    !
    DO i = 1, 10
      at(i) = ap(i)
    END DO
    !
    CALL CHPFA(at,n,ipvt,info)
    IF ( info/=0 ) THEN
      WRITE (Lun,99002) kprog(1:4), kfail(1:4)
      Nerr = Nerr + 1
    END IF
    !
    !     TEST CHPSL
    !
    CALL CHPSL(at,n,ipvt,bt)
    indx = 0
    DO i = 1, n
      IF ( CABS1(c(i)-bt(i))>.0001 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(11:14), kfail(12:19)
      Nerr = Nerr + 1
    END IF
    !
    !     FORM AT FOR CHPCO, TEST CHPCO
    !
    DO i = 1, 10
      at(i) = ap(i)
    END DO
    !
    CALL CHPCO(at,n,ipvt,rcond,z)
    r = ABS(rcnd-rcond)
    IF ( r>=.0001 ) THEN
      WRITE (Lun,99002) kprog(6:9), kfail(6:10)
      Nerr = Nerr + 1
    END IF
    !
    !     TEST CHPDI FOR JOB=111
    !
    CALL CHPDI(at,n,ipvt,det,inert,z,111)
    indx = 0
    DO i = 1, 2
      IF ( ABS(dc(i)-det(i))>.0001 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(16:19), kfail(21:31)
      Nerr = Nerr + 1
    END IF
    !
    indx = 0
    DO i = 1, 10
      IF ( CABS1(ainv(i)-at(i))>.0001 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(16:19), kfail(33:39)
      Nerr = Nerr + 1
    END IF
    !
    indx = 0
    DO i = 1, 3
      IF ( (inert(i)-irt(i))/=0 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(16:19), kfail(41:47)
      Nerr = Nerr + 1
    END IF
    !
    IF ( Kprint>=2.OR.Nerr/=0 ) WRITE (Lun,99001) Nerr
    !
    99001 FORMAT (/' * CHPQC - TEST FOR CHPFA, CHPCO, CHPSL AND CHPDI FOUND ',I1,&
      ' ERRORS.'/)
    RETURN
    99002 FORMAT (/' *** C',A,' FAILURE - ERROR IN ',A)
  END SUBROUTINE CHPQC
  !** CPBQC
  SUBROUTINE CPBQC(Lun,Kprint,Nerr)
    !>
    !  Quick check for CPBFA, CPBCO, CPBSL and CPBDI.
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
    !    ALL FAILURES DETECTED BY CPBQC.
    !
    !***
    ! **Routines called:**  CPBCO, CPBDI, CPBFA, CPBSL

    !* REVISION HISTORY  (YYMMDD)
    !   801020  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Restructured using IF-THEN-ELSE-ENDIF and cleaned up
    !           FORMATs.  (RWC)
    USE slatec, ONLY : CPBCO, CPBDI, CPBFA, CPBSL, CABS1
    INTEGER Kprint, Lun
    COMPLEX at(3,4), bt(4), z(4)
    REAL r, rcond, det(2)
    INTEGER lda, n, info, i, j, indx, Nerr, m
    COMPLEX, PARAMETER :: abd(2,4) = RESHAPE( [ (0.E0,0.E0), (2.E0,0.E0), &
      (0.E0,-1.E0), (2.E0,0.E0),   (0.E0,0.E0), (3.E0,0.E0), &
      (0.E0,-1.E0), (4.E0,0.E0) ], [2,4] )
    COMPLEX, PARAMETER :: b(4) = [ (3.E0,2.E0), (-1.E0,3.E0), (0.E0,-4.E0), (5.E0,0.E0) ]
    COMPLEX, PARAMETER :: c(4) = [ (1.E0,1.E0), (0.E0,1.E0), (0.E0,-1.E0), (1.E0,0.E0) ]
    REAL, PARAMETER :: dc(2) = [ 3.3E0, 1.0E0 ]
    CHARACTER(19), PARAMETER :: kprog = 'PBFA PBCO PBSL PBDI'
    CHARACTER(39), PARAMETER :: kfail = 'INFO RCOND SOLUTION DETERMINANT INVERSE'
    REAL, PARAMETER :: rcnd = .24099E0
    !* FIRST EXECUTABLE STATEMENT  CPBQC
    lda = 3
    n = 4
    m = 1
    Nerr = 0
    !
    !     FORM AT FOR CPBFA AND BT FOR CPBSL, TEST CPBFA
    !
    DO j = 1, n
      bt(j) = b(j)
      DO i = 1, 2
        at(i,j) = abd(i,j)
      END DO
    END DO
    !
    CALL CPBFA(at,lda,n,m,info)
    IF ( info/=0 ) THEN
      WRITE (Lun,99002) kprog(1:4), kfail(1:4)
      Nerr = Nerr + 1
    END IF
    !
    !     TEST CPBSL
    !
    CALL CPBSL(at,lda,n,m,bt)
    indx = 0
    DO i = 1, n
      IF ( CABS1(c(i)-bt(i))>.0001 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(11:14), kfail(12:19)
      Nerr = Nerr + 1
    END IF
    !
    !     FORM AT FOR CPBCO, TEST CPBCO
    !
    DO j = 1, n
      DO i = 1, 2
        at(i,j) = abd(i,j)
      END DO
    END DO
    !
    CALL CPBCO(at,lda,n,m,rcond,z,info)
    r = ABS(rcnd-rcond)
    IF ( r>=.0001 ) THEN
      WRITE (Lun,99002) kprog(6:9), kfail(6:10)
      Nerr = Nerr + 1
    END IF
    !
    IF ( info/=0 ) THEN
      WRITE (Lun,99002) kprog(6:9), kfail(1:4)
      Nerr = Nerr + 1
    END IF
    !
    !     TEST CPBDI
    !
    CALL CPBDI(at,lda,n,m,det)
    indx = 0
    DO i = 1, 2
      IF ( ABS(dc(i)-det(i))>.0001 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(16:19), kfail(21:31)
      Nerr = Nerr + 1
    END IF
    !
    IF ( Kprint>=2.OR.Nerr/=0 ) WRITE (Lun,99001) Nerr
    !
    99001 FORMAT (/' * CPBQC - TEST FOR CPBFA, CPBCO, CPBSL AND CPBDI FOUND ',I1,&
      ' ERRORS.'/)
    RETURN
    99002 FORMAT (/' *** C',A,' FAILURE - ERROR IN ',A)
  END SUBROUTINE CPBQC
  !** CPOQC
  SUBROUTINE CPOQC(Lun,Kprint,Nerr)
    !>
    !  Quick check for CPOFA, CPOCO, CPOSL and CPODI.
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
    USE slatec, ONLY : CPOCO, CPODI, CPOFA, CPOSL, CABS1
    INTEGER Kprint, Lun
    COMPLEX at(5,4), bt(4), z(4)
    REAL r, rcond, det(2)
    INTEGER lda, n, info, i, j, indx, Nerr
    COMPLEX, PARAMETER :: a(4,4) = RESHAPE( [ &
      (2.E0,0.E0), (0.E0,1.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,-1.E0), (2.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (3.E0,0.E0), (0.E0,1.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (0.E0,-1.E0), (4.E0,0.E0) ], [4,4] )
    COMPLEX, PARAMETER :: b(4) = [ (3.E0,2.E0), (-1.E0,3.E0), (0.E0,-4.E0), (5.E0,0.E0) ]
    COMPLEX, PARAMETER :: c(4) = [ (1.E0,1.E0), (0.E0,1.E0), (0.E0,-1.E0), (1.E0,0.E0) ]
    COMPLEX, PARAMETER :: ainv(4,4) = RESHAPE( [ &
      (.66667E0,0.E0), (0.E0,1.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,.33333E0), (.66667E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (.36364E0,0.E0), (0.E0,1.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (0.E0,.09091E0), (.27273E0,0.E0) ], [4,4] )
    REAL, PARAMETER :: dc(2) = [ 3.3E0, 1.0E0 ]
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
    IF ( info/=0 ) THEN
      WRITE (Lun,99002) kprog(1:4), kfail(1:4)
      Nerr = Nerr + 1
    END IF
    !
    !     TEST CPOSL
    !
    CALL CPOSL(at,lda,n,bt)
    indx = 0
    DO i = 1, n
      IF ( CABS1(c(i)-bt(i))>.0001 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
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
    IF ( r>=.0001 ) THEN
      WRITE (Lun,99002) kprog(6:9), kfail(6:10)
      Nerr = Nerr + 1
    END IF
    !
    IF ( info/=0 ) THEN
      WRITE (Lun,99002) kprog(6:9), kfail(1:4)
      Nerr = Nerr + 1
    END IF
    !
    !     TEST CPODI FOR JOB=11
    !
    CALL CPODI(at,lda,n,det,11)
    indx = 0
    DO i = 1, 2
      IF ( ABS(dc(i)-det(i))>.0001 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(16:19), kfail(21:31)
      Nerr = Nerr + 1
    END IF
    !
    indx = 0
    DO i = 1, n
      DO j = 1, n
        IF ( CABS1(ainv(i,j)-at(i,j))>.0001 ) indx = indx + 1
      END DO
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(16:19), kfail(33:39)
      Nerr = Nerr + 1
    END IF
    !
    !
    IF ( Kprint>=2.OR.Nerr/=0 ) WRITE (Lun,99001) Nerr
    !
    99001 FORMAT (/' * CPOQC - TEST FOR CPOFA, CPOCO, CPOSL AND CPODI FOUND ',I1,&
      ' ERRORS.'/)
    RETURN
    99002 FORMAT (/' *** C',A,' FAILURE - ERROR IN ',A)
  END SUBROUTINE CPOQC
  !** CPPQC
  SUBROUTINE CPPQC(Lun,Kprint,Nerr)
    !>
    !  Quick check for CPPFA, CPPCO, CPPSL and CPPDI.
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
    !    ALL FAILURES DETECTED BY CPPQC.
    !
    !***
    ! **Routines called:**  CPPCO, CPPDI, CPPFA, CPPSL

    !* REVISION HISTORY  (YYMMDD)
    !   801016  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Restructured using IF-THEN-ELSE-ENDIF and cleaned up
    !           FORMATs.  (RWC)
    USE slatec, ONLY : CPPCO, CPPDI, CPPFA, CPPSL, CABS1
    INTEGER Kprint, Lun
    COMPLEX at(10), bt(4), z(4)
    REAL r, rcond, det(2)
    INTEGER n, info, i, j, indx, Nerr
    COMPLEX, PARAMETER :: ap(10) = [ (2.E0,0.E0), (0.E0,-1.E0), (2.E0,0.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (3.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,-1.E0), (4.E0,0.E0) ]
    COMPLEX, PARAMETER :: b(4) = [ (3.E0,2.E0), (-1.E0,3.E0), (0.E0,-4.E0), (5.E0,0.E0) ]
    COMPLEX, PARAMETER :: c(4) = [ (1.E0,1.E0), (0.E0,1.E0), (0.E0,-1.E0), (1.E0,0.E0) ]
    COMPLEX, PARAMETER :: ainv(10) = [ (.66667E0,0.E0), (0.E0,.33333E0), &
      (.66667E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), (.36364E0,0.E0), (0.E0,0.E0), &
      (0.E0,0.E0), (0.E0,.09091E0), (.27273E0,0.E0) ]
    REAL, PARAMETER :: dc(2) = [ 3.3E0, 1.0E0 ]
    CHARACTER(19), PARAMETER :: kprog = 'PPFA PPCO PPSL PPDI'
    CHARACTER(39), PARAMETER :: kfail = 'INFO RCOND SOLUTION DETERMINANT INVERSE'
    REAL, PARAMETER :: rcnd = .24099E0
    !* FIRST EXECUTABLE STATEMENT  CPPQC
    n = 4
    Nerr = 0
    !
    !     FORM AT FOR CPPFA AND BT FOR CPPSL, TEST CPPFA
    !
    DO j = 1, n
      bt(j) = b(j)
    END DO
    !
    DO i = 1, 10
      at(i) = ap(i)
    END DO
    !
    CALL CPPFA(at,n,info)
    IF ( info/=0 ) THEN
      WRITE (Lun,99002) kprog(1:4), kfail(1:4)
      Nerr = Nerr + 1
    END IF
    !
    !     TEST CPPSL
    !
    CALL CPPSL(at,n,bt)
    indx = 0
    DO i = 1, n
      IF ( CABS1(c(i)-bt(i))>.0001 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(11:14), kfail(12:19)
      Nerr = Nerr + 1
    END IF
    !
    !     FORM AT FOR CPPCO, TEST CPPCO
    !
    DO i = 1, 10
      at(i) = ap(i)
    END DO
    !
    CALL CPPCO(at,n,rcond,z,info)
    r = ABS(rcnd-rcond)
    IF ( r>=.0001 ) THEN
      WRITE (Lun,99002) kprog(6:9), kfail(6:10)
      Nerr = Nerr + 1
    END IF
    !
    IF ( info/=0 ) THEN
      WRITE (Lun,99002) kprog(6:9), kfail(1:4)
      Nerr = Nerr + 1
    END IF
    !
    !     TEST CPPDI FOR JOB=11
    !
    CALL CPPDI(at,n,det,11)
    indx = 0
    DO i = 1, 2
      IF ( ABS(dc(i)-det(i))>.0001 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(16:19), kfail(21:31)
      Nerr = Nerr + 1
    END IF
    !
    indx = 0
    DO i = 1, 10
      IF ( CABS1(ainv(i)-at(i))>.0001 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(16:19), kfail(33:39)
      Nerr = Nerr + 1
    END IF
    !
    IF ( Kprint>=2.OR.Nerr/=0 ) WRITE (Lun,99001) Nerr
    !
    99001 FORMAT (/' * CPPQC - TEST FOR CPPFA, CPPCO, CPPSL AND CPPDI FOUND ',I1,&
      ' ERRORS.'/)
    RETURN
    99002 FORMAT (/' *** C',A,' FAILURE - ERROR IN ',A)
  END SUBROUTINE CPPQC
  !** CPTQC
  SUBROUTINE CPTQC(Lun,Kprint,Nerr)
    !>
    !  Quick check for CPTSL.
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
    !    OF THE PROPER TYPE FOR THE LINPACK SUBROUTINE BEING TESTED.
    !    THE VALUES OF  A  AND  B  AND THE PRE-COMPUTED VALUES OF  CX
    !    (THE SOLUTION VECTOR) ARE ENTERED WITH DATA STATEMENTS.
    !
    !    THE COMPUTED VALUES OF  X  ARE COMPARED TO THE STORED
    !    PRE-COMPUTED VALUES OF CX.  FAILURE OF THE TEST OCCURS WHEN
    !    AGREEMENT TO 3 SIGNIFICANT DIGITS IS NOT ACHIEVED AND AN
    !    ERROR MESSAGE IS PRINTED.  A SUMMARY LINE IS ALWAYS PRINTED.
    !
    !    NO INPUT ARGUMENTS ARE REQUIRED.
    !    ON RETURN,  NERR  (INTEGER TYPE) CONTAINS THE TOTAL COUNT
    !    OF ALL FAILURES DETECTED BY CPTQC.
    !
    !***
    ! **Routines called:**  CPTSL

    !* REVISION HISTORY  (YYMMDD)
    !   801024  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Restructured using IF-THEN-ELSE-ENDIF and cleaned up
    !           FORMATs.  (RWC)
    USE slatec, ONLY : CPTSL, CABS1
    INTEGER Kprint, Lun
    COMPLEX dt(4), et(4), bt(4)
    INTEGER n, i, indx, Nerr
    REAL delx
    COMPLEX, PARAMETER :: d(4) = [ (2.E0,0.E0), (2.E0,0.E0), (3.E0,0.E0), (4.E0,0.E0) ]
    COMPLEX, PARAMETER :: e(4) = [ (0.E0,-1.E0), (0.E0,0.E0), (0.E0,-1.E0), (0.E0,0.E0) ]
    COMPLEX, PARAMETER :: b(4) = [ (3.E0,2.E0), (-1.E0,3.E0), (0.E0,-4.E0), (5.E0,0.E0) ]
    COMPLEX, PARAMETER :: cx(4) = [ (1.E0,1.E0), (0.E0,1.E0), (0.E0,-1.E0), (1.E0,0.E0) ]
    !* FIRST EXECUTABLE STATEMENT  CPTQC
    n = 4
    Nerr = 0
    DO i = 1, n
      dt(i) = d(i)
      et(i) = e(i)
      bt(i) = b(i)
    END DO
    !
    CALL CPTSL(n,dt,et,bt)
    indx = 0
    DO i = 1, n
      delx = CABS1(bt(i)-cx(i))
      IF ( delx>.0001 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99001)
      99001 FORMAT (/' *** CPTSL FAILURE - ERROR IN SOLUTION')
      Nerr = Nerr + 1
    END IF
    !
    IF ( Kprint>=2.OR.Nerr/=0 ) WRITE (Lun,99002) Nerr
    !
    99002 FORMAT (/' * CPTQC - TEST FOR CPTSL FOUND ',I1,' ERRORS.'/)
    RETURN
  END SUBROUTINE CPTQC
  !** CQRQC
  SUBROUTINE CQRQC(Lun,Kprint,Nerr)
    !>
    !  Quick check for CQRDC and CQRSL.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Voorhees, E. A., (LANL)
    !***
    ! **Description:**
    !
    !    THE RETURNED FLOATING POINT VALUES FROM CQRDC AND CQRSL FOR
    !    FACTORED X, QRAUX, QY, QTY, B, RSD, AND XB ARE COMPARED TO
    !    THEIR CORRESPONDING STORED PRE-COMPUTED VALUES (ENTERED
    !    WITH DATA STATEMENTS).  FAILURE OF THE TEST OCCURS WHEN
    !    AGREEMENT TO 3 SIGNIFICANT DIGITS IS NOT ACHIEVED AND AN
    !    ERROR MESSAGE IS THEN PRINTED.
    !
    !    THE RETURNED INTEGER VALUES OF JPVT AND INFO ARE ALSO CHECKED.
    !    LACK OF AGREEMENT RESULTS IN AN ERROR MESSAGE.  A SUMMARY
    !    LINE IS ALWAYS PRINTED.
    !
    !    NO INPUT ARGUMENTS ARE REQUIRED.  ON RETURN, NERR (INTEGER
    !    TYPE) CONTAINS THE TOTAL COUNT OF ALL FAILURES DETECTED.
    !
    !***
    ! **Routines called:**  CQRDC, CQRSL

    !* REVISION HISTORY  (YYMMDD)
    !   801029  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Restructured using IF-THEN-ELSE-ENDIF, moved an ARITHMETIC
    !           STATEMENT FUNCTION ahead of the FIRST EXECUTABLE STATEMENT
    !           record and cleaned up FORMATs.  (RWC)
    USE slatec, ONLY : CQRDC, CQRSL, CABS1
    INTEGER Kprint, Lun
    COMPLEX qraux(4), work(4), qy(4), qty(4), b(4), rsd(4), xb(4), at(5,4)
    INTEGER ldx, n, p, job, k, infO, jpvtt(4), i, j, indx(5), Nerr, l
    COMPLEX, PARAMETER :: a(4,4) = RESHAPE( [ &
      (2.E0,0.E0), (0.E0,1.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,-1.E0), (2.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,0.E0),  (0.E0,0.E0), (3.E0,0.E0), (0.E0,1.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (0.E0,-1.E0), (4.E0,0.E0) ], [4,4] )
    INTEGER, PARAMETER :: jpvt(4) = [ 0, -1, 1, 0 ]
    COMPLEX, PARAMETER :: y(4) = [ (3.E0,2.E0), (-1.E0,3.E0), (0.E0,-4.E0), (5.E0,0.E0) ]
    COMPLEX, PARAMETER :: ac(4,4) = RESHAPE( [ &
      (-3.16228E0,0.E0), (0.E0,0.E0), (.94868E0,0.E0), (0.E0,.31623E0), &
      (0.E0,2.21359E0), (-3.47851E0,0.E0), (0.E0,.31623E0), (.94868E0,0.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (2.23607E0,0.E0), (0.E0,.70711E0), &
      (0.E0,0.E0), (0.E0,0.E0), (0.E0,-1.78885E0), (-1.34164E0,0.E0) ], [4,4] )
    COMPLEX, PARAMETER :: qrauxc(4) = [ (1.E0,0.E0), (1.E0,0.E0), &
      (1.70711E0,0.E0), (0.E0,0.E0) ]
    INTEGER, PARAMETER :: jpvtc(4) = [ 3, 4, 1, 2 ]
    COMPLEX, PARAMETER :: qyc(4) = [ (0.E0,-5.81378E0), (-2.68328E0,0.E0), &
      (-1.89737E0,-1.58114E0), (1.58114E0,-3.79473E0) ]
    COMPLEX, PARAMETER :: qtyc(4) = [ (0.E0,5.37587E0), (-3.47851E0,0.E0), &
      (4.02492E0,2.23607E0), (0.E0,-1.34164E0) ]
    COMPLEX, PARAMETER :: bc(4) = [ (0.E0,-1.E0), (1.E0,0.E0), (1.E0,1.E0), (0.E0,1.E0) ]
    COMPLEX, PARAMETER :: rsdc(4) = [ (0.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0) ]
    COMPLEX, PARAMETER :: xbc(4) = [ (3.E0,2.E0), (-1.E0,3.E0), (0.E0,-4.E0), (5.E0,0.E0) ]
    CHARACTER(9), PARAMETER :: kprog = 'QRDC QRSL'
    CHARACTER(75), PARAMETER :: kfail = 'FACTOR QRAUX  JPVT  QY        QTY       &
      &SOLUTION  RSD        XB        INFO'
    !* FIRST EXECUTABLE STATEMENT  CQRQC
    ldx = 5
    n = 4
    p = 4
    k = 4
    Nerr = 0
    !
    !     FORM AT AND JPVTT
    !
    DO j = 1, n
      jpvtt(j) = jpvt(j)
      DO i = 1, n
        at(i,j) = a(i,j)
      END DO
    END DO
    !
    !     TEST CQRDC (FACTOR, QRAUX, JPVT)
    !
    job = 1
    CALL CQRDC(at,ldx,n,p,qraux,jpvtt,work,job)
    indx(1) = 0
    DO j = 1, n
      DO i = 1, n
        IF ( CABS1(at(i,j)-ac(i,j))>.0001 ) indx(1) = indx(1) + 1
      END DO
    END DO
    !
    IF ( indx(1)/=0 ) THEN
      WRITE (Lun,99002) kprog(1:4), kfail(1:6)
      Nerr = Nerr + 1
    END IF
    !
    DO i = 1, 2
      indx(i) = 0
    END DO
    !
    DO i = 1, n
      IF ( CABS1(qraux(i)-qrauxc(i))>.0001 ) indx(1) = indx(1) + 1
      IF ( jpvtt(i)/=jpvtc(i) ) indx(2) = indx(2) + 1
    END DO
    !
    DO i = 1, 2
      l = 7*i + 1
      IF ( indx(i)/=0 ) THEN
        WRITE (Lun,99002) kprog(1:4), kfail(l:l+4)
        Nerr = Nerr + 1
      END IF
    END DO
    !
    !     TEST CQRSL (QY, QTY, SOLUTION, RSD, XB, INFO)
    !
    job = 11111
    DO i = 1, 5
      indx(i) = 0
    END DO
    !
    CALL CQRSL(at,ldx,n,k,qraux,y,qy,qty,b,rsd,xb,job,info)
    DO i = 1, n
      IF ( CABS1(qy(i)-qyc(i))>.0001 ) indx(1) = indx(1) + 1
      IF ( CABS1(qty(i)-qtyc(i))>.0001 ) indx(2) = indx(2) + 1
      IF ( CABS1(b(i)-bc(i))>.0001 ) indx(3) = indx(3) + 1
      IF ( CABS1(rsd(i)-rsdc(i))>.0001 ) indx(4) = indx(4) + 1
      IF ( CABS1(xb(i)-xbc(i))>.0001 ) indx(5) = indx(5) + 1
    END DO
    !
    DO i = 1, 5
      l = 10*i + 11
      IF ( indx(i)/=0 ) THEN
        WRITE (Lun,99002) kprog(6:9), kfail(l:l+8)
        Nerr = Nerr + 1
      END IF
    END DO
    !
    IF ( info/=0 ) THEN
      WRITE (Lun,99002) kprog(6:9), kfail(71:74)
      Nerr = Nerr + 1
    END IF
    !
    IF ( Kprint>=2.OR.Nerr/=0 ) WRITE (Lun,99001) Nerr
    !
    99001 FORMAT (/' * CQRQC - TEST FOR CQRDC AND CQRSL FOUND ',I1,' ERRORS.'/)
    RETURN
    99002 FORMAT (/' *** C',A,' FAILURE - ERROR IN ',A)
  END SUBROUTINE CQRQC
  !** CSIQC
  SUBROUTINE CSIQC(Lun,Kprint,Nerr)
    !>
    !  Quick check for CSIFA, CSICO, CSISL and CSIDI.
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
    !    ALL FAILURES DETECTED BY CSIQC.
    !
    !***
    ! **Routines called:**  CSICO, CSIDI, CSIFA, CSISL

    !* REVISION HISTORY  (YYMMDD)
    !   801021  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Restructured using IF-THEN-ELSE-ENDIF and cleaned up
    !           FORMATs.  (RWC)
    USE slatec, ONLY : CSICO, CSIDI, CSIFA, CSISL, CABS1
    INTEGER Kprint, Lun
    COMPLEX at(5,4), bt(4), det(2), z(4)
    REAL r, rcond
    INTEGER lda, n, ipvt(4), info, i, j, indx, Nerr
    COMPLEX, PARAMETER :: a(4,4) = RESHAPE( [ &
      (2.E0,0.E0), (0.E0,1.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,-1.E0), (2.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (3.E0,0.E0), (0.E0,1.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (0.E0,-1.E0), (4.E0,0.E0) ], [4,4] )
    COMPLEX, PARAMETER :: b(4) = [ (3.E0,2.E0), (1.E0,1.E0), (0.E0,-4.E0), (3.E0,0.E0) ]
    COMPLEX, PARAMETER :: c(4) = [ (1.E0,1.E0), (0.E0,1.E0), (0.E0,-1.E0), (1.E0,0.E0) ]
    COMPLEX, PARAMETER :: ainv(4,4) = RESHAPE( [ &
      (.40000E0,0.E0), (0.E0,1.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,.20000E0), (.40000E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (.30769E0,0.E0), (0.E0,1.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (0.E0,.07692E0), (.23077E0,0.E0) ], [4,4] )
    COMPLEX, PARAMETER :: dc(2) = [ (6.5E0,0.E0), (1.0E0,0.E0) ]
    CHARACTER(19), PARAMETER :: kprog = 'SIFA SICO SISL SIDI'
    CHARACTER(39), PARAMETER :: kfail = 'INFO RCOND SOLUTION DETERMINANT INVERSE'
    REAL, PARAMETER :: rcnd = .58692E0
    !* FIRST EXECUTABLE STATEMENT  CSIQC
    lda = 5
    n = 4
    Nerr = 0
    !
    !     FORM AT FOR CSIFA AND BT FOR CSISL, TEST CSIFA
    !
    DO j = 1, n
      bt(j) = b(j)
      DO i = 1, n
        at(i,j) = a(i,j)
      END DO
    END DO
    !
    CALL CSIFA(at,lda,n,ipvt,info)
    IF ( info/=0 ) THEN
      WRITE (Lun,99002) kprog(1:4), kfail(1:4)
      Nerr = Nerr + 1
    END IF
    !
    !     TEST CSISL
    !
    CALL CSISL(at,lda,n,ipvt,bt)
    indx = 0
    DO i = 1, n
      IF ( CABS1(c(i)-bt(i))>.0001 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(11:14), kfail(12:19)
      Nerr = Nerr + 1
    END IF
    !
    !     FORM AT FOR CSICO, TEST CSICO
    !
    DO j = 1, n
      DO i = 1, n
        at(i,j) = a(i,j)
      END DO
    END DO
    !
    CALL CSICO(at,lda,n,ipvt,rcond,z)
    r = ABS(rcnd-rcond)
    IF ( r>=.0001 ) THEN
      WRITE (Lun,99002) kprog(6:9), kfail(6:10)
      Nerr = Nerr + 1
    END IF
    !
    !     TEST CSIDI FOR JOB=11
    !
    CALL CSIDI(at,lda,n,ipvt,det,z,11)
    indx = 0
    DO i = 1, 2
      IF ( CABS1(dc(i)-det(i))>.0001 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(16:19), kfail(21:31)
      Nerr = Nerr + 1
    END IF
    !
    indx = 0
    DO i = 1, n
      DO j = 1, n
        IF ( CABS1(ainv(i,j)-at(i,j))>.0001 ) indx = indx + 1
      END DO
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(16:19), kfail(33:39)
      Nerr = Nerr + 1
    END IF
    !
    IF ( Kprint>=2.OR.Nerr/=0 ) WRITE (Lun,99001) Nerr
    !
    99001 FORMAT (/' * CSIQC - TEST FOR CSIFA, CSICO, CSISL AND CSIDI FOUND ',I1,&
      ' ERRORS.'/)
    RETURN
    99002 FORMAT (/' *** C',A,' FAILURE - ERROR IN ',A)
  END SUBROUTINE CSIQC
  !** CSPQC
  SUBROUTINE CSPQC(Lun,Kprint,Nerr)
    !>
    !  Quick check for CSPFA, CSPCO, CSPSL and CSPDI.
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
    !    ALL FAILURES DETECTED BY CSPQC.
    !
    !***
    ! **Routines called:**  CSPCO, CSPDI, CSPFA, CSPSL

    !* REVISION HISTORY  (YYMMDD)
    !   801021  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Restructured using IF-THEN-ELSE-ENDIF and cleaned up
    !           FORMATs.  (RWC)
    USE slatec, ONLY : CSPCO, CSPDI, CSPFA, CSPSL, CABS1
    INTEGER Kprint, Lun
    COMPLEX at(10), bt(4), det(2), z(4)
    REAL r, rcond
    INTEGER n, ipvt(4), info, i, j, indx, Nerr
    COMPLEX, PARAMETER :: ap(10) = [ (2.E0,0.E0), (0.E0,-1.E0), (2.E0,0.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (3.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,-1.E0), (4.E0,0.E0) ]
    COMPLEX, PARAMETER :: b(4) = [ (3.E0,2.E0), (1.E0,1.E0), (0.E0,-4.E0), (3.E0,0.E0) ]
    COMPLEX, PARAMETER :: c(4) = [ (1.E0,1.E0), (0.E0,1.E0), (0.E0,-1.E0), (1.E0,0.E0) ]
    COMPLEX, PARAMETER :: ainv(10) = [ (.4E0,0.E0), (0.E0,.2E0), (.4E0,0.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (.30769E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,.07692E0), (.23077E0,0.E0) ]
    COMPLEX, PARAMETER :: dc(2) = [ (6.5E0,0.E0), (1.0E0,0.E0) ]
    CHARACTER(19), PARAMETER :: kprog = 'SPFA SPCO SPSL SPDI'
    CHARACTER(39), PARAMETER :: kfail = 'INFO RCOND SOLUTION DETERMINANT INVERSE'
    REAL, PARAMETER :: rcnd = .58692E0
    !* FIRST EXECUTABLE STATEMENT  CSPQC
    n = 4
    Nerr = 0
    !
    !     FORM AT FOR CSPFA AND BT FOR CSPSL, TEST CSPFA
    !
    DO j = 1, n
      bt(j) = b(j)
    END DO
    !
    DO i = 1, 10
      at(i) = ap(i)
    END DO
    !
    CALL CSPFA(at,n,ipvt,info)
    IF ( info/=0 ) THEN
      WRITE (Lun,99002) kprog(1:4), kfail(1:4)
      Nerr = Nerr + 1
    END IF
    !
    !     TEST CSPSL
    !
    CALL CSPSL(at,n,ipvt,bt)
    indx = 0
    DO i = 1, n
      IF ( CABS1(c(i)-bt(i))>.0001 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(11:14), kfail(12:19)
      Nerr = Nerr + 1
    END IF
    !
    !     FORM AT FOR CSPCO, TEST CSPCO
    !
    DO i = 1, 10
      at(i) = ap(i)
    END DO
    !
    CALL CSPCO(at,n,ipvt,rcond,z)
    r = ABS(rcnd-rcond)
    IF ( r>=.0001 ) THEN
      WRITE (Lun,99002) kprog(6:9), kfail(6:10)
      Nerr = Nerr + 1
    END IF
    !
    !     TEST CSPDI FOR JOB=11
    !
    CALL CSPDI(at,n,ipvt,det,z,11)
    indx = 0
    DO i = 1, 2
      IF ( CABS1(dc(i)-det(i))>.0001 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(16:19), kfail(21:31)
      Nerr = Nerr + 1
    END IF
    !
    indx = 0
    DO i = 1, 10
      IF ( CABS1(ainv(i)-at(i))>.0001 ) indx = indx + 1
    END DO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(16:19), kfail(33:39)
      Nerr = Nerr + 1
    END IF
    !
    IF ( Kprint>=2.OR.Nerr/=0 ) WRITE (Lun,99001) Nerr
    !
    99001 FORMAT (/' * CSPQC - TEST FOR CSPFA, CSPCO, CSPSL AND CSPDI FOUND ',I1,&
      ' ERRORS.'/)
    RETURN
    99002 FORMAT (/'*** C',A,' FAILURE - ERROR IN ',A)
  END SUBROUTINE CSPQC
  !** CSVQC
  SUBROUTINE CSVQC(Lun,Kprint,Nerr)
    !>
    !  Quick check for CSVDC.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Voorhees, E. A., (LANL)
    !***
    ! **Description:**
    !
    !    THE RETURNED FLOATING POINT VALUES FROM CSVDC FOR
    !    S, E, U, AND  V  ARE COMPARED TO THEIR
    !    CORRESPONDING STORED PRE-COMPUTED VALUES (ENTERED
    !    WITH DATA STATEMENTS).  FAILURE OF THE TEST OCCURS WHEN
    !    AGREEMENT TO 3 SIGNIFICANT DIGITS IS NOT ACHIEVED AND
    !    AN ERROR MESSAGE IS THEN PRINTED.
    !
    !    THE RETURNED INTEGER VALUE OF INFO IS ALSO CHECKED.
    !    LACK OF AGREEMENT RESULTS IN AN ERROR MESSAGE.  A SUMMARY
    !    LINE IS ALWAYS PRINTED.
    !
    !    NO INPUT ARGUMENTS ARE REQUIRED.  ON RETURN, NERR (INTEGER
    !    TYPE) CONTAINS THE TOTAL COUNT OF ALL FAILURES DETECTED.
    !
    !***
    ! **Routines called:**  CSVDC

    !* REVISION HISTORY  (YYMMDD)
    !   801031  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Restructured using IF-THEN-ELSE-ENDIF, moved an ARITHMETIC
    !           STATEMENT FUNCTION ahead of the FIRST EXECUTABLE STATEMENT
    !           record and cleaned up FORMATs.  (RWC)
    USE slatec, ONLY : CSVDC, CABS1
    INTEGER kone, Kprint, Lun, Nerr
    COMPLEX work(4), s(4), e(4), u(4,4), v(4,4), at(5,4)
    INTEGER ldx, n, p, ldu, ldv, job, info
    INTEGER i, j, indx(4)
    COMPLEX, PARAMETER :: a(4,4) = RESHAPE( [ &
      (2.E0,0.E0), (0.E0,1.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,-1.E0), (2.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (3.E0,0.E0), (0.E0,1.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (0.E0,-1.E0), (4.E0,0.E0) ], [4,4] )
    CHARACTER(12), PARAMETER :: kfail = 'S E U V INFO'
    COMPLEX, PARAMETER :: sc(4) = [ (4.61803E0,0.E0), (3.0E0,0.E0), &
      (2.38197E0,0.E0), (1.E0,0.E0) ]
    COMPLEX, PARAMETER :: ec(4) = [ (0.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0) ]
    COMPLEX, PARAMETER :: uvc(4,4) = RESHAPE( [ &
      (0.E0,0.E0), (0.E0,0.E0), (-.52573E0,0.E0), (0.E0,-.85065E0), &
      (.70711E0,0.E0), (0.E0,.70711E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (-.85065E0,0.E0), (0.E0,.52573E0), &
      (-.70711E0,0.E0), (0.E0,.70711E0), (0.E0,0.E0), (0.E0,0.E0) ], [4,4] )
    !* FIRST EXECUTABLE STATEMENT  CSVQC
    n = 4
    p = 4
    ldx = 5
    ldu = 4
    ldv = 4
    Nerr = 0
    job = 11
    !
    !     FORM AT
    !
    DO j = 1, n
      DO i = 1, n
        at(i,j) = a(i,j)
      END DO
    END DO
    !
    !     TEST CSVDC  (S, E, U, V, INFO)
    !
    DO i = 1, 4
      indx(i) = 0
    END DO
    !
    CALL CSVDC(at,ldx,n,p,s,e,u,ldu,v,ldv,work,job,info)
    DO j = 1, n
      IF ( CABS1(s(j)-sc(j))>.0001 ) indx(1) = indx(1) + 1
      IF ( CABS1(e(j)-ec(j))>.0001 ) indx(2) = indx(2) + 1
      DO i = 1, n
        IF ( CABS1(u(i,j)-uvc(i,j))>.0001 ) indx(3) = indx(3) + 1
        IF ( CABS1(v(i,j)-uvc(i,j))>.0001 ) indx(4) = indx(4) + 1
      END DO
    END DO
    !
    DO i = 1, 4
      kone = 2*i - 1
      IF ( indx(i)/=0 ) THEN
        WRITE (Lun,99002) kfail(kone:kone)
        Nerr = Nerr + 1
      END IF
    END DO
    !
    IF ( info/=0 ) THEN
      WRITE (Lun,99002) kfail(9:12)
      Nerr = Nerr + 1
    END IF
    !
    IF ( Kprint>=2.OR.Nerr/=0 ) WRITE (Lun,99001) Nerr
    !
    99001 FORMAT (/' * CSVQC - TEST FOR CSVDC FOUND ',I1,' ERRORS.'/)
    RETURN
    99002 FORMAT (/' *** CSVQC FAILURE - ERROR IN ',A)
  END SUBROUTINE CSVQC
  !** CTRQC
  SUBROUTINE CTRQC(Lun,Kprint,Nerr)
    !>
    !  Quick check for CTRFA, CTRCO, CTRSL and CTRDI.
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
    !    ALL FAILURES DETECTED BY CTRQC.
    !
    !***
    ! **Routines called:**  CTRCO, CTRDI, CTRSL

    !* REVISION HISTORY  (YYMMDD)
    !   801023  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Restructured using IF-THEN-ELSE-ENDIF and cleaned up
    !           FORMATs.  (RWC)
    USE slatec, ONLY : CTRCO, CTRDI, CTRSL, CABS1
    INTEGER Kprint, Lun
    COMPLEX at(5,4), bt(4), det(2), z(4)
    REAL r, rcond
    INTEGER lda, n, info, i, j, indx, Nerr
    INTEGER job, k, kk
    COMPLEX, PARAMETER :: a(4,4) = RESHAPE( [ &
      (2.E0,0.E0), (0.E0,1.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,-1.E0), (2.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (3.E0,0.E0), (0.E0,1.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (0.E0,-1.E0), (4.E0,0.E0) ], [4,4] )
    COMPLEX, PARAMETER :: b(4,2) = RESHAPE( [ &
      (2.E0,2.E0), (-1.E0,3.E0), (0.E0,-3.E0), (5.E0,0.E0), &
      (3.E0,2.E0), (0.E0,2.E0), (0.E0,-4.E0), (4.E0,0.E0) ], [4,2] )
    COMPLEX, PARAMETER :: c(4) = [ (1.E0,1.E0), (0.E0,1.E0), (0.E0,-1.E0), (1.E0,0.E0) ]
    COMPLEX, PARAMETER :: ainv(4,4,2) = RESHAPE( [ &
      (.50000E0,0.E0), (0.E0,-.25000E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,-1.00000E0), (.50000E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (.33333E0,0.E0), (0.E0,-.083333E0), &
      (0.E0,0.E0), (0.E0,0.E0), (0.E0,-1.00000E0), (.25000E0,0.E0), &
      (.50000E0,0.E0), (0.E0,1.00000E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,.25000E0), (.50000E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), &
      (0.E0,0.E0), (0.E0,0.E0), (.33333E0,0.E0), (0.E0,1.00000E0), &
      (0.E0,0.E0), (0.E0,0.E0), (0.E0,.083333E0), (.25000E0,0.E0) ], [4,4,2] )
    COMPLEX, PARAMETER :: dc(2) = [ (4.8E0,0.E0), (1.0E0,0.E0) ]
    CHARACTER(19), PARAMETER :: kprog = 'TRFA TRCO TRSL TRDI'
    CHARACTER(39), PARAMETER :: kfail  = 'INFO RCOND SOLUTION DETERMINANT INVERSE'
    REAL, PARAMETER :: rcnd(2) = [ .45695E0, .37047E0 ]
    !* FIRST EXECUTABLE STATEMENT  CTRQC
    lda = 5
    n = 4
    Nerr = 0
    !
    !     K=1 FOR LOWER, K=2 FOR UPPER
    !
    DO k = 1, 2
      !
      !        FORM AT FOR CTRCO AND BT FOR CTRSL, TEST CTRCO
      !
      DO j = 1, n
        bt(j) = b(j,k)
        DO i = 1, n
          at(i,j) = a(i,j)
        END DO
      END DO
      !
      job = k - 1
      CALL CTRCO(at,lda,n,rcond,z,job)
      r = ABS(rcnd(k)-rcond)
      IF ( r>=.0001 ) THEN
        WRITE (Lun,99002) kprog(6:9), kfail(6:10)
        Nerr = Nerr + 1
      END IF
      !
      !        TEST CTRSL FOR JOB= 0 OR 1
      !
      CALL CTRSL(at,lda,n,bt,job,info)
      IF ( info/=0 ) THEN
        WRITE (Lun,99002) kprog(11:14), kfail(1:4)
        Nerr = Nerr + 1
      END IF
      !
      indx = 0
      DO i = 1, n
        IF ( CABS1(c(i)-bt(i))>.0001 ) indx = indx + 1
      END DO
      !
      IF ( indx/=0 ) THEN
        WRITE (Lun,99002) kprog(11:14), kfail(12:19)
        Nerr = Nerr + 1
      END IF
      !
      !        FORM BT FOR CTRSL
      !
      kk = 3 - k
      DO j = 1, n
        bt(j) = b(j,kk)
      END DO
      !
      !        TEST CTRSL FOR JOB EQUAL TO 10 OR 11
      !
      job = 9 + k
      CALL CTRSL(at,lda,n,bt,job,info)
      IF ( info/=0 ) THEN
        WRITE (Lun,99002) kprog(11:14), kfail(1:4)
        Nerr = Nerr + 1
      END IF
      !
      indx = 0
      DO i = 1, n
        IF ( CABS1(c(i)-bt(i))>.0001 ) indx = indx + 1
      END DO
      !
      IF ( indx/=0 ) THEN
        WRITE (Lun,99002) kprog(11:14), kfail(12:19)
        Nerr = Nerr + 1
      END IF
      !
      !        TEST CTRDI FOR JOB= 110 OR 111
      !
      job = 109 + k
      CALL CTRDI(at,lda,n,det,job,info)
      IF ( info/=0 ) THEN
        WRITE (Lun,99002) kprog(16:19), kfail(1:4)
        Nerr = Nerr + 1
      END IF
      !
      indx = 0
      DO i = 1, 2
        IF ( CABS1(dc(i)-det(i))>.0001 ) indx = indx + 1
      END DO
      !
      IF ( indx/=0 ) THEN
        WRITE (Lun,99002) kprog(16:19), kfail(21:31)
        Nerr = Nerr + 1
      END IF
      !
      indx = 0
      DO i = 1, n
        DO j = 1, n
          IF ( CABS1(ainv(i,j,k)-at(i,j))>.0001 ) indx = indx + 1
        END DO
      END DO
      !
      IF ( indx/=0 ) THEN
        WRITE (Lun,99002) kprog(16:19), kfail(33:39)
        Nerr = Nerr + 1
      END IF
    END DO
    !
    IF ( Kprint>=2.OR.Nerr/=0 ) WRITE (Lun,99001) Nerr
    !
    99001 FORMAT (/' * CTRQC - TEST FOR CTRCO, CTRSL AND CTRDI FOUND ',I2,&
      ' ERRORS.'/)
    RETURN
    99002 FORMAT (/' *** C',A,' FAILURE - ERROR IN ',A)
  END SUBROUTINE CTRQC
END MODULE TEST23_MOD
!** TEST23
PROGRAM TEST23
  USE TEST23_MOD, ONLY : CCHQC, CGBQC, CGECK, CGTQC, CHIQC, CHPQC, CPBQC, CPOQC, &
    CPPQC, CPTQC, CQRQC, CSIQC, CSPQC, CSVQC, CTRQC
  USE slatec, ONLY : I1MACH, XSETF, XSETUN, XERMAX
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !>
  !  Driver for testing SLATEC subprograms
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
  INTEGER kprint, lin, lun, nerr, nfail
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
  IF ( kprint<=1 ) THEN
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
  CALL CPPQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CPBQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CSIQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CSPQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CHIQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CHPQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CTRQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CGTQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CPTQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CCHQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CQRQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CSVQC(lun,kprint,nerr)
  nfail = nfail + nerr
  !
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST23 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST23 *************')
  END IF
  STOP
END PROGRAM TEST23
