!*==SQCK.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK SQCK
SUBROUTINE SQCK(Lun,Kprint,Nerr)
  IMPLICIT NONE
  !*--SQCK5
  !*** Start of declarations inserted by SPAG
  INTEGER Kprint , Lun
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  SQCK
  !***PURPOSE  Quick check for SPOFS, SPOIR, SNBFS and SNBIR.
  !***LIBRARY   SLATEC
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Voorhees, E. A., (LANL)
  !***DESCRIPTION
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
  !***ROUTINES CALLED  R1MACH, SNBFS, SNBIR, SPOFS, SPOIR
  !***REVISION HISTORY  (YYMMDD)
  !   800930  DATE WRITTEN
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891009  Removed unreferenced statement label.  (WRB)
  !   891009  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901009  Routine writes illegal character to column 1, fixed.
  !           Editorial changes made, code fixed to test all four
  !           routines.  (RWC)
  !   901009  Restructured using IF-THEN-ELSE-ENDIF, cleaned up FORMATs,
  !           including removing an illegal character from column 1, and
  !           fixed code to test all four routines.  (RWC)
  !***END PROLOGUE  SQCK
  REAL a(4,4) , at(5,4) , abe(5,7) , abet(5,7) , b(4) , bt(4) , c(4) ,&
    work(35) , r , delx , delmax , sign , R1MACH
  CHARACTER(4) :: list(4)
  INTEGER lda , n , ml , mu , ind , iwork(4) , Nerr , i , j , j1 , j2 , jd ,&
    mlp , k , kcase , kprog
  DATA a/5.0E0 , 4.0E0 , 1.0E0 , 1.0E0 , 4.0E0 , 5.0E0 , 1.0E0 , 1.0E0 ,&
    1.0E0 , 1.0E0 , 4.0E0 , 2.0E0 , 1.0E0 , 1.0E0 , 2.0E0 , 4.0E0/
  DATA list/'POFS' , 'POIR' , 'NBFS' , 'NBIR'/
  !***FIRST EXECUTABLE STATEMENT  SQCK
  IF ( Kprint>=3 ) WRITE (Lun,99001)
  !
  99001 FORMAT (/' *    SQCK - QUICK CHECK FOR SPOFS, SPOIR, SNBFS AND ','SNBIR'/)
  lda = 5
  n = 4
  ml = 2
  mu = 1
  jd = 2*ml + mu + 1
  Nerr = 0
  r = R1MACH(4)**0.8E0
  !
  !     COMPUTE C VECTOR.
  !
  sign = 1.0E0
  DO i = 1 , n
    c(i) = sign/i
    sign = -sign
  ENDDO
  !
  !     CASE 1 FOR WELL-CONDITIONED MATRIX, CASE 2 FOR SINGULAR MATRIX.
  !
  DO kcase = 1 , 2
    DO kprog = 1 , 4
      !           SET VECTOR B TO ZERO.
      DO i = 1 , n
        b(i) = 0.0E0
      ENDDO
      !
      !           FORM VECTOR B FOR NON-BANDED.
      !
      IF ( kprog<=2 ) THEN
        DO i = 1 , n
          DO j = 1 , n
            b(i) = b(i) + a(i,j)*c(j)
          ENDDO
        ENDDO
      ELSE
        !
        !              FORM ABE(NB ARRAY) FROM MATRIX A
        !              AND FORM VECTOR B FOR BANDED.
        !
        DO j = 1 , jd
          DO i = 1 , n
            abe(i,j) = 0.0E0
          ENDDO
        ENDDO
        !
        mlp = ml + 1
        DO i = 1 , n
          j1 = MAX(1,i-ml)
          j2 = MIN(n,i+mu)
          DO j = j1 , j2
            k = j - i + mlp
            abe(i,k) = a(i,j)
            b(i) = b(i) + (a(i,j)*c(j))
          ENDDO
        ENDDO
      ENDIF
      !
      !           FORM BT FROM B, AT FROM A, AND ABET FROM ABE.
      !
      DO i = 1 , n
        bt(i) = b(i)
        DO j = 1 , n
          at(i,j) = a(i,j)
        ENDDO
      ENDDO
      !
      DO j = 1 , jd
        DO i = 1 , n
          abet(i,j) = abe(i,j)
        ENDDO
      ENDDO
      !
      !           MAKE AT AND ABET SINGULAR FOR CASE  =  2
      !
      IF ( kcase==2 ) THEN
        DO j = 1 , n
          at(1,j) = 0.0E0
        ENDDO
        !
        DO j = 1 , jd
          abet(1,j) = 0.0E0
        ENDDO
      ENDIF
      !
      !           SOLVE FOR X
      !
      IF ( kprog==1 ) CALL SPOFS(at,lda,n,bt,1,ind,work)
      IF ( kprog==2 ) CALL SPOIR(at,lda,n,bt,1,ind,work)
      IF ( kprog==3 ) CALL SNBFS(abet,lda,n,ml,mu,bt,1,ind,work,iwork)
      IF ( kprog==4 ) CALL SNBIR(abet,lda,n,ml,mu,bt,1,ind,work,iwork)
      !
      !           COMPARE EXACT AND COMPUTED SOLUTIONS FOR CASE 1
      !
      IF ( kcase==1 ) THEN
        delmax = 0.0E0
        DO i = 1 , n
          delx = ABS(bt(i)-c(i))
          delmax = MAX(delmax,delx)
        ENDDO
        !
        IF ( r<=delmax ) THEN
          Nerr = Nerr + 1
          WRITE (Lun,99002) list(kprog) , kcase , delmax
          99002         FORMAT ('   PROBLEM WITH S',A,', CASE ',I1,'.  MAX ABS ERROR OF',&
            E11.4/)
        ENDIF
        !
        !              CHECK CONTROL FOR SINGULAR MATRIX FOR CASE 2
        !
      ELSEIF ( ind/=-4 ) THEN
        Nerr = Nerr + 1
        WRITE (Lun,99003) list(kprog) , kcase , ind
        99003       FORMAT ('   PROBLEM WITH S',A,', CASE ',I1,'.  IND = ',I2,&
          ' INSTEAD OF -4'/)
      ENDIF
    ENDDO
  ENDDO
  !
  !     SUMMARY PRINT
  !
  IF ( Nerr/=0 ) WRITE (Lun,99004) Nerr
  99004 FORMAT (/' **** SQCK DETECTED A TOTAL OF ',I2,' PROBLEMS. ****'/)
  IF ( Kprint>=2.AND.Nerr==0 ) WRITE (Lun,99005)
  99005 FORMAT ('     SQCK DETECTED NO PROBLEMS.'/)
  RETURN
END SUBROUTINE SQCK
