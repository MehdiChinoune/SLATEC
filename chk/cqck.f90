!*==CQCK.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CQCK
SUBROUTINE CQCK(Lun,Kprint,Nerr)
  IMPLICIT NONE
  !*--CQCK5
  !*** Start of declarations inserted by SPAG
  INTEGER Kprint , Lun
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  CQCK
  !***PURPOSE  Quick check for CPOFS, CPOIR, CNBFS and CNBIR.
  !***LIBRARY   SLATEC
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Voorhees, E. A., (LANL)
  !***DESCRIPTION
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
  !***ROUTINES CALLED  CNBFS, CNBIR, CPOFS, CPOIR, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   801002  DATE WRITTEN
  !   891009  Removed unreferenced statement labels.  (WRB)
  !   891009  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901009  Restructured using IF-THEN-ELSE-ENDIF, cleaned up FORMATs,
  !           including removing an illegal character from column 1, and
  !           editorial changes.  (RWC)
  !***END PROLOGUE  CQCK
  REAL r , delx , delmax , R1MACH
  COMPLEX a(4,4) , at(5,4) , abe(5,7) , abet(5,7) , b(4) , bt(4) , c(4) ,&
    work(35)
  CHARACTER(4) :: list(4)
  INTEGER lda , n , ml , mu , ind , iwork(4) , Nerr , i , j , j1 , j2 , jd ,&
    mlp , k , kcase , kprog
  DATA a/(2.E0,0.E0) , (0.E0,1.E0) , (0.E0,0.E0) , (0.E0,0.E0) ,&
    (0.E0,-1.E0) , (2.E0,0.E0) , (0.E0,0.E0) , (0.E0,0.E0) , (0.E0,0.E0)&
    , (0.E0,0.E0) , (3.E0,0.E0) , (0.E0,1.E0) , (0.E0,0.E0) , (0.E0,0.E0)&
    , (0.E0,-1.E0) , (4.E0,0.E0)/
  DATA c/(1.E0,1.E0) , (0.E0,1.E0) , (0.E0,-1.E0) , (1.E0,0.E0)/
  DATA b/(3.E0,2.E0) , (-1.E0,3.E0) , (0.E0,-4.E0) , (5.E0,0.E0)/
  DATA list/'POFS' , 'POIR' , 'NBFS' , 'NBIR'/
  !***FIRST EXECUTABLE STATEMENT  CQCK
  IF ( Kprint>=3 ) WRITE (Lun,99001)
  !
  99001 FORMAT (/' *    CQCK - QUICK CHECK FOR CPOFS, CPOIR, CNBFS AND ','CNBIR'/)
  lda = 5
  n = 4
  ml = 2
  mu = 1
  jd = 2*ml + mu + 1
  Nerr = 0
  r = R1MACH(4)**0.8E0
  !
  !     FORM ABE(NB ARRAY) FROM MATRIX A.
  !
  DO j = 1 , jd
    DO i = 1 , n
      abe(i,j) = (0.0E0,0.0E0)
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
    ENDDO
  ENDDO
  !
  !     CASE 1 FOR WELL-CONDITIONED MATRIX, CASE 2 FOR SINGULAR MATRIX
  !
  DO kcase = 1 , 2
    DO kprog = 1 , 4
      !           FORM BT FROM B, AT FROM A, AND ABET FROM ABE.
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
          at(1,j) = (0.0E0,0.0E0)
        ENDDO
        !
        DO j = 1 , jd
          abet(1,j) = (0.0E0,0.0E0)
        ENDDO
      ENDIF
      !
      !           SOLVE FOR X
      !
      IF ( kprog==1 ) CALL CPOFS(at,lda,n,bt,1,ind,work)
      IF ( kprog==2 ) CALL CPOIR(at,lda,n,bt,1,ind,work)
      IF ( kprog==3 ) CALL CNBFS(abet,lda,n,ml,mu,bt,1,ind,work,iwork)
      IF ( kprog==4 ) CALL CNBIR(abet,lda,n,ml,mu,bt,1,ind,work,iwork)
      !
      !           COMPARE EXACT AND COMPUTED SOLUTIONS FOR CASE 1
      !
      IF ( kcase==1 ) THEN
        delmax = 0.0E0
        DO i = 1 , n
          delx = ABS(REAL(bt(i))-REAL(c(i)))
          delmax = MAX(delmax,delx)
          delx = ABS(AIMAG(bt(i))-AIMAG(c(i)))
          delmax = MAX(delmax,delx)
        ENDDO
        !
        IF ( r<=delmax ) THEN
          Nerr = Nerr + 1
          WRITE (Lun,99002) list(kprog) , kcase , delmax
          99002         FORMAT ('   PROBLEM WITH C',A,', CASE ',I1,'.  MAX ABS ERROR OF',&
            E11.4/)
        ENDIF
        !              CHECK CONTROL FOR SINGULAR MATRIX FOR CASE 2
        !
      ELSEIF ( ind/=-4 ) THEN
        Nerr = Nerr + 1
        WRITE (Lun,99003) list(kprog) , kcase , ind
        99003       FORMAT ('   PROBLEM WITH C',A,', CASE ',I1,'.  IND = ',I2,&
          ' INSTEAD OF -4'/)
      ENDIF
    ENDDO
  ENDDO
  !
  !     SUMMARY PRINT
  !
  IF ( Nerr/=0 ) WRITE (Lun,99004) Nerr
  99004 FORMAT (/' **** CQCK DETECTED A TOTAL OF ',I2,' PROBLEMS. ****'/)
  IF ( Kprint>=2.AND.Nerr==0 ) WRITE (Lun,99005)
  99005 FORMAT ('     CQCK DETECTED NO PROBLEMS.'/)
  RETURN
END SUBROUTINE CQCK
