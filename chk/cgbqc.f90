!*==CGBQC.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CGBQC
SUBROUTINE CGBQC(Lun,Kprint,Nerr)
  IMPLICIT NONE
  !*--CGBQC5
  !*** Start of declarations inserted by SPAG
  INTEGER Kprint , Lun
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  CGBQC
  !***PURPOSE  Quick check for CGBFA, CGBCO, CGBSL and CGBDI.
  !***LIBRARY   SLATEC
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Voorhees, E. A., (LANL)
  !***DESCRIPTION
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
  !***ROUTINES CALLED  CGBCO, CGBDI, CGBFA, CGBSL
  !***REVISION HISTORY  (YYMMDD)
  !   801015  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901010  Restructured using IF-THEN-ELSE-ENDIF, moved an ARITHMETIC
  !           STATEMENT FUNCTION ahead of the FIRST EXECUTABLE STATEMENT
  !           record and cleaned up FORMATs.  (RWC)
  !***END PROLOGUE  CGBQC
  COMPLEX abd(6,4) , at(7,4) , b(4) , bt(4) , c(4) , det(2) , dc(2) , z(4) ,&
    xa , xb
  REAL r , rcond , rcnd , DELX
  CHARACTER kfail*39 , kprog*19
  INTEGER lda , n , ipvt(4) , info , i , j , indx , Nerr
  INTEGER ml , mu
  DATA abd/(0.E0,0.E0) , (0.E0,0.E0) , (0.E0,0.E0) , (0.E0,0.E0) ,&
    (2.E0,0.E0) , (0.E0,1.E0) , (0.E0,0.E0) , (0.E0,0.E0) , (0.E0,0.E0) ,&
    (0.E0,-1.E0) , (2.E0,0.E0) , (0.E0,0.E0) , (0.E0,0.E0) , (0.E0,0.E0)&
    , (0.E0,0.E0) , (0.E0,0.E0) , (3.E0,0.E0) , (0.E0,1.E0) , (0.E0,0.E0)&
    , (0.E0,0.E0) , (0.E0,0.E0) , (0.E0,-1.E0) , (4.E0,0.E0) ,&
    (0.E0,0.E0)/
  DATA b/(3.E0,2.E0) , (-1.E0,3.E0) , (0.E0,-4.E0) , (5.E0,0.E0)/
  DATA c/(1.E0,1.E0) , (0.E0,1.E0) , (0.E0,-1.E0) , (1.E0,0.E0)/
  DATA dc/(3.3E0,0.E0) , (1.0E0,0.E0)/
  DATA kprog/'GBFA GBCO GBSL GBDI'/
  DATA kfail/'INFO RCOND SOLUTION DETERMINANT INVERSE'/
  DATA rcnd/.24099E0/
  !***FIRST EXECUTABLE STATEMENT  CGBQC
  lda = 7
  n = 4
  ml = 1
  mu = 3
  Nerr = 0
  !
  !     FORM AT FOR CGBFA AND BT FOR CGBSL, TEST CGBFA
  !
  DO j = 1 , n
    bt(j) = b(j)
    DO i = 1 , 6
      at(i,j) = abd(i,j)
    ENDDO
  ENDDO
  !
  CALL CGBFA(at,lda,n,ml,mu,ipvt,info)
  IF ( info/=0 ) THEN
    WRITE (Lun,99002) kprog(1:4) , kfail(1:4)
    Nerr = Nerr + 1
  ENDIF
  !
  !     TEST CGBSL FOR JOB=0
  !
  CALL CGBSL(at,lda,n,ml,mu,ipvt,bt,0)
  indx = 0
  DO i = 1 , n
    IF ( DELX(c(i),bt(i))>.0001 ) indx = indx + 1
  ENDDO
  !
  IF ( indx/=0 ) THEN
    WRITE (Lun,99002) kprog(11:14) , kfail(12:19)
    Nerr = Nerr + 1
  ENDIF
  !
  !     FORM AT FOR CGBCO AND BT FOR CGBSL, TEST CGBCO
  !
  DO j = 1 , n
    bt(j) = b(j)
    DO i = 1 , 6
      at(i,j) = abd(i,j)
    ENDDO
  ENDDO
  !
  CALL CGBCO(at,lda,n,ml,mu,ipvt,rcond,z)
  r = ABS(rcnd-rcond)
  IF ( r>=.0001 ) THEN
    WRITE (Lun,99002) kprog(6:9) , kfail(6:10)
    Nerr = Nerr + 1
  ENDIF
  !
  !     TEST CGBSL FOR JOB NOT EQUAL TO 0
  !
  CALL CGBSL(at,lda,n,ml,mu,ipvt,bt,1)
  indx = 0
  DO i = 1 , n
    IF ( DELX(c(i),bt(i))>.0001 ) indx = indx + 1
  ENDDO
  !
  IF ( indx/=0 ) THEN
    WRITE (Lun,99002) kprog(11:14) , kfail(12:19)
    Nerr = Nerr + 1
  ENDIF
  !
  !     TEST CGBDI
  !
  CALL CGBDI(at,lda,n,ml,mu,ipvt,det)
  indx = 0
  DO i = 1 , 2
    IF ( DELX(dc(i),det(i))>.0001 ) indx = indx + 1
  ENDDO
  !
  IF ( indx/=0 ) THEN
    WRITE (Lun,99002) kprog(16:19) , kfail(21:31)
    Nerr = Nerr + 1
  ENDIF
  !
  IF ( Kprint>=2.OR.Nerr/=0 ) WRITE (Lun,99001) Nerr
  !
  99001 FORMAT (/' * CGBQC - TEST FOR CGBFA, CGBCO, CGBSL AND CGBDI FOUND ',I1,&
    ' ERRORS.'/)
  RETURN
  99002 FORMAT (/' *** C',A,' FAILURE - ERROR IN ',A)
END SUBROUTINE CGBQC
