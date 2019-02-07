!*==CTRQC.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CTRQC
SUBROUTINE CTRQC(Lun,Kprint,Nerr)
  IMPLICIT NONE
  !*--CTRQC5
  !*** Start of declarations inserted by SPAG
  INTEGER Kprint , Lun
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  CTRQC
  !***PURPOSE  Quick check for CTRFA, CTRCO, CTRSL and CTRDI.
  !***LIBRARY   SLATEC
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Voorhees, E. A., (LANL)
  !***DESCRIPTION
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
  !***ROUTINES CALLED  CTRCO, CTRDI, CTRSL
  !***REVISION HISTORY  (YYMMDD)
  !   801023  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901010  Restructured using IF-THEN-ELSE-ENDIF and cleaned up
  !           FORMATs.  (RWC)
  !***END PROLOGUE  CTRQC
  COMPLEX a(4,4) , at(5,4) , b(4,2) , bt(4) , c(4) , ainv(4,4,2) , det(2) ,&
    dc(2) , z(4) , xa , xb
  REAL r , rcond , rcnd(2) , DELX
  CHARACTER kprog*19 , kfail*39
  INTEGER lda , n , info , i , j , indx , Nerr
  INTEGER job , k , kk
  DATA a/(2.E0,0.E0) , (0.E0,1.E0) , (0.E0,0.E0) , (0.E0,0.E0) ,&
    (0.E0,-1.E0) , (2.E0,0.E0) , (0.E0,0.E0) , (0.E0,0.E0) , (0.E0,0.E0)&
    , (0.E0,0.E0) , (3.E0,0.E0) , (0.E0,1.E0) , (0.E0,0.E0) , (0.E0,0.E0)&
    , (0.E0,-1.E0) , (4.E0,0.E0)/
  DATA b/(2.E0,2.E0) , (-1.E0,3.E0) , (0.E0,-3.E0) , (5.E0,0.E0) ,&
    (3.E0,2.E0) , (0.E0,2.E0) , (0.E0,-4.E0) , (4.E0,0.E0)/
  DATA c/(1.E0,1.E0) , (0.E0,1.E0) , (0.E0,-1.E0) , (1.E0,0.E0)/
  DATA ainv/(.50000E0,0.E0) , (0.E0,-.25000E0) , (0.E0,0.E0) , (0.E0,0.E0) ,&
    (0.E0,-1.00000E0) , (.50000E0,0.E0) , (0.E0,0.E0) , (0.E0,0.E0) ,&
    (0.E0,0.E0) , (0.E0,0.E0) , (.33333E0,0.E0) , (0.E0,-.083333E0) ,&
    (0.E0,0.E0) , (0.E0,0.E0) , (0.E0,-1.00000E0) , (.25000E0,0.E0) ,&
    (.50000E0,0.E0) , (0.E0,1.00000E0) , (0.E0,0.E0) , (0.E0,0.E0) ,&
    (0.E0,.25000E0) , (.50000E0,0.E0) , (0.E0,0.E0) , (0.E0,0.E0) ,&
    (0.E0,0.E0) , (0.E0,0.E0) , (.33333E0,0.E0) , (0.E0,1.00000E0) ,&
    (0.E0,0.E0) , (0.E0,0.E0) , (0.E0,.083333E0) , (.25000E0,0.E0)/
  DATA dc/(4.8E0,0.E0) , (1.0E0,0.E0)/
  DATA kprog/'TRFA TRCO TRSL TRDI'/
  DATA kfail/'INFO RCOND SOLUTION DETERMINANT INVERSE'/
  DATA rcnd/.45695E0 , .37047E0/
  !***FIRST EXECUTABLE STATEMENT  CTRQC
  lda = 5
  n = 4
  Nerr = 0
  !
  !     K=1 FOR LOWER, K=2 FOR UPPER
  !
  DO k = 1 , 2
    !
    !        FORM AT FOR CTRCO AND BT FOR CTRSL, TEST CTRCO
    !
    DO j = 1 , n
      bt(j) = b(j,k)
      DO i = 1 , n
        at(i,j) = a(i,j)
      ENDDO
    ENDDO
    !
    job = k - 1
    CALL CTRCO(at,lda,n,rcond,z,job)
    r = ABS(rcnd(k)-rcond)
    IF ( r>=.0001 ) THEN
      WRITE (Lun,99002) kprog(6:9) , kfail(6:10)
      Nerr = Nerr + 1
    ENDIF
    !
    !        TEST CTRSL FOR JOB= 0 OR 1
    !
    CALL CTRSL(at,lda,n,bt,job,info)
    IF ( info/=0 ) THEN
      WRITE (Lun,99002) kprog(11:14) , kfail(1:4)
      Nerr = Nerr + 1
    ENDIF
    !
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
    !        FORM BT FOR CTRSL
    !
    kk = 3 - k
    DO j = 1 , n
      bt(j) = b(j,kk)
    ENDDO
    !
    !        TEST CTRSL FOR JOB EQUAL TO 10 OR 11
    !
    job = 9 + k
    CALL CTRSL(at,lda,n,bt,job,info)
    IF ( info/=0 ) THEN
      WRITE (Lun,99002) kprog(11:14) , kfail(1:4)
      Nerr = Nerr + 1
    ENDIF
    !
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
    !        TEST CTRDI FOR JOB= 110 OR 111
    !
    job = 109 + k
    CALL CTRDI(at,lda,n,det,job,info)
    IF ( info/=0 ) THEN
      WRITE (Lun,99002) kprog(16:19) , kfail(1:4)
      Nerr = Nerr + 1
    ENDIF
    !
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
    indx = 0
    DO i = 1 , n
      DO j = 1 , n
        IF ( DELX(ainv(i,j,k),at(i,j))>.0001 ) indx = indx + 1
      ENDDO
    ENDDO
    !
    IF ( indx/=0 ) THEN
      WRITE (Lun,99002) kprog(16:19) , kfail(33:39)
      Nerr = Nerr + 1
    ENDIF
  ENDDO
  !
  IF ( Kprint>=2.OR.Nerr/=0 ) WRITE (Lun,99001) Nerr
  !
  99001 FORMAT (/' * CTRQC - TEST FOR CTRCO, CTRSL AND CTRDI FOUND ',I2,&
    ' ERRORS.'/)
  RETURN
  99002 FORMAT (/' *** C',A,' FAILURE - ERROR IN ',A)
END SUBROUTINE CTRQC
