!*==CPOQC.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CPOQC
      SUBROUTINE CPOQC(Lun,Kprint,Nerr)
      IMPLICIT NONE
!*--CPOQC5
!*** Start of declarations inserted by SPAG
      INTEGER Kprint , Lun
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  CPOQC
!***PURPOSE  Quick check for CPOFA, CPOCO, CPOSL and CPODI.
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
!    ALL FAILURES DETECTED BY CPOQC.
!
!***ROUTINES CALLED  CPOCO, CPODI, CPOFA, CPOSL
!***REVISION HISTORY  (YYMMDD)
!   801016  DATE WRITTEN
!   890618  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   901010  Restructured using IF-THEN-ELSE-ENDIF and cleaned up
!           FORMATs.  (RWC)
!***END PROLOGUE  CPOQC
      COMPLEX a(4,4) , at(5,4) , b(4) , bt(4) , c(4) , ainv(4,4) , z(4) , xa ,
     &        xb
      REAL r , rcond , rcnd , DELX , det(2) , dc(2)
      CHARACTER kprog*19 , kfail*39
      INTEGER lda , n , info , i , j , indx , Nerr
      DATA a/(2.E0,0.E0) , (0.E0,1.E0) , (0.E0,0.E0) , (0.E0,0.E0) ,
     &     (0.E0,-1.E0) , (2.E0,0.E0) , (0.E0,0.E0) , (0.E0,0.E0) , (0.E0,0.E0)
     &     , (0.E0,0.E0) , (3.E0,0.E0) , (0.E0,1.E0) , (0.E0,0.E0) , (0.E0,0.E0)
     &     , (0.E0,-1.E0) , (4.E0,0.E0)/
      DATA b/(3.E0,2.E0) , (-1.E0,3.E0) , (0.E0,-4.E0) , (5.E0,0.E0)/
      DATA c/(1.E0,1.E0) , (0.E0,1.E0) , (0.E0,-1.E0) , (1.E0,0.E0)/
      DATA ainv/(.66667E0,0.E0) , (0.E0,1.E0) , (0.E0,0.E0) , (0.E0,0.E0) ,
     &     (0.E0,.33333E0) , (.66667E0,0.E0) , (0.E0,0.E0) , (0.E0,0.E0) ,
     &     (0.E0,0.E0) , (0.E0,0.E0) , (.36364E0,0.E0) , (0.E0,1.E0) ,
     &     (0.E0,0.E0) , (0.E0,0.E0) , (0.E0,.09091E0) , (.27273E0,0.E0)/
      DATA dc/3.3E0 , 1.0E0/
      DATA kprog/'POFA POCO POSL PODI'/
      DATA kfail/'INFO RCOND SOLUTION DETERMINANT INVERSE'/
      DATA rcnd/.24099E0/
!***FIRST EXECUTABLE STATEMENT  CPOQC
      lda = 5
      n = 4
      Nerr = 0
!
!     FORM AT FOR CPOFA AND BT FOR CPOSL, TEST CPOFA
!
      DO j = 1 , n
        bt(j) = b(j)
        DO i = 1 , n
          at(i,j) = a(i,j)
        ENDDO
      ENDDO
!
      CALL CPOFA(at,lda,n,info)
      IF ( info/=0 ) THEN
        WRITE (Lun,99002) kprog(1:4) , kfail(1:4)
        Nerr = Nerr + 1
      ENDIF
!
!     TEST CPOSL
!
      CALL CPOSL(at,lda,n,bt)
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
!     FORM AT FOR CPOCO, TEST CPOCO
!
      DO j = 1 , n
        DO i = 1 , n
          at(i,j) = a(i,j)
        ENDDO
      ENDDO
!
      CALL CPOCO(at,lda,n,rcond,z,info)
      r = ABS(rcnd-rcond)
      IF ( r>=.0001 ) THEN
        WRITE (Lun,99002) kprog(6:9) , kfail(6:10)
        Nerr = Nerr + 1
      ENDIF
!
      IF ( info/=0 ) THEN
        WRITE (Lun,99002) kprog(6:9) , kfail(1:4)
        Nerr = Nerr + 1
      ENDIF
!
!     TEST CPODI FOR JOB=11
!
      CALL CPODI(at,lda,n,det,11)
      indx = 0
      DO i = 1 , 2
        IF ( ABS(dc(i)-det(i))>.0001 ) indx = indx + 1
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
          IF ( DELX(ainv(i,j),at(i,j))>.0001 ) indx = indx + 1
        ENDDO
      ENDDO
!
      IF ( indx/=0 ) THEN
        WRITE (Lun,99002) kprog(16:19) , kfail(33:39)
        Nerr = Nerr + 1
      ENDIF
!
!
      IF ( Kprint>=2.OR.Nerr/=0 ) WRITE (Lun,99001) Nerr
!
99001 FORMAT (/' * CPOQC - TEST FOR CPOFA, CPOCO, CPOSL AND CPODI FOUND ',I1,
     &        ' ERRORS.'/)
      RETURN
99002 FORMAT (/' *** C',A,' FAILURE - ERROR IN ',A)
      END SUBROUTINE CPOQC
