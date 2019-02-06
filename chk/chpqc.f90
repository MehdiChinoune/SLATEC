!*==CHPQC.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CHPQC
      SUBROUTINE CHPQC(Lun,Kprint,Nerr)
      IMPLICIT NONE
!*--CHPQC5
!*** Start of declarations inserted by SPAG
      INTEGER Kprint , Lun
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  CHPQC
!***PURPOSE  Quick check for CHPFA, CHPCO, CHPSL and CHPDI.
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
!    ALL FAILURES DETECTED BY CHPQC.
!
!***ROUTINES CALLED  CHPCO, CHPDI, CHPFA, CHPSL
!***REVISION HISTORY  (YYMMDD)
!   801022  DATE WRITTEN
!   890618  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   901010  Restructured using IF-THEN-ELSE-ENDIF and cleaned up
!           FORMATs.  (RWC)
!***END PROLOGUE  CHPQC
      COMPLEX ap(10) , at(10) , b(4) , bt(4) , c(4) , ainv(10) , z(4) , xa , xb
      REAL r , rcond , rcnd , DELX , det(2) , dc(2)
      CHARACTER kprog*19 , kfail*47
      INTEGER n , ipvt(4) , info , i , j , indx , Nerr
      INTEGER inert(3) , irt(3)
      DATA ap/(2.E0,0.E0) , (0.E0,-1.E0) , (2.E0,0.E0) , (0.E0,0.E0) ,
     &     (0.E0,0.E0) , (3.E0,0.E0) , (0.E0,0.E0) , (0.E0,0.E0) , (0.E0,-1.E0)
     &     , (4.E0,0.E0)/
      DATA b/(3.E0,2.E0) , (-1.E0,3.E0) , (0.E0,-4.E0) , (5.E0,0.E0)/
      DATA c/(1.E0,1.E0) , (0.E0,1.E0) , (0.E0,-1.E0) , (1.E0,0.E0)/
      DATA ainv/(.66667E0,0.E0) , (0.E0,.33333E0) , (.66667E0,0.E0) ,
     &     (0.E0,0.E0) , (0.E0,0.E0) , (.36364E0,0.E0) , (0.E0,0.E0) ,
     &     (0.E0,0.E0) , (0.E0,.09091E0) , (.27273E0,0.E0)/
      DATA dc/3.3E0 , 1.0E0/
      DATA kprog/'HPFA HPCO HPSL HPDI'/
      DATA kfail/'INFO RCOND SOLUTION DETERMINANT INVERSE INERTIA'/
      DATA rcnd/.24099E0/
      DATA irt/4 , 0 , 0/
!***FIRST EXECUTABLE STATEMENT  CHPQC
      n = 4
      Nerr = 0
!
!     FORM AT FOR CHPFA AND BT FOR CHPSL, TEST CHPFA
!
      DO j = 1 , n
        bt(j) = b(j)
      ENDDO
!
      DO i = 1 , 10
        at(i) = ap(i)
      ENDDO
!
      CALL CHPFA(at,n,ipvt,info)
      IF ( info/=0 ) THEN
        WRITE (Lun,99002) kprog(1:4) , kfail(1:4)
        Nerr = Nerr + 1
      ENDIF
!
!     TEST CHPSL
!
      CALL CHPSL(at,n,ipvt,bt)
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
!     FORM AT FOR CHPCO, TEST CHPCO
!
      DO i = 1 , 10
        at(i) = ap(i)
      ENDDO
!
      CALL CHPCO(at,n,ipvt,rcond,z)
      r = ABS(rcnd-rcond)
      IF ( r>=.0001 ) THEN
        WRITE (Lun,99002) kprog(6:9) , kfail(6:10)
        Nerr = Nerr + 1
      ENDIF
!
!     TEST CHPDI FOR JOB=111
!
      CALL CHPDI(at,n,ipvt,det,inert,z,111)
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
      DO i = 1 , 10
        IF ( DELX(ainv(i),at(i))>.0001 ) indx = indx + 1
      ENDDO
!
      IF ( indx/=0 ) THEN
        WRITE (Lun,99002) kprog(16:19) , kfail(33:39)
        Nerr = Nerr + 1
      ENDIF
!
      indx = 0
      DO i = 1 , 3
        IF ( (inert(i)-irt(i))/=0 ) indx = indx + 1
      ENDDO
!
      IF ( indx/=0 ) THEN
        WRITE (Lun,99002) kprog(16:19) , kfail(41:47)
        Nerr = Nerr + 1
      ENDIF
!
      IF ( Kprint>=2.OR.Nerr/=0 ) WRITE (Lun,99001) Nerr
!
99001 FORMAT (/' * CHPQC - TEST FOR CHPFA, CHPCO, CHPSL AND CHPDI FOUND ',I1,
     &        ' ERRORS.'/)
      RETURN
99002 FORMAT (/' *** C',A,' FAILURE - ERROR IN ',A)
      END SUBROUTINE CHPQC
