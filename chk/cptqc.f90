!*==CPTQC.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CPTQC
      SUBROUTINE CPTQC(Lun,Kprint,Nerr)
      IMPLICIT NONE
!*--CPTQC5
!*** Start of declarations inserted by SPAG
      INTEGER Kprint , Lun
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  CPTQC
!***PURPOSE  Quick check for CPTSL.
!***LIBRARY   SLATEC
!***KEYWORDS  QUICK CHECK
!***AUTHOR  Voorhees, E. A., (LANL)
!***DESCRIPTION
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
!***ROUTINES CALLED  CPTSL
!***REVISION HISTORY  (YYMMDD)
!   801024  DATE WRITTEN
!   890618  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   901010  Restructured using IF-THEN-ELSE-ENDIF and cleaned up
!           FORMATs.  (RWC)
!***END PROLOGUE  CPTQC
      COMPLEX d(4) , e(4) , b(4) , cx(4) , dt(4) , et(4) , bt(4)
      INTEGER n , i , indx , Nerr
      REAL delx
      DATA d/(2.E0,0.E0) , (2.E0,0.E0) , (3.E0,0.E0) , (4.E0,0.E0)/
      DATA e/(0.E0,-1.E0) , (0.E0,0.E0) , (0.E0,-1.E0) , (0.E0,0.E0)/
      DATA b/(3.E0,2.E0) , (-1.E0,3.E0) , (0.E0,-4.E0) , (5.E0,0.E0)/
      DATA cx/(1.E0,1.E0) , (0.E0,1.E0) , (0.E0,-1.E0) , (1.E0,0.E0)/
!***FIRST EXECUTABLE STATEMENT  CPTQC
      n = 4
      Nerr = 0
      DO i = 1 , n
        dt(i) = d(i)
        et(i) = e(i)
        bt(i) = b(i)
      ENDDO
!
      CALL CPTSL(n,dt,et,bt)
      indx = 0
      DO i = 1 , n
        delx = ABS(REAL(bt(i)-cx(i))) + ABS(AIMAG(bt(i)-cx(i)))
        IF ( delx>.0001 ) indx = indx + 1
      ENDDO
!
      IF ( indx/=0 ) THEN
        WRITE (Lun,99001)
99001   FORMAT (/' *** CPTSL FAILURE - ERROR IN SOLUTION')
        Nerr = Nerr + 1
      ENDIF
!
      IF ( Kprint>=2.OR.Nerr/=0 ) WRITE (Lun,99002) Nerr
!
99002 FORMAT (/' * CPTQC - TEST FOR CPTSL FOUND ',I1,' ERRORS.'/)
      RETURN
      END SUBROUTINE CPTQC
