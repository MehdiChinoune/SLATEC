!*==CGTQC.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CGTQC
SUBROUTINE CGTQC(Lun,Kprint,Nerr)
  IMPLICIT NONE
  !*--CGTQC5
  !*** Start of declarations inserted by SPAG
  INTEGER Kprint, Lun
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  CGTQC
  !***PURPOSE  Quick check for CGTSL.
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
  !    OF ALL FAILURES DETECTED BY CGTQC.
  !
  !***ROUTINES CALLED  CGTSL
  !***REVISION HISTORY  (YYMMDD)
  !   801024  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901010  Restructured using IF-THEN-ELSE-ENDIF, moved an ARITHMETIC
  !           STATEMENT FUNCTION ahead of the FIRST EXECUTABLE STATEMENT
  !           record and cleaned up FORMATs.  (RWC)
  !***END PROLOGUE  CGTQC
  COMPLEX c(4), d(4), e(4), b(4), cx(4), ct(4), dt(4), et(4), bt(4)
  CHARACTER kfail*13
  INTEGER n, info, i, indx, Nerr
  REAL delx
  DATA c/(0.E0,0.E0), (0.E0,1.E0), (0.E0,0.E0), (0.E0,1.E0)/
  DATA d/(2.E0,0.E0), (2.E0,0.E0), (3.E0,0.E0), (4.E0,0.E0)/
  DATA e/(0.E0,-1.E0), (0.E0,0.E0), (0.E0,-1.E0), (0.E0,0.E0)/
  DATA b/(3.E0,2.E0), (-1.E0,3.E0), (0.E0,-4.E0), (5.E0,0.E0)/
  DATA cx/(1.E0,1.E0), (0.E0,1.E0), (0.E0,-1.E0), (1.E0,0.E0)/
  DATA kfail/'INFO SOLUTION'/
  !***FIRST EXECUTABLE STATEMENT  CGTQC
  n = 4
  Nerr = 0
  DO i = 1, n
    ct(i) = c(i)
    dt(i) = d(i)
    et(i) = e(i)
    bt(i) = b(i)
  ENDDO
  !
  CALL CGTSL(n,ct,dt,et,bt,info)
  IF ( info/=0 ) THEN
    WRITE (Lun,99002) kfail(1:4)
    Nerr = Nerr + 1
  ENDIF
  !
  indx = 0
  DO i = 1, n
    delx = ABS(REAL(bt(i)-cx(i))) + ABS(AIMAG(bt(i)-cx(i)))
    IF ( delx>.0001 ) indx = indx + 1
  ENDDO
  !
  IF ( indx/=0 ) THEN
    WRITE (Lun,99002) kfail(6:13)
    Nerr = Nerr + 1
  ENDIF
  !
  IF ( Kprint>=2.OR.Nerr/=0 ) WRITE (Lun,99001) Nerr
  !
  99001 FORMAT (/' * CGTQC - TEST FOR CGTSL FOUND ',I1,' ERRORS.'/)
  RETURN
  99002 FORMAT (/' *** CGTSL FAILURE - ERROR IN ',A)
END SUBROUTINE CGTQC
