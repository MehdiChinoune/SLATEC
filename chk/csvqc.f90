!DECK CSVQC
SUBROUTINE CSVQC(Lun,Kprint,Nerr)
  IMPLICIT NONE
  INTEGER kone, Kprint, Lun, Nerr
  !***BEGIN PROLOGUE  CSVQC
  !***PURPOSE  Quick check for CSVDC.
  !***LIBRARY   SLATEC
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Voorhees, E. A., (LANL)
  !***DESCRIPTION
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
  !***ROUTINES CALLED  CSVDC
  !***REVISION HISTORY  (YYMMDD)
  !   801031  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901010  Restructured using IF-THEN-ELSE-ENDIF, moved an ARITHMETIC
  !           STATEMENT FUNCTION ahead of the FIRST EXECUTABLE STATEMENT
  !           record and cleaned up FORMATs.  (RWC)
  !***END PROLOGUE  CSVQC
  COMPLEX a(4,4), work(4), s(4), e(4), u(4,4), v(4,4)
  COMPLEX at(5,4), sc(4), ec(4), uvc(4,4), x1, x2
  INTEGER ldx, n, p, ldu, ldv, job, info
  CHARACTER kfail*12
  INTEGER i, j, indx(4)
  REAL DELX
  DATA a/(2.E0,0.E0), (0.E0,1.E0), (0.E0,0.E0), (0.E0,0.E0) ,&
    (0.E0,-1.E0), (2.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0)&
    , (0.E0,0.E0), (3.E0,0.E0), (0.E0,1.E0), (0.E0,0.E0), (0.E0,0.E0)&
    , (0.E0,-1.E0), (4.E0,0.E0)/
  DATA kfail/'S E U V INFO'/
  DATA sc/(4.61803E0,0.E0), (3.0E0,0.E0), (2.38197E0,0.E0), (1.E0,0.E0)/
  DATA ec/(0.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0)/
  DATA uvc/(0.E0,0.E0), (0.E0,0.E0), (-.52573E0,0.E0), (0.E0,-.85065E0) ,&
    (.70711E0,0.E0), (0.E0,.70711E0), (0.E0,0.E0), (0.E0,0.E0) ,&
    (0.E0,0.E0), (0.E0,0.E0), (-.85065E0,0.E0), (0.E0,.52573E0) ,&
    (-.70711E0,0.E0), (0.E0,.70711E0), (0.E0,0.E0), (0.E0,0.E0)/
  !***FIRST EXECUTABLE STATEMENT  CSVQC
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
    ENDDO
  ENDDO
  !
  !     TEST CSVDC  (S, E, U, V, INFO)
  !
  DO i = 1, 4
    indx(i) = 0
  ENDDO
  !
  CALL CSVDC(at,ldx,n,p,s,e,u,ldu,v,ldv,work,job,info)
  DO j = 1, n
    IF ( DELX(s(j),sc(j))>.0001 ) indx(1) = indx(1) + 1
    IF ( DELX(e(j),ec(j))>.0001 ) indx(2) = indx(2) + 1
    DO i = 1, n
      IF ( DELX(u(i,j),uvc(i,j))>.0001 ) indx(3) = indx(3) + 1
      IF ( DELX(v(i,j),uvc(i,j))>.0001 ) indx(4) = indx(4) + 1
    ENDDO
  ENDDO
  !
  DO i = 1, 4
    kone = 2*i - 1
    IF ( indx(i)/=0 ) THEN
      WRITE (Lun,99002) kfail(kone:kone)
      Nerr = Nerr + 1
    ENDIF
  ENDDO
  !
  IF ( info/=0 ) THEN
    WRITE (Lun,99002) kfail(9:12)
    Nerr = Nerr + 1
  ENDIF
  !
  IF ( Kprint>=2.OR.Nerr/=0 ) WRITE (Lun,99001) Nerr
  !
  99001 FORMAT (/' * CSVQC - TEST FOR CSVDC FOUND ',I1,' ERRORS.'/)
  RETURN
  99002 FORMAT (/' *** CSVQC FAILURE - ERROR IN ',A)
END SUBROUTINE CSVQC
