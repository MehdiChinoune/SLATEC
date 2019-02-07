!*==CQRQC.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CQRQC
SUBROUTINE CQRQC(Lun,Kprint,Nerr)
  IMPLICIT NONE
  !*--CQRQC5
  !*** Start of declarations inserted by SPAG
  INTEGER Kprint , Lun
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  CQRQC
  !***PURPOSE  Quick check for CQRDC and CQRSL.
  !***LIBRARY   SLATEC
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Voorhees, E. A., (LANL)
  !***DESCRIPTION
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
  !***ROUTINES CALLED  CQRDC, CQRSL
  !***REVISION HISTORY  (YYMMDD)
  !   801029  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901010  Restructured using IF-THEN-ELSE-ENDIF, moved an ARITHMETIC
  !           STATEMENT FUNCTION ahead of the FIRST EXECUTABLE STATEMENT
  !           record and cleaned up FORMATs.  (RWC)
  !***END PROLOGUE  CQRQC
  COMPLEX a(4,4) , qraux(4) , work(4) , y(4) , qy(4) , qty(4) , b(4) ,&
    rsd(4) , xb(4)
  COMPLEX at(5,4) , ac(4,4) , qrauxc(4) , qyc(4) , qtyc(4) , bc(4) , rsdc(4)&
    , xbc(4) , x1 , x2
  CHARACTER kprog*9 , kfail*75
  INTEGER ldx , n , p , jpvt(4) , job , k , info
  INTEGER jpvtt(4) , jpvtc(4) , i , j , indx(5) , Nerr , l
  REAL DELX
  DATA a/(2.E0,0.E0) , (0.E0,1.E0) , (0.E0,0.E0) , (0.E0,0.E0) ,&
    (0.E0,-1.E0) , (2.E0,0.E0) , (0.E0,0.E0) , (0.E0,0.E0) , (0.E0,0.E0)&
    , (0.E0,0.E0) , (3.E0,0.E0) , (0.E0,1.E0) , (0.E0,0.E0) , (0.E0,0.E0)&
    , (0.E0,-1.E0) , (4.E0,0.E0)/
  DATA jpvt/0 , -1 , 1 , 0/
  DATA y/(3.E0,2.E0) , (-1.E0,3.E0) , (0.E0,-4.E0) , (5.E0,0.E0)/
  DATA ac/(-3.16228E0,0.E0) , (0.E0,0.E0) , (.94868E0,0.E0) ,&
    (0.E0,.31623E0) , (0.E0,2.21359E0) , (-3.47851E0,0.E0) ,&
    (0.E0,.31623E0) , (.94868E0,0.E0) , (0.E0,0.E0) , (0.E0,0.E0) ,&
    (2.23607E0,0.E0) , (0.E0,.70711E0) , (0.E0,0.E0) , (0.E0,0.E0) ,&
    (0.E0,-1.78885E0) , (-1.34164E0,0.E0)/
  DATA qrauxc/(1.E0,0.E0) , (1.E0,0.E0) , (1.70711E0,0.E0) , (0.E0,0.E0)/
  DATA jpvtc/3 , 4 , 1 , 2/
  DATA qyc/(0.E0,-5.81378E0) , (-2.68328E0,0.E0) , (-1.89737E0,-1.58114E0) ,&
    (1.58114E0,-3.79473E0)/
  DATA qtyc/(0.E0,5.37587E0) , (-3.47851E0,0.E0) , (4.02492E0,2.23607E0) ,&
    (0.E0,-1.34164E0)/
  DATA bc/(0.E0,-1.E0) , (1.E0,0.E0) , (1.E0,1.E0) , (0.E0,1.E0)/
  DATA rsdc/(0.E0,0.E0) , (0.E0,0.E0) , (0.E0,0.E0) , (0.E0,0.E0)/
  DATA xbc/(3.E0,2.E0) , (-1.E0,3.E0) , (0.E0,-4.E0) , (5.E0,0.E0)/
  DATA kprog/'QRDC QRSL'/
  DATA kfail/&
    'FACTOR QRAUX  JPVT  QY        QTY       SOLUTION  RSD        XB        INFO'/
  !***FIRST EXECUTABLE STATEMENT  CQRQC
  ldx = 5
  n = 4
  p = 4
  k = 4
  Nerr = 0
  !
  !     FORM AT AND JPVTT
  !
  DO j = 1 , n
    jpvtt(j) = jpvt(j)
    DO i = 1 , n
      at(i,j) = a(i,j)
    ENDDO
  ENDDO
  !
  !     TEST CQRDC (FACTOR, QRAUX, JPVT)
  !
  job = 1
  CALL CQRDC(at,ldx,n,p,qraux,jpvtt,work,job)
  indx(1) = 0
  DO j = 1 , n
    DO i = 1 , n
      IF ( DELX(at(i,j),ac(i,j))>.0001 ) indx(1) = indx(1) + 1
    ENDDO
  ENDDO
  !
  IF ( indx(1)/=0 ) THEN
    WRITE (Lun,99002) kprog(1:4) , kfail(1:6)
    Nerr = Nerr + 1
  ENDIF
  !
  DO i = 1 , 2
    indx(i) = 0
  ENDDO
  !
  DO i = 1 , n
    IF ( DELX(qraux(i),qrauxc(i))>.0001 ) indx(1) = indx(1) + 1
    IF ( jpvtt(i)/=jpvtc(i) ) indx(2) = indx(2) + 1
  ENDDO
  !
  DO i = 1 , 2
    l = 7*i + 1
    IF ( indx(i)/=0 ) THEN
      WRITE (Lun,99002) kprog(1:4) , kfail(l:l+4)
      Nerr = Nerr + 1
    ENDIF
  ENDDO
  !
  !     TEST CQRSL (QY, QTY, SOLUTION, RSD, XB, INFO)
  !
  job = 11111
  DO i = 1 , 5
    indx(i) = 0
  ENDDO
  !
  CALL CQRSL(at,ldx,n,k,qraux,y,qy,qty,b,rsd,xb,job,info)
  DO i = 1 , n
    IF ( DELX(qy(i),qyc(i))>.0001 ) indx(1) = indx(1) + 1
    IF ( DELX(qty(i),qtyc(i))>.0001 ) indx(2) = indx(2) + 1
    IF ( DELX(b(i),bc(i))>.0001 ) indx(3) = indx(3) + 1
    IF ( DELX(rsd(i),rsdc(i))>.0001 ) indx(4) = indx(4) + 1
    IF ( DELX(xb(i),xbc(i))>.0001 ) indx(5) = indx(5) + 1
  ENDDO
  !
  DO i = 1 , 5
    l = 10*i + 11
    IF ( indx(i)/=0 ) THEN
      WRITE (Lun,99002) kprog(6:9) , kfail(l:l+8)
      Nerr = Nerr + 1
    ENDIF
  ENDDO
  !
  IF ( info/=0 ) THEN
    WRITE (Lun,99002) kprog(6:9) , kfail(71:74)
    Nerr = Nerr + 1
  ENDIF
  !
  IF ( Kprint>=2.OR.Nerr/=0 ) WRITE (Lun,99001) Nerr
  !
  99001 FORMAT (/' * CQRQC - TEST FOR CQRDC AND CQRSL FOUND ',I1,' ERRORS.'/)
  RETURN
  99002 FORMAT (/' *** C',A,' FAILURE - ERROR IN ',A)
END SUBROUTINE CQRQC
