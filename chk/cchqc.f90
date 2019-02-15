!DECK CCHQC
SUBROUTINE CCHQC(Lun,Kprint,Nerr)
  IMPLICIT NONE
  INTEGER Kprint, Lun, Nerr
  !***BEGIN PROLOGUE  CCHQC
  !***PURPOSE  Quick check for CCHDC.
  !***LIBRARY   SLATEC
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Voorhees, E. A., (LANL)
  !***DESCRIPTION
  !
  !    QUICK CHECK FOR LINPACK SUBROUTINE CCHDC.
  !
  !    THE CHOLESKY FACTORIZATION OF MATRIX  A  IS COMPARED TO
  !    THE STORED PRE-COMPUTED FACTORIZATION OF  A  (ENTERED
  !    WITH A DATA STATEMENT).  FAILURE OF THE TEST OCCURS WHEN
  !    AGREEMENT TO 3 SIGNIFICANT DIGITS IS NOT ACHIEVED AND AN
  !    ERROR MESSAGE IS PRINTED.
  !
  !    THE INTEGER VALUES OF JPVT AND INFO ARE SIMILARLY TESTED.
  !    LACK OF AGREEMENT RESULTS IN AN ERROR MESSAGE.  A SUMMARY
  !    LINE IS ALWAYS PRINTED.
  !
  !    NO INPUT ARGUMENTS ARE REQUIRED.  ON RETURN, NERR (INTEGER
  !    TYPE) CONTAINS THE TOTAL COUNT OF ALL FAILURES DETECTED.
  !
  !***ROUTINES CALLED  CCHDC
  !***REVISION HISTORY  (YYMMDD)
  !   801027  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901010  Restructured using IF-THEN-ELSE-ENDIF and cleaned up
  !           FORMATs.  (RWC)
  !***END PROLOGUE  CCHQC
  COMPLEX a(4,4), work(4), at(5,4), af(4,4)
  INTEGER lda, p, jpvt(4), job, info, jpvtt(4), i, j, infoc ,&
    jpvtc(4)
  CHARACTER(20) :: kfail
  INTEGER indx
  REAL delx
  DATA a/(2.E0,0.E0), (0.E0,1.E0), (0.E0,0.E0), (0.E0,0.E0) ,&
    (0.E0,-1.E0), (2.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0), (0.E0,0.E0)&
    , (0.E0,0.E0), (3.E0,0.E0), (0.E0,1.E0), (0.E0,0.E0), (0.E0,0.E0)&
    , (0.E0,-1.E0), (4.E0,0.E0)/
  DATA jpvt/0, -1, 1, 0/
  DATA af/(1.73205E0,0.E0), (0.E0,1.E0), (0.E0,0.E0), (0.E0,0.E0) ,&
    (0.E0,-.57735E0), (1.91485E0,0.E0), (0.E0,0.E0), (0.E0,0.E0) ,&
    (0.E0,0.E0), (0.E0,0.E0), (1.41421E0,0.E0), (0.E0,1.E0) ,&
    (0.E0,0.E0), (0.E0,0.E0), (0.E0,-.70711E0), (1.22475E0,0.E0)/
  DATA infoc/4/
  DATA jpvtc/3, 4, 1, 2/
  DATA kfail/'FACTORING JPVT INFO '/
  !***FIRST EXECUTABLE STATEMENT  CCHQC
  job = 1
  lda = 5
  p = 4
  Nerr = 0
  !
  !     FORM AT AND JPVTT.
  !
  DO j = 1, p
    jpvtt(j) = jpvt(j)
    DO i = 1, p
      at(i,j) = a(i,j)
    ENDDO
  ENDDO
  !
  !     TEST CCHDC.
  !
  CALL CCHDC(at,lda,p,work,jpvtt,job,info)
  indx = 0
  DO j = 1, p
    DO i = 1, p
      delx = ABS(REAL(at(i,j)-af(i,j))) + ABS(AIMAG(at(i,j)-af(i,j)))
      IF ( delx>.0001 ) indx = indx + 1
    ENDDO
  ENDDO
  !
  IF ( indx/=0 ) THEN
    WRITE (Lun,99002) kfail(1:9)
    Nerr = Nerr + 1
  ENDIF
  !
  indx = 0
  DO i = 1, p
    IF ( jpvtt(i)/=jpvtc(i) ) indx = indx + 1
  ENDDO
  !
  IF ( indx/=0 ) THEN
    WRITE (Lun,99002) kfail(11:14)
    Nerr = Nerr + 1
  ENDIF
  !
  IF ( info/=infoc ) THEN
    WRITE (Lun,99002) kfail(16:19)
    Nerr = Nerr + 1
  ENDIF
  !
  IF ( Kprint>=2.OR.Nerr/=0 ) WRITE (Lun,99001) Nerr
  !
  99001 FORMAT (/' * CCHQC - TEST FOR CCHDC FOUND ',I1,' ERRORS.'/6X,&
    '(NO TEST FOR CCHUD, CCHDD OR CCHEX)'/)
  RETURN
  99002 FORMAT (/' *** CCHDC FAILURE - ERROR IN ',A)
END SUBROUTINE CCHQC
