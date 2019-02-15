!DECK DGEQC
SUBROUTINE DGEQC(Lun,Kprint,Nerr)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DGEQC
  !***PURPOSE  Quick check for DGEFS.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (SGEQC-S, DGEQC-D, CGEQC-C)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Jacobsen, Nancy, (LANL)
  !***DESCRIPTION
  !
  !   Let A*X=B be a DOUBLE PRECISION linear system where the
  !   matrix is of the proper type for the Linpack subroutines
  !   being called.  The values of A and B and the pre-computed
  !   values of BXEX (the solution vector) are given in DATA
  !   statements.  The computed test results for X are compared to
  !   the stored pre-computed values.  Failure of the test occurs
  !   when there is less than 80% agreement between the absolute
  !   values.  There are 2 tests - one for the normal case and one
  !   for the singular case.  A message is printed indicating
  !   whether each subroutine has passed or failed for each case.
  !
  !   On return, NERR (INTEGER type) contains the total count of
  !   all failures detected.
  !
  !***ROUTINES CALLED  D1MACH, DGEFS
  !***REVISION HISTORY  (YYMMDD)
  !   801022  DATE WRITTEN
  !   891009  Removed unreferenced statement label.  (WRB)
  !   891009  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920601  Code reworked and TYPE section added.  (RWC, WRB)
  !***END PROLOGUE  DGEQC
  !     .. Scalar Arguments ..
  INTEGER Kprint, Lun, Nerr
  !     .. Local Scalars ..
  REAL(8) :: errcmp, errmax
  INTEGER i, ind, itask, j, kprog, lda, n
  !     .. Local Arrays ..
  REAL(8) :: a(5,4), atemp(5,4), b(4), btemp(4), bxex(4), &
    work(20)
  INTEGER iwork(4)
  CHARACTER list(2)*4
  !     .. External Functions ..
  REAL(8) :: D1MACH
  EXTERNAL D1MACH
  !     .. External Subroutines ..
  EXTERNAL DGEFS
  !     .. Intrinsic Functions ..
  INTRINSIC ABS, MAX
  !     .. Data statements ..
  DATA a/5.0D0, 1.0D0, 0.3D0, 2.1D0, 0.0D0, -1.0D0, -0.5D0, 1.0D0, &
    1.0D0, 0.0D0, 4.5D0, -1.0D0, -1.7D0, 2.0D0, 0.0D0, 0.5D0, &
    2.0D0, 0.6D0, 1.3D0, 0.0D0/
  DATA b/0.0D0, 3.5D0, 3.6D0, 2.4D0/
  DATA bxex/0.10D+01, 0.10D+01, -0.10D+01, 0.10D+01/
  DATA list/'GEFS', 'GEIR'/
  !***FIRST EXECUTABLE STATEMENT  DGEQC
  n = 4
  lda = 5
  Nerr = 0
  errcmp = D1MACH(4)**0.8D0
  IF ( Kprint>=2 ) WRITE (Lun,99001)
  !
  99001 FORMAT (//,2X,'DGEFS Quick Check'/)
  !
  kprog = 1
  !
  !     First test case - normal
  !
  itask = 1
  DO i = 1, n
    btemp(i) = b(i)
  ENDDO
  DO j = 1, n
    DO i = 1, n
      atemp(i,j) = a(i,j)
    ENDDO
  ENDDO
  CALL DGEFS(atemp,lda,n,btemp,itask,ind,work,iwork)
  IF ( ind<0 ) THEN
    IF ( Kprint>=2 ) WRITE (Lun,FMT=99005) list(kprog), ind
    Nerr = Nerr + 1
  ENDIF
  !
  !     Calculate error for first test
  !
  errmax = 0.0D0
  !
  DO i = 1, n
    errmax = MAX(errmax,ABS(btemp(i)-bxex(i)))
  ENDDO
  IF ( errcmp>errmax ) THEN
    IF ( Kprint>=3 ) WRITE (Lun,FMT=99004) list(kprog)
  ELSE
    IF ( Kprint>=2 ) WRITE (Lun,FMT=99005) list(kprog), errmax
    Nerr = Nerr + 1
  ENDIF
  !
  !     Second test case - singular matrix
  !
  itask = 1
  DO i = 1, n
    btemp(i) = b(i)
  ENDDO
  DO j = 1, n
    DO i = 1, n
      atemp(i,j) = a(i,j)
    ENDDO
  ENDDO
  DO j = 1, n
    atemp(1,j) = 0.0D0
  ENDDO
  CALL DGEFS(atemp,lda,n,btemp,itask,ind,work,iwork)
  IF ( ind==-4 ) THEN
    IF ( Kprint>=3 ) WRITE (Lun,FMT=99006) list(kprog)
  ELSE
    IF ( Kprint>=2 ) WRITE (Lun,FMT=99007) list(kprog), ind
    Nerr = Nerr + 1
  ENDIF
  !
  IF ( Kprint>=3.AND.Nerr==0 ) WRITE (Lun,99002)
  99002 FORMAT (/,2X,'DGEFS Quick Check PASSED'/)
  IF ( Kprint>=2.AND.Nerr/=0 ) WRITE (Lun,99003)
  99003 FORMAT (/,2X,'SGEFS and SGEIR Quick Check FAILED'/)
  RETURN
  99004 FORMAT (/,5X,'D',A,' Normal test PASSED')
  99005 FORMAT (/,5X,'D',A,' Test FAILED, MAX ABS(ERROR) is',E13.5)
  99006 FORMAT (/,5X,'D',A,' Singular test PASSED')
  99007 FORMAT (/,5X,'D',A,' Singular test FAILED, IND=',I3)
END SUBROUTINE DGEQC
