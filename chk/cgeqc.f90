!*==CGEQC.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CGEQC
SUBROUTINE CGEQC(Lun,Kprint,Nerr)
  IMPLICIT NONE
  !*--CGEQC5
  !***BEGIN PROLOGUE  CGEQC
  !***PURPOSE  Quick check for CGEFS and CGEIR.
  !***LIBRARY   SLATEC
  !***TYPE      COMPLEX (SGEQC-S, DGEQC-D, CGEQC-C)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Jacobsen, Nancy, (LANL)
  !***DESCRIPTION
  !
  !   Let A*X=B be a COMPLEX linear system where the
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
  !***ROUTINES CALLED  CGEFS, CGEIR
  !***REVISION HISTORY  (YYMMDD)
  !   801029  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920601  Code reworked and TYPE section added.  (RWC, WRB)
  !***END PROLOGUE  CGEQC
  !     .. Scalar Arguments ..
  INTEGER Kprint , Lun , Nerr
  !     .. Local Scalars ..
  COMPLEX xa , xb
  INTEGER i , ind , indx , itask , j , kprog , lda , n
  !     .. Local Arrays ..
  COMPLEX a(3,3) , atemp(5,3) , b(3) , btemp(3) , bxex(3) , work(12)
  INTEGER iwork(3)
  CHARACTER list(2)*4
  !     .. External Subroutines ..
  EXTERNAL CGEFS , CGEIR
  !     .. Intrinsic Functions ..
  INTRINSIC ABS , AIMAG , REAL
  REAL DELX
  !     .. Data statements ..
  DATA a/(2.,3.) , (1.,1.) , (1.,2.) , (2.,0.) , (1.,-1.) , (0.,0.) ,&
    (0.,0.) , (2.,5.) , (3.,2.)/
  DATA b/(-1.,1.) , (-5.,4.) , (-4.,7.)/
  DATA bxex/(.21459E-01,.209012E+01) , (.261373E+01,-.162231E+01) ,&
    (.785407E+00,.109871E+01)/
  DATA list/'GEFS' , 'GEIR'/
  !***FIRST EXECUTABLE STATEMENT  CGEQC
  n = 3
  lda = 5
  Nerr = 0
  IF ( Kprint>=2 ) WRITE (Lun,99001)
  !
  99001 FORMAT (//,2X,'CGEFS and CGEIR Quick Check'/)
  !
  DO kprog = 1 , 2
    !
    !     First test case - normal
    !
    itask = 1
    DO i = 1 , n
      btemp(i) = b(i)
    ENDDO
    DO j = 1 , n
      DO i = 1 , n
        atemp(i,j) = a(i,j)
      ENDDO
    ENDDO
    IF ( kprog==1 ) THEN
      CALL CGEFS(atemp,lda,n,btemp,itask,ind,work,iwork)
    ELSE
      CALL CGEIR(atemp,lda,n,btemp,itask,ind,work,iwork)
    ENDIF
    IF ( ind<0 ) THEN
      IF ( Kprint>=2 ) WRITE (Lun,FMT=99005) list(kprog) , ind
      Nerr = Nerr + 1
    ENDIF
    !
    !       Calculate error for first test
    !
    indx = 0
    DO i = 1 , n
      IF ( DELX(bxex(i),btemp(i))>.0001 ) indx = indx + 1
    ENDDO
    IF ( indx==0 ) THEN
      IF ( Kprint>=3 ) WRITE (Lun,FMT=99004) list(kprog)
    ELSE
      IF ( Kprint>=2 ) WRITE (Lun,FMT=99005) list(kprog)
      Nerr = Nerr + 1
    ENDIF
    !
    !       Second test case - singular matrix
    !
    itask = 1
    DO i = 1 , n
      btemp(i) = b(i)
    ENDDO
    DO j = 1 , n
      DO i = 1 , n
        atemp(i,j) = a(i,j)
      ENDDO
    ENDDO
    DO j = 1 , n
      atemp(1,j) = (0.E0,0.E0)
    ENDDO
    IF ( kprog==1 ) THEN
      CALL CGEFS(atemp,lda,n,btemp,itask,ind,work,iwork)
    ELSE
      CALL CGEIR(atemp,lda,n,btemp,itask,ind,work,iwork)
    ENDIF
    IF ( ind==-4 ) THEN
      IF ( Kprint>=3 ) WRITE (Lun,FMT=99006) list(kprog)
    ELSE
      IF ( Kprint>=2 ) WRITE (Lun,FMT=99007) list(kprog) , ind
      Nerr = Nerr + 1
    ENDIF
  ENDDO
  !
  IF ( Kprint>=3.AND.Nerr==0 ) WRITE (Lun,99002)
  99002 FORMAT (/,2X,'CGEFS and CGEIR Quick Check PASSED'/)
  IF ( Kprint>=2.AND.Nerr/=0 ) WRITE (Lun,99003)
  99003 FORMAT (/,2X,'CGEFS and CGEIR Quick Check FAILED'/)
  RETURN
  99004 FORMAT (/,5X,'C',A,' Normal test PASSED')
  99005 FORMAT (/,5X,'C',A,' Test FAILED')
  99006 FORMAT (/,5X,'C',A,' Singular test PASSED')
  99007 FORMAT (/,5X,'C',A,' Singular test FAILED, IND=',I3)
END SUBROUTINE CGEQC
