!*==EISQX1.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK EISQX1
SUBROUTINE EISQX1(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--EISQX15
  !*** Start of declarations inserted by SPAG
  INTEGER info
  REAL R1MACH
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  EISQX1
  !***PURPOSE  Quick check for SGEEV and CGEEV.
  !***LIBRARY   SLATEC
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     THIS QUICK CHECK ROUTINE IS WRITTEN FOR EISPACK DRIVERS
  !     SGEEV AND CGEEV.  THE EIGENVALUES OF INPUT MATRIX A(.,.)
  !     ARE STORED IN EK(.).  RELERR IS THE RELATIVE ACCURACY
  !     REQUIRED FOR THEM TO PASS.
  !
  !***ROUTINES CALLED  CGEEV, R1MACH, SGEEV
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900405  CALL to XERROR replaced by message to LUN.  (WRB)
  !***END PROLOGUE  EISQX1
  INTEGER Kprint, Ipass, Lun
  INTEGER lda, n, ldv, job, i, j, id
  REAL a(3,3), ek(3), w(9)
  REAL err, erri, relerr, recj
  COMPLEX ac(3,3), ec(3), vc(3,3)
  DATA lda, n, ldv/3*3/
  DATA a/1., -2., 6., -1., 0., -3., 2., 5., 6./
  DATA ek/ - 1., 3., 5./
  !***FIRST EXECUTABLE STATEMENT  EISQX1
  Ipass = 1
  relerr = SQRT(R1MACH(4))
  DO j = 1, n
    DO i = 1, n
      ac(i,j) = CMPLX(a(i,j),0.)
    ENDDO
  ENDDO
  job = 1
  CALL CGEEV(ac,lda,n,ec,vc,ldv,w,job,info)
  IF ( info/=0 ) THEN
    IF ( Kprint>=2 ) WRITE (Lun,99003) 'CGEEV', info
    Ipass = 0
  ENDIF
  DO j = 1, n
    err = ABS(AIMAG(ec(j)))
    IF ( err>=relerr ) Ipass = 0
    recj = REAL(ec(j))
    err = ABS(recj-ek(1))
    id = 1
    DO i = 2, n
      erri = ABS(recj-ek(i))
      IF ( erri<err ) id = i
      err = MIN(erri,err)
    ENDDO
    IF ( ABS(recj-ek(id))/ABS(ek(id))>=relerr ) Ipass = 0
  ENDDO
  job = 0
  CALL SGEEV(a,lda,n,ec,vc,ldv,w,job,info)
  IF ( info/=0 ) THEN
    IF ( Kprint>=2 ) WRITE (Lun,99003) 'SGEEV', info
    Ipass = 0
  ENDIF
  DO j = 1, n
    err = ABS(AIMAG(ec(j)))
    IF ( err>=relerr ) Ipass = 0
    recj = REAL(ec(j))
    err = ABS(recj-ek(1))
    id = 1
    DO i = 2, n
      erri = ABS(recj-ek(i))
      IF ( erri<err ) id = i
      err = MIN(erri,err)
    ENDDO
    IF ( ABS(recj-ek(id))/ABS(ek(id))>=relerr ) Ipass = 0
  ENDDO
  IF ( Kprint>=2.AND.Ipass/=0 ) WRITE (Lun,99001)
  99001 FORMAT (' EISQX1 PASSES ALL TESTS.')
  IF ( Kprint>=1.AND.Ipass==0 ) WRITE (Lun,99002)
  99002 FORMAT (' EISQX1 FAILS SOME TESTS.')
  99003 FORMAT (1X,'Eigenvalue iteration failed to converge in ',A5,', INFO = ',&
    I4)
END SUBROUTINE EISQX1
