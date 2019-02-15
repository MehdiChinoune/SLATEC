!DECK EISQX2
SUBROUTINE EISQX2(Lun,Kprint,Ipass)
  IMPLICIT NONE
  INTEGER info
  REAL R1MACH
  !***BEGIN PROLOGUE  EISQX2
  !***PURPOSE  Quick check for SSIEV, CHIEV and SSPEV.
  !***LIBRARY   SLATEC
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Kahaner, D. K., (NBS)
  !***DESCRIPTION
  !
  !     THIS QUICK CHECK ROUTINE IS WRITTEN FOR EISPACK DRIVERS
  !     SSIEV, CHIEV AND SSPEV.  THE EIGENVALUES OF INPUT MATRIX
  !     A(.,.) ARE STORED IN EK(.).  RELERR IS THE RELATIVE
  !     ACCURACY REQUIRED FOR THEM TO PASS.
  !
  !***ROUTINES CALLED  CHIEV, R1MACH, SSIEV, SSPEV
  !***REVISION HISTORY  (YYMMDD)
  !   800808  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900405  CALL to XERROR replaced by message to LUN.  (WRB)
  !***END PROLOGUE  EISQX2
  INTEGER Kprint, Ipass, Lun
  INTEGER lda, n, ldv, job, i, j, id
  REAL a1(4,4), a2(10), ap(10), e(4), v(4,4), ek(4), w(16)
  REAL err, erri, relerr
  COMPLEX ac(4,4), vc(4,4)
  EQUIVALENCE (v,vc)
  DATA lda, n, ldv/3*4/
  DATA ap/5., 4., 5., 1., 1., 4., 1., 1., 2., 4./
  DATA ek/1., 2., 5., 10./
  !***FIRST EXECUTABLE STATEMENT  EISQX2
  Ipass = 1
  relerr = SQRT(R1MACH(4))
  id = 0
  DO j = 1, n
    DO i = 1, j
      id = id + 1
      a1(i,j) = ap(id)
      a2(id) = ap(id)
      ac(i,j) = CMPLX(ap(id),0.)
    ENDDO
  ENDDO
  job = 1
  CALL CHIEV(ac,lda,n,e,vc,ldv,w,job,info)
  IF ( info/=0 ) THEN
    IF ( Kprint>=2 ) WRITE (Lun,99003) 'CHIEV', info
    Ipass = 0
  ENDIF
  DO j = 1, n
    err = ABS(e(j)-ek(1))
    id = 1
    DO i = 2, n
      erri = ABS(e(j)-ek(i))
      IF ( erri<err ) id = i
      err = MIN(erri,err)
    ENDDO
    IF ( ABS(e(j)-ek(id))/ABS(ek(id))>=relerr ) Ipass = 0
  ENDDO
  CALL SSIEV(a1,lda,n,e,w,job,info)
  IF ( info/=0 ) THEN
    IF ( Kprint>=2 ) WRITE (Lun,99003) 'SSIEV', info
    Ipass = 0
  ENDIF
  DO j = 1, n
    err = ABS(e(j)-ek(1))
    id = 1
    DO i = 2, n
      erri = ABS(e(j)-ek(i))
      IF ( erri<err ) id = i
      err = MIN(erri,err)
    ENDDO
    IF ( ABS(e(j)-ek(id))/ABS(ek(id))>=relerr ) Ipass = 0
  ENDDO
  job = 0
  CALL SSPEV(a2,n,e,v,ldv,w,job,info)
  IF ( info/=0 ) THEN
    IF ( Kprint>=2 ) WRITE (Lun,99003) 'SSPEV', info
    Ipass = 0
  ENDIF
  DO j = 1, n
    err = ABS(e(j)-ek(1))
    id = 1
    DO i = 2, n
      erri = ABS(e(j)-ek(i))
      IF ( erri<err ) id = i
      err = MIN(erri,err)
    ENDDO
    IF ( ABS(e(j)-ek(id))/ABS(ek(id))>=relerr ) Ipass = 0
  ENDDO
  IF ( Kprint>=2.AND.Ipass/=0 ) WRITE (Lun,99001)
  99001 FORMAT (' EISQX2 PASSES ALL TESTS.')
  IF ( Kprint>=1.AND.Ipass==0 ) WRITE (Lun,99002)
  99002 FORMAT (' EISQX2 FAILS SOME TESTS.')
  99003 FORMAT (1X,'Eigenvalue iteration failed to converge in ',A5,', INFO = ',&
    I4)
END SUBROUTINE EISQX2
