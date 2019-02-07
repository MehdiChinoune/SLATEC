!*==CPRPQX.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CPRPQX
SUBROUTINE CPRPQX(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--CPRPQX5
  !*** Start of declarations inserted by SPAG
  REAL R1MACH
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  CPRPQX
  !***PURPOSE  Quick check for CPZERO and RPZERO.
  !***LIBRARY   SLATEC
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Kahaner, D. K., (NBS)
  !***DESCRIPTION
  !
  !     THIS QUICK CHECK ROUTINE IS WRITTEN FOR CPZERO AND RPZERO.
  !     THE ZEROS OF POLYNOMIAL WITH COEFFICIENTS A(.) ARE STORED
  !     IN ZK(.).  RELERR IS THE RELATIVE ACCURACY REQUIRED FOR
  !     THEM TO PASS.
  !
  !***ROUTINES CALLED  CPZERO, R1MACH, RPZERO
  !***REVISION HISTORY  (YYMMDD)
  !   810223  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  CPRPQX
  INTEGER Kprint, Ipass, Lun
  INTEGER ideg, idegp1, info, i, j, id
  REAL a(6), err, erri, relerr
  COMPLEX ac(6), z(5), zk(5), w(21)
  DATA ideg/5/
  DATA a/1., -3.7, 7.4, -10.8, 10.8, -6.8/
  DATA zk/(1.7,0.), (1.,1.), (1.,-1.), (0.,1.4142135623730950488), &
    (0.,-1.4142135623730950488)/
  !***FIRST EXECUTABLE STATEMENT  CPRPQX
  Ipass = 1
  idegp1 = ideg + 1
  relerr = SQRT(R1MACH(4))
  DO j = 1, idegp1
    ac(j) = CMPLX(a(j),0.)
  ENDDO
  info = 0
  CALL CPZERO(ideg,ac,z,w(4),info,w)
  IF ( info/=0 ) THEN
    Ipass = 0
    IF ( info==1.AND.Kprint>=1 ) WRITE (Lun,99001)
    !
    99001   FORMAT (' CPZERO TEST FAILS: LEADING COEFFICIENT OR DEGREE OF',&
      ' POLYNOMIAL IS ZERO')
    IF ( info==2.AND.Kprint>=1 ) WRITE (Lun,99002)
    99002   FORMAT (' CPZERO TEST FAILS: NON-CONVERGENCE IN 125 ITERATIONS')
  ENDIF
  DO j = 1, ideg
    err = ABS(z(j)-zk(1))
    id = 1
    DO i = 2, ideg
      erri = ABS(z(j)-zk(i))
      IF ( erri<err ) id = i
      err = MIN(erri,err)
    ENDDO
    IF ( ABS(z(j)-zk(id))/ABS(zk(id))>=relerr ) Ipass = 0
  ENDDO
  info = 0
  CALL RPZERO(ideg,a,z,w(4),info,w)
  IF ( info/=0 ) THEN
    Ipass = 0
    IF ( info==1.AND.Kprint>=1 ) WRITE (Lun,99003)
    99003   FORMAT (' RPZERO TEST FAILS: LEADING COEFFICIENT OR DEGREE OF',&
      ' POLYNOMIAL IS ZERO')
    IF ( info==2.AND.Kprint>=1 ) WRITE (Lun,99004)
    99004   FORMAT (' RPZERO TEST FAILS: NON-CONVERGENCE IN 125 ITERATIONS')
  ENDIF
  DO j = 1, ideg
    err = ABS(z(j)-zk(1))
    id = 1
    DO i = 2, ideg
      erri = ABS(z(j)-zk(i))
      IF ( erri<err ) id = i
      err = MIN(erri,err)
    ENDDO
    IF ( ABS(z(j)-zk(id))/ABS(zk(id))>=relerr ) Ipass = 0
  ENDDO
  IF ( Kprint>=2.AND.Ipass/=0 ) WRITE (Lun,99005)
  99005 FORMAT (' CPRPQX PASSES ALL TESTS.')
  IF ( Kprint>=1.AND.Ipass==0 ) WRITE (Lun,99006)
  99006 FORMAT (' CPRPQX FAILS SOME TESTS.')
  RETURN
END SUBROUTINE CPRPQX
