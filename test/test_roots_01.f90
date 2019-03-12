MODULE TEST34_MOD
  IMPLICIT NONE

CONTAINS
  !DECK CPRPQX
  SUBROUTINE CPRPQX(Lun,Kprint,Ipass)
    IMPLICIT NONE
    REAL R1MACH
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
      99001 FORMAT (' CPZERO TEST FAILS: LEADING COEFFICIENT OR DEGREE OF',&
        ' POLYNOMIAL IS ZERO')
      IF ( info==2.AND.Kprint>=1 ) WRITE (Lun,99002)
      99002 FORMAT (' CPZERO TEST FAILS: NON-CONVERGENCE IN 125 ITERATIONS')
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
      99003 FORMAT (' RPZERO TEST FAILS: LEADING COEFFICIENT OR DEGREE OF',&
        ' POLYNOMIAL IS ZERO')
      IF ( info==2.AND.Kprint>=1 ) WRITE (Lun,99004)
      99004 FORMAT (' RPZERO TEST FAILS: NON-CONVERGENCE IN 125 ITERATIONS')
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
  !DECK FZTEST
  SUBROUTINE FZTEST(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  FZTEST
    !***PURPOSE  Quick check for FZERO.
    !***LIBRARY   SLATEC
    !***TYPE      SINGLE PRECISION (FZTEST-S, DFZTST-D)
    !***AUTHOR  (UNKNOWN)
    !***ROUTINES CALLED  FZERO, R1MACH, XERCLR, XGETF, XSETF
    !***REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Changed usage of R1MACH(3) to R1MACH(4).  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    !   920212  Code completely restructured to test IFLAG for all values
    !           of KPRINT.  (WRB)
    !***END PROLOGUE  FZTEST
    !     .. Scalar Arguments ..
    INTEGER Ipass, Kprint, Lun
    !     .. Local Scalars ..
    INTEGER iflag, kontrl
    REAL ae, b, c, pi, r, re, tol
    LOGICAL fatal
    !     .. External Functions ..
    REAL R1MACH
    EXTERNAL R1MACH
    !     .. External Subroutines ..
    EXTERNAL FZERO, XERCLR, XGETF, XSETF
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, ATAN, MAX, SIN, SQRT, TAN
    !***FIRST EXECUTABLE STATEMENT  FZTEST
    IF ( Kprint>=2 ) WRITE (Lun,99001)
    99001 FORMAT ('1'/' FZERO QUICK CHECK')
    Ipass = 1
    pi = 4.0E0*ATAN(1.0E0)
    re = 1.0E-6
    ae = 1.0E-6
    tol = MAX(1.0E-5,SQRT(R1MACH(4)))
    !
    !     Set up and solve example problem
    !
    b = 0.1E0
    c = 4.0E0
    r = c - b
    CALL FZERO(SIN,b,c,r,re,ae,iflag)
    !
    !     See if test was passed.
    !
    IF ( ABS(b-pi)<=tol.AND.ABS(c-pi)<=tol ) THEN
      IF ( Kprint>=3 ) WRITE (Lun,99007) 'PASSED', b, c, iflag
    ELSE
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (Lun,99007) 'FAILED', b, c, iflag
    ENDIF
    !
    !     Trigger 2 error conditions
    !
    CALL XGETF(kontrl)
    IF ( Kprint<=2 ) THEN
      CALL XSETF(0)
    ELSE
      CALL XSETF(1)
    ENDIF
    fatal = .FALSE.
    CALL XERCLR
    !
    IF ( Kprint>=3 ) WRITE (Lun,99002)
    99002 FORMAT (/' IFLAG 3 and 4 tests')
    b = 1.0E0
    !
    !     IFLAG=3 (Singular point)
    !
    c = 2.0E0
    r = 0.5E0*(b+c)
    CALL FZERO(TAN,b,c,b,re,ae,iflag)
    IF ( iflag/=3 ) THEN
      Ipass = 0
      fatal = .TRUE.
      IF ( Kprint>=2 ) WRITE (Lun,99008) iflag, 2
    ENDIF
    !
    !     IFLAG=4 (No sign change)
    !
    b = -3.0E0
    c = -0.1E0
    r = 0.5E0*(b+c)
    CALL FZERO(SIN,b,c,r,re,ae,iflag)
    IF ( iflag/=4 ) THEN
      Ipass = 0
      fatal = .TRUE.
      IF ( Kprint>=2 ) WRITE (Lun,99008) iflag, 4
    ENDIF
    !
    CALL XERCLR
    !
    CALL XSETF(kontrl)
    IF ( fatal ) THEN
      IF ( Kprint>=2 ) THEN
        WRITE (Lun,99003)
        99003 FORMAT (/' At least IFLAG test failed')
      ENDIF
    ELSEIF ( Kprint>=3 ) THEN
      WRITE (Lun,99004)
      99004 FORMAT (/' All IFLAG tests passed')
    ENDIF
    !
    IF ( Ipass==1.AND.Kprint>=2 ) WRITE (Lun,99005)
    99005 FORMAT (/' ***************FZERO PASSED ALL TESTS**************')
    IF ( Ipass==0.AND.Kprint>=1 ) WRITE (Lun,99006)
    99006 FORMAT (/' ***************FZERO FAILED SOME TESTS*************')
    RETURN
    99007 FORMAT (' Accuracy test ',&
      A/' Example problem results:  (answer = PI),  B =',F20.14,' C =',&
      F20.14/' IFLAG =',I2)
    99008 FORMAT (/' IFLAG test FAILED.  IFLAG =',I2,', but should ','have been',I2)
  END SUBROUTINE FZTEST
  !DECK DFZTST
  SUBROUTINE DFZTST(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  DFZTST
    !***PURPOSE  Quick check for DFZERO.
    !***LIBRARY   SLATEC
    !***TYPE      DOUBLE PRECISION (FZTEST-S, DFZTST-D)
    !***AUTHOR  Boland, W. Robert, (LANL)
    !***ROUTINES CALLED  D1MACH, DFZERO, XERCLR, XGETF, XSETF
    !***REVISION HISTORY  (YYMMDD)
    !   920212  DATE WRITTEN
    !***END PROLOGUE  DFZTST
    !     .. Scalar Arguments ..
    INTEGER Ipass, Kprint, Lun
    !     .. Local Scalars ..
    INTEGER iflag, kontrl
    REAL(8) :: ae, b, c, pi, r, re, tol
    LOGICAL fatal
    !     .. External Functions ..
    REAL(8) :: D1MACH
    EXTERNAL D1MACH
    !     .. External Subroutines ..
    EXTERNAL DFZERO, XERCLR, XGETF, XSETF
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, ATAN, DSIN, DTAN, MAX, SQRT
    !***FIRST EXECUTABLE STATEMENT  DFZTST
    IF ( Kprint>=2 ) WRITE (Lun,99001)
    99001 FORMAT ('1'/' DFZERO QUICK CHECK')
    Ipass = 1
    pi = 4.0D0*ATAN(1.0D0)
    re = 1.0D-10
    ae = 1.0D-10
    tol = MAX(1.0D-9,SQRT(D1MACH(4)))
    !
    !     Set up and solve example problem
    !
    b = 0.1D0
    c = 4.0D0
    r = c - b
    CALL DFZERO(DSIN,b,c,r,re,ae,iflag)
    !
    !     See if test was passed.
    !
    IF ( ABS(b-pi)<=tol.AND.ABS(c-pi)<=tol ) THEN
      IF ( Kprint>=3 ) WRITE (Lun,99007) 'PASSED', b, c, iflag
    ELSE
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (Lun,99007) 'FAILED', b, c, iflag
    ENDIF
    !
    !     Trigger 2 error conditions
    !
    CALL XGETF(kontrl)
    IF ( Kprint<=2 ) THEN
      CALL XSETF(0)
    ELSE
      CALL XSETF(1)
    ENDIF
    fatal = .FALSE.
    CALL XERCLR
    !
    IF ( Kprint>=3 ) WRITE (Lun,99002)
    99002 FORMAT (/' IFLAG 3 and 4 tests')
    b = 1.0D0
    !
    !     IFLAG=3 (Singular point)
    !
    c = 2.0D0
    r = 0.5D0*(b+c)
    CALL DFZERO(DTAN,b,c,b,re,ae,iflag)
    IF ( iflag/=3 ) THEN
      Ipass = 0
      fatal = .TRUE.
      IF ( Kprint>=2 ) WRITE (Lun,99008) iflag, 2
    ENDIF
    !
    !     IFLAG=4 (No sign change)
    !
    b = -3.0D0
    c = -0.1D0
    r = 0.5D0*(b+c)
    CALL DFZERO(DSIN,b,c,r,re,ae,iflag)
    IF ( iflag/=4 ) THEN
      Ipass = 0
      fatal = .TRUE.
      IF ( Kprint>=2 ) WRITE (Lun,99008) iflag, 4
    ENDIF
    !
    CALL XERCLR
    !
    CALL XSETF(kontrl)
    IF ( fatal ) THEN
      IF ( Kprint>=2 ) THEN
        WRITE (Lun,99003)
        99003 FORMAT (/' At least IFLAG test failed')
      ENDIF
    ELSEIF ( Kprint>=3 ) THEN
      WRITE (Lun,99004)
      99004 FORMAT (/' All IFLAG tests passed')
    ENDIF
    !
    IF ( Ipass==1.AND.Kprint>=2 ) WRITE (Lun,99005)
    99005 FORMAT (/' ***************DFZERO PASSED ALL TESTS**************')
    IF ( Ipass==0.AND.Kprint>=1 ) WRITE (Lun,99006)
    99006 FORMAT (/' ***************DFZERO FAILED SOME TESTS*************')
    RETURN
    99007 FORMAT (' Accuracy test ',&
      A/' Example problem results:  (answer = PI),  B =',F20.14,' C =',&
      F20.14/' IFLAG =',I2)
    99008 FORMAT (/' IFLAG test FAILED.  IFLAG =',I2,', but should ','have been',I2)
  END SUBROUTINE DFZTST
  !DECK RQRTST
  SUBROUTINE RQRTST(Lun,Kprint,Ipass)
    IMPLICIT NONE
    REAL beta, R1MACH, tol, work
    INTEGER i, ierr, Ipass, j, kontrl, Kprint, Lun, nerr, NUMXER
    !***BEGIN PROLOGUE  RQRTST
    !***PURPOSE  Quick check for RPQR79.
    !***LIBRARY   SLATEC
    !***TYPE      SINGLE PRECISION (RQRTST-S, CQRTST-C)
    !***AUTHOR  (UNKNOWN)
    !***ROUTINES CALLED  NUMXER, PASS, R1MACH, RPQR79, XERCLR, XGETF, XSETF
    !***REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Restructured using IF-THEN-ELSE-ENDIF, cleaned up FORMATs
    !           and changed TOL from sqrt R1MACH(3) to sqrt R1MACH(4) for
    !           the IBM 370 mainframes.  (RWC)
    !   911010  Code reworked and simplified.  (RWC and WRB)
    !***END PROLOGUE  RQRTST
    INTEGER itmp(7)
    COMPLEX root(7), chk(7)
    DIMENSION work(63)
    REAL coef(8)
    LOGICAL fatal
    !
    DATA chk/(1.4142135623731,1.4142135623731), &
      (1.4142135623731,-1.4142135623731), (0.0,2.0), (0.0,-2.0), &
      (-2.0,0.0), (-1.4142135623731,1.4142135623731), &
      (-1.4142135623731,-1.4142135623731)/
    !***FIRST EXECUTABLE STATEMENT  RQRTST
    IF ( Kprint>=2 ) WRITE (Lun,99001)
    !
    99001 FORMAT ('1',/,' RPQR79 QUICK CHECK')
    tol = SQRT(R1MACH(4))
    Ipass = 1
    !
    !     Initialize variables for testing.
    !
    beta = 0.0078125
    DO j = 1, 8
      coef(j) = beta
      beta = 2.0*beta
    ENDDO
    !
    CALL RPQR79(7,coef,root,ierr,work)
    !
    !     Check to see if test passed.
    !
    DO i = 1, 7
      itmp(i) = 0
    ENDDO
    !
    !     Check for roots in any order.
    !
    DO i = 1, 7
      DO j = 1, 7
        IF ( ABS(root(i)-chk(j))<=tol ) THEN
          itmp(j) = 1
          EXIT
        ENDIF
      ENDDO
    ENDDO
    !
    !     Check that we found all 7 roots.
    !
    Ipass = 1
    DO i = 1, 7
      Ipass = Ipass*itmp(i)
    ENDDO
    !
    !     Print test results.
    !
    IF ( Kprint>=3.OR.(Kprint>=2.AND.Ipass==0) ) THEN
      WRITE (Lun,99002)
      99002 FORMAT (/,' CHECK REAL AND IMAGINARY PARTS OF ROOT'/' COEFFICIENTS')
      WRITE (Lun,99003) (j,coef(j),j=1,8)
      99003 FORMAT (/(I6,3X,1P,E22.14))
      WRITE (Lun,99004)
      99004 FORMAT (//25X,'TABLE of ROOTS   ROOT         REAL  PART',12X,&
        'IMAG  PART'/'  NUMBER',8X,2(' of  ZERO ',12X))
      WRITE (Lun,99005) (j,root(j),j=1,7)
      99005 FORMAT (I6,3X,1P,2E22.14)
    ENDIF
    IF ( Kprint>=2 ) CALL PASS(Lun,1,Ipass)
    !
    !     Trigger 2 error conditions
    !
    CALL XGETF(kontrl)
    IF ( Kprint<=2 ) THEN
      CALL XSETF(0)
    ELSE
      CALL XSETF(1)
    ENDIF
    fatal = .FALSE.
    CALL XERCLR
    IF ( Kprint>=3 ) WRITE (Lun,99006)
    99006 FORMAT (//' TRIGGER 2 ERROR CONDITIONS'//)
    !
    !     CALL RPQR79 with 0 degree polynomial.
    !
    CALL RPQR79(0,coef,root,ierr,work)
    IF ( NUMXER(nerr)/=3 ) fatal = .TRUE.
    CALL XERCLR
    !
    !     CALL RPQR79 with zero leading coefficient.
    !
    coef(1) = 0.0
    CALL RPQR79(2,coef,root,ierr,work)
    IF ( NUMXER(nerr)/=2 ) fatal = .TRUE.
    CALL XERCLR
    !
    CALL XSETF(kontrl)
    IF ( fatal ) THEN
      Ipass = 0
      IF ( Kprint>=2 ) THEN
        WRITE (Lun,99007)
        99007 FORMAT (/' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
      ENDIF
    ELSEIF ( Kprint>=3 ) THEN
      WRITE (Lun,99008)
      99008 FORMAT (/' ALL INCORRECT ARGUMENT TESTS PASSED')
    ENDIF
    !
    IF ( Ipass==1.AND.Kprint>1 ) WRITE (Lun,99009)
    99009 FORMAT (/' **************RPQR79 PASSED ALL TESTS**************')
    IF ( Ipass==0.AND.Kprint/=0 ) WRITE (Lun,99010)
    99010 FORMAT (/' **************RPQR79 FAILED SOME TESTS*************')
    RETURN
  END SUBROUTINE RQRTST
  !DECK CQRTST
  SUBROUTINE CQRTST(Lun,Kprint,Ipass)
    IMPLICIT NONE
    INTEGER i, ierr, Ipass, j, kontrl, Kprint, Lun, nerr, NUMXER
    REAL R1MACH, tol
    !***BEGIN PROLOGUE  CQRTST
    !***PURPOSE  Quick check for CPQR79.
    !***LIBRARY   SLATEC
    !***TYPE      COMPLEX (RQRTST-S, CQRTST-C)
    !***AUTHOR  (UNKNOWN)
    !***ROUTINES CALLED  CPQR79, NUMXER, PASS, R1MACH, XERCLR, XGETF, XSETF
    !***REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Changed usage of R1MACH(3) to R1MACH(4).  (RWC)
    !   911010  Code reworked and simplified.  (RWC and WRB)
    !***END PROLOGUE  CQRTST
    INTEGER itest(2), itmp(7)
    REAL work(144)
    COMPLEX coeff1(9), coeff2(2), coeff3(2), root(8), chk1(8), chk2
    LOGICAL fatal
    !
    DATA coeff1/(1.0,0.0), (-7.0,-2.0), (8.0,6.0), (28.0,8.0), &
      (-49.0,-24.0), (7.0,2.0), (-8.0,-6.0), (-28.0,-8.0), (48.0,24.0)/
    DATA coeff2/(1.0,1.0), (1.0,3.0)/
    DATA coeff3/(0.0,0.0), (1.0,3.0)/
    DATA chk1/(4.0,2.0), (3.0,0.0), (-2.0,0.0), (2.0,0.0), (0.0,-1.0), &
      (-1.0,0.0), (0.0,1.0), (1.0,0.0)/
    DATA chk2/(-2.0,-1.0)/
    !***FIRST EXECUTABLE STATEMENT  CQRTST
    IF ( Kprint>=2 ) WRITE (Lun,99001)
    !
    99001 FORMAT ('1',/,' CPQR79 QUICK CHECK')
    tol = SQRT(R1MACH(4))
    Ipass = 1
    !
    !     First test.
    !
    CALL CPQR79(8,coeff1,root,ierr,work)
    !
    !     Check to see if test passed.
    !
    DO i = 1, 7
      itmp(i) = 0
    ENDDO
    !
    !     Check for roots in any order.
    !
    DO i = 1, 7
      DO j = 1, 7
        IF ( ABS(root(i)-chk1(j))<=tol ) THEN
          itmp(j) = 1
          EXIT
        ENDIF
      ENDDO
    ENDDO
    !
    !     Check that we found all 7 roots.
    !
    itest(1) = 1
    DO i = 1, 7
      itest(1) = itest(1)*itmp(i)
    ENDDO
    !
    !     Print test results.
    !
    IF ( Kprint>=3.OR.(Kprint>=2.AND.itest(1)==0) ) THEN
      WRITE (Lun,99008)
      WRITE (Lun,99009) (j,coeff1(j),j=1,9)
      WRITE (Lun,99010)
      WRITE (Lun,99011) (j,root(j),j=1,7)
    ENDIF
    IF ( Kprint>=2 ) CALL PASS(Lun,1,itest(1))
    !
    !     Set up next problem.
    !
    CALL CPQR79(1,coeff2,root,ierr,work)
    !
    !     Check to see if test passed.
    !
    itest(2) = 1
    IF ( ABS(root(1)-chk2)>tol ) itest(2) = 0
    !
    !     Print test results for second test.
    !
    IF ( Kprint>=3.OR.(Kprint>=2.AND.itest(1)==0) ) THEN
      WRITE (Lun,99002)
      99002 FORMAT (/,' TEST SUBSEQUENT RELATED CALL')
      WRITE (Lun,99008)
      WRITE (Lun,99009) (j,coeff2(j),j=1,2)
      WRITE (Lun,99010)
      WRITE (Lun,99011) (j,root(j),j=1,1)
    ENDIF
    IF ( Kprint>=2 ) CALL PASS(Lun,2,itest(2))
    !
    !     Trigger 2 error conditions
    !
    CALL XGETF(kontrl)
    IF ( Kprint<=2 ) THEN
      CALL XSETF(0)
    ELSE
      CALL XSETF(1)
    ENDIF
    fatal = .FALSE.
    CALL XERCLR
    IF ( Kprint>=3 ) WRITE (Lun,99003)
    99003 FORMAT (//' TRIGGER 2 ERROR CONDITIONS'//)
    !
    !     CALL CPQR79 with 0 degree polynomial.
    !
    CALL CPQR79(0,coeff2,root,ierr,work)
    IF ( NUMXER(nerr)/=3 ) fatal = .TRUE.
    CALL XERCLR
    !
    !     CALL CPQR79 with zero leading coefficient.
    !
    CALL CPQR79(2,coeff3,root,ierr,work)
    IF ( NUMXER(nerr)/=2 ) fatal = .TRUE.
    CALL XERCLR
    !
    CALL XSETF(kontrl)
    IF ( fatal ) THEN
      Ipass = 0
      IF ( Kprint>=2 ) THEN
        WRITE (Lun,99004)
        99004 FORMAT (/' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
      ENDIF
    ELSEIF ( Kprint>=3 ) THEN
      WRITE (Lun,99005)
      99005 FORMAT (/' ALL INCORRECT ARGUMENT TESTS PASSED')
    ENDIF
    !
    !     See if all tests passed.
    !
    Ipass = Ipass*itest(1)*itest(2)
    !
    IF ( Ipass==1.AND.Kprint>1 ) WRITE (Lun,99006)
    99006 FORMAT (/' **************CPQR79 PASSED ALL TESTS**************')
    IF ( Ipass==0.AND.Kprint/=0 ) WRITE (Lun,99007)
    99007 FORMAT (/' **************CPQR79 FAILED SOME TESTS*************')
    RETURN
    99008 FORMAT (/,' CHECK REAL AND IMAGINARY PARTS OF ROOT'/' COEFFICIENTS')
    99009 FORMAT (/(I6,3X,1P,2E22.14))
    99010 FORMAT (//25X,'TABLE of ROOTS   ROOT         REAL  PART',12X,&
      'IMAG  PART'/'  NUMBER',8X,2(' of  ZERO ',12X))
    99011 FORMAT (I6,3X,1P,2E22.14)
  END SUBROUTINE CQRTST
END MODULE TEST34_MOD
!DECK TEST34
PROGRAM TEST34
  USE TEST34_MOD
  IMPLICIT NONE
  INTEGER I1MACH
  !***BEGIN PROLOGUE  TEST34
  !***PURPOSE  Driver for testing SLATEC subprograms
  !***LIBRARY   SLATEC
  !***CATEGORY  F1A
  !***TYPE      ALL (TEST34-A)
  !***KEYWORDS  QUICK CHECK DRIVER
  !***AUTHOR  SLATEC Common Mathematical Library Committee
  !***DESCRIPTION
  !
  ! *Usage:
  !     One input data record is required
  !         READ (LIN, '(I1)') KPRINT
  !
  ! *Arguments:
  !     KPRINT = 0  Quick checks - No printing.
  !                 Driver       - Short pass or fail message printed.
  !              1  Quick checks - No message printed for passed tests,
  !                                short message printed for failed tests.
  !                 Driver       - Short pass or fail message printed.
  !              2  Quick checks - Print short message for passed tests,
  !                                fuller information for failed tests.
  !                 Driver       - Pass or fail message printed.
  !              3  Quick checks - Print complete quick check results.
  !                 Driver       - Pass or fail message printed.
  !
  ! *Description:
  !     Driver for testing SLATEC subprograms
  !        RPZERO   CPZERO
  !        FZERO    DFZERO
  !        RPQR79   CPQR79
  !
  !***REFERENCES  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***ROUTINES CALLED  CPRPQX, CQRTST, DFZTST, FZTEST, I1MACH, RQRTST,
  !                    XERMAX, XSETF, XSETUN
  !***REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  !***END PROLOGUE  TEST34
  INTEGER ipass, kprint, lin, lun, nfail
  !***FIRST EXECUTABLE STATEMENT  TEST34
  lun = I1MACH(2)
  lin = I1MACH(1)
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  CALL XERMAX(1000)
  CALL XSETUN(lun)
  IF ( kprint<=1 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  ENDIF
  !
  !     Test CPZERO and RPZERO
  !
  CALL CPRPQX(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test FZERO
  !
  CALL FZTEST(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test DFZERO
  !
  CALL DFZTST(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test RPQR79
  !
  CALL RQRTST(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test CPQR79
  !
  CALL CQRTST(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST34 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST34 *************')
  ENDIF
  STOP
END PROGRAM TEST34
