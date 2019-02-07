!*==RQRTST.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK RQRTST
SUBROUTINE RQRTST(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--RQRTST5
  !*** Start of declarations inserted by SPAG
  REAL beta, R1MACH, tol, work
  INTEGER i, ierr, Ipass, j, kontrl, Kprint, Lun, nerr, NUMXER
  !*** End of declarations inserted by SPAG
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
    99002   FORMAT (/,' CHECK REAL AND IMAGINARY PARTS OF ROOT'/' COEFFICIENTS')
    WRITE (Lun,99003) (j,coef(j),j=1,8)
    99003   FORMAT (/(I6,3X,1P,E22.14))
    WRITE (Lun,99004)
    99004   FORMAT (//25X,'TABLE of ROOTS'//'   ROOT         REAL  PART',12X,&
      'IMAG  PART'/'  NUMBER',8X,2(' of  ZERO ',12X))
    WRITE (Lun,99005) (j,root(j),j=1,7)
    99005   FORMAT (I6,3X,1P,2E22.14)
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
      99007     FORMAT (/' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
    ENDIF
  ELSEIF ( Kprint>=3 ) THEN
    WRITE (Lun,99008)
    99008   FORMAT (/' ALL INCORRECT ARGUMENT TESTS PASSED')
  ENDIF
  !
  IF ( Ipass==1.AND.Kprint>1 ) WRITE (Lun,99009)
  99009 FORMAT (/' **************RPQR79 PASSED ALL TESTS**************')
  IF ( Ipass==0.AND.Kprint/=0 ) WRITE (Lun,99010)
  99010 FORMAT (/' **************RPQR79 FAILED SOME TESTS*************')
  RETURN
END SUBROUTINE RQRTST
