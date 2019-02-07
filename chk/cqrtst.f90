!*==CQRTST.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CQRTST
SUBROUTINE CQRTST(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--CQRTST5
  !*** Start of declarations inserted by SPAG
  INTEGER i , ierr , Ipass , j , kontrl , Kprint , Lun , nerr , NUMXER
  REAL R1MACH , tol
  !*** End of declarations inserted by SPAG
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
  INTEGER itest(2) , itmp(7)
  REAL work(144)
  COMPLEX coeff1(9) , coeff2(2) , coeff3(2) , root(8) , chk1(8) , chk2
  LOGICAL fatal
  !
  DATA coeff1/(1.0,0.0) , (-7.0,-2.0) , (8.0,6.0) , (28.0,8.0) , &
    (-49.0,-24.0) , (7.0,2.0) , (-8.0,-6.0) , (-28.0,-8.0) , (48.0,24.0)/
  DATA coeff2/(1.0,1.0) , (1.0,3.0)/
  DATA coeff3/(0.0,0.0) , (1.0,3.0)/
  DATA chk1/(4.0,2.0) , (3.0,0.0) , (-2.0,0.0) , (2.0,0.0) , (0.0,-1.0) , &
    (-1.0,0.0) , (0.0,1.0) , (1.0,0.0)/
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
  DO i = 1 , 7
    itmp(i) = 0
  ENDDO
  !
  !     Check for roots in any order.
  !
  DO i = 1 , 7
    DO j = 1 , 7
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
  DO i = 1 , 7
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
    99002   FORMAT (/,' TEST SUBSEQUENT RELATED CALL')
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
      99004     FORMAT (/' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
    ENDIF
  ELSEIF ( Kprint>=3 ) THEN
    WRITE (Lun,99005)
    99005   FORMAT (/' ALL INCORRECT ARGUMENT TESTS PASSED')
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
  99010 FORMAT (//25X,'TABLE of ROOTS'//'   ROOT         REAL  PART',12X,&
    'IMAG  PART'/'  NUMBER',8X,2(' of  ZERO ',12X))
  99011 FORMAT (I6,3X,1P,2E22.14)
END SUBROUTINE CQRTST
