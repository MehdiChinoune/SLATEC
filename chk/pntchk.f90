!*==PNTCHK.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK PNTCHK
SUBROUTINE PNTCHK(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--PNTCHK5
  !***BEGIN PROLOGUE  PNTCHK
  !***PURPOSE  Quick check for POLINT, POLCOF and POLYVL
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (PNTCHK-S, DPNTCK-D)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  NUMXER, POLCOF, POLINT, POLYVL, R1MACH, XERCLR,
  !                    XGETF, XSETF
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901205  Changed usage of R1MACH(3) to R1MACH(4).  (RWC)
  !   910708  Minor modifications in use of KPRINT.  (WRB)
  !   920212  Code completely restructured to test errors for all values
  !           of KPRINT.  (WRB)
  !***END PROLOGUE  PNTCHK
  !     .. Scalar Arguments ..
  INTEGER Ipass , Kprint , Lun
  !     .. Local Scalars ..
  REAL tol , yf
  INTEGER i , ierr , kontrl , n , nerr
  LOGICAL fatal
  !     .. Local Arrays ..
  REAL c(6) , d(6) , dchk(6) , w(12) , x(6) , xchk(6) , y(6)
  !     .. External Functions ..
  REAL R1MACH
  INTEGER NUMXER
  EXTERNAL R1MACH , NUMXER
  !     .. External Subroutines ..
  EXTERNAL POLCOF , POLINT , POLYVL , XERCLR , XGETF , XSETF
  !     .. Intrinsic Functions ..
  INTRINSIC ABS , SQRT
  !     .. Data statements ..
  DATA x/1.0E0 , 2.0E0 , 3.0E0 , -1.0E0 , -2.0E0 , -3.0E0/
  DATA y/0.0E0 , 9.0E0 , 64.0E0 , 0.0E0 , 9.0E0 , 64.0E0/
  DATA xchk/1.0E0 , 0.0E0 , -2.0E0 , 0.0E0 , 1.0E0 , 0.0E0/
  DATA dchk/1.0E0 , 0.0E0 , -4.0E0 , 0.0E0 , 24.0E0 , 0.0E0/
  !***FIRST EXECUTABLE STATEMENT  PNTCHK
  IF ( Kprint>=2 ) WRITE (Lun,99001)
  !
  99001 FORMAT ('1'/' Test POLINT, POLCOF and POLYVL')
  !
  !     Initialize variables for tests.
  !
  tol = SQRT(R1MACH(4))
  Ipass = 1
  n = 6
  !
  !     Set up polynomial test.
  !
  CALL POLINT(n,x,y,c)
  CALL POLCOF(0.0E0,n,x,c,d,w)
  !
  !     Check to see if POLCOF test passed.
  !
  fatal = .FALSE.
  DO i = 1 , n
    IF ( ABS(d(i)-xchk(i))>tol ) THEN
      Ipass = 0
      fatal = .TRUE.
    ENDIF
  ENDDO
  IF ( fatal ) THEN
    IF ( Kprint>=2 ) WRITE (Lun,99007) 'FAILED' , (d(i),i=1,n)
  ELSE
    IF ( Kprint>=3 ) WRITE (Lun,99007) 'PASSED' , (d(i),i=1,n)
  ENDIF
  !
  !     Test POLYVL.
  !
  CALL POLYVL(5,0.0E0,yf,d,n,x,c,w,ierr)
  IF ( ABS(dchk(1)-yf)<=tol ) THEN
    IF ( Kprint>=3 ) WRITE (Lun,99008) 'PASSED' , yf , (d(i),i=1,5)
  ELSE
    Ipass = 0
    IF ( Kprint>=2 ) WRITE (Lun,99008) 'FAILED' , yf , (d(i),i=1,5)
  ENDIF
  !
  fatal = .FALSE.
  DO i = 1 , 5
    IF ( ABS(dchk(i+1)-d(i))>tol ) THEN
      Ipass = 0
      fatal = .TRUE.
    ENDIF
  ENDDO
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
  99002 FORMAT (/' 2 Error messages expected')
  CALL POLINT(0,x,y,c)
  IF ( NUMXER(nerr)/=2 ) THEN
    Ipass = 0
    fatal = .TRUE.
  ENDIF
  CALL XERCLR
  !
  x(1) = -1.0E0
  CALL POLINT(n,x,y,c)
  IF ( NUMXER(nerr)/=2 ) THEN
    Ipass = 0
    fatal = .TRUE.
  ENDIF
  CALL XERCLR
  !
  CALL XSETF(kontrl)
  IF ( fatal ) THEN
    IF ( Kprint>=2 ) THEN
      WRITE (Lun,99003)
      99003     FORMAT (/' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
    ENDIF
  ELSEIF ( Kprint>=3 ) THEN
    WRITE (Lun,99004)
    99004   FORMAT (/' ALL INCORRECT ARGUMENT TESTS PASSED')
  ENDIF
  !
  IF ( Ipass==1.AND.Kprint>=2 ) WRITE (Lun,99005)
  99005 FORMAT (/' ****************POLINT PASSED ALL TESTS**************')
  IF ( Ipass==0.AND.Kprint>=1 ) WRITE (Lun,99006)
  99006 FORMAT (/' ***************POLINT FAILED SOME TESTS**************')
  RETURN
  99007 FORMAT (/'POLCOF ',A,&
    ' test'/' Taylor coefficients for the quintic should be'/6X,&
    '1.000',5X,'0.000',4X,'-2.000',5X,'0.000',5X,'1.000',5X,&
    '0.000'/' Taylor coefficients from POLCOF are'/1X,6F10.3/)
  99008 FORMAT (' Derivative test ',&
    A/' The derivatives of the polynomial at zero as ',&
    'computed by POLYVL are'/1X,6F10.3/)
END SUBROUTINE PNTCHK
