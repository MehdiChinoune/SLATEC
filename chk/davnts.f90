!*==DAVNTS.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DAVNTS
SUBROUTINE DAVNTS(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--DAVNTS5
  !*** Start of declarations inserted by SPAG
  INTEGER kontrl
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DAVNTS
  !***PURPOSE  Quick check for DAVINT.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (AVNTST-S, DAVNTS-D)
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  D1MACH, DAVINT, XERCLR, XGETF, XSETF
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901205  Changed usage of D1MACH(3) to D1MACH(4).  (RWC)
  !   910501  Added PURPOSE and TYPE records.  (WRB)
  !   910708  Minor modifications in use of KPRINT.  (WRB)
  !   920210  Code restructured and revised to test error returns for all
  !           values of KPRINT.  (WRB)
  !***END PROLOGUE  DAVNTS
  REAL(8) :: D1MACH
  INTEGER i , ierr , Ipass , Kprint , Lun , n
  REAL(8) :: a , ans , b , del , rn1 , sqb , tol , tol1 , x(501) , &
    xint , y(501)
  LOGICAL fatal
  !***FIRST EXECUTABLE STATEMENT  DAVNTS
  IF ( Kprint>=2 ) WRITE (Lun,99001)
  99001 FORMAT ('1'/' DAVINT Quick Check')
  Ipass = 1
  tol = MAX(.0001D0,SQRT(D1MACH(4)))
  tol1 = 1.0D-2*tol
  !
  !     Perform first accuracy test.
  !
  a = 0.0D0
  b = 5.0D0
  xint = EXP(5.0D0) - 1.0D0
  n = 500
  rn1 = n - 1
  sqb = SQRT(b)
  del = 0.4D0*(b-a)/(n-1)
  DO i = 1 , n
    x(i) = sqb*SQRT(a+(i-1)*(b-a)/rn1) + del
    y(i) = EXP(x(i))
  ENDDO
  CALL DAVINT(x,y,n,a,b,ans,ierr)
  !
  !     See if test was passed.
  !
  IF ( ABS(ans-xint)>tol ) THEN
    Ipass = 0
    IF ( Kprint>=3 ) WRITE (Lun,99009) ierr , ans , xint
  ENDIF
  !
  !     Perform second accuracy test.
  !
  x(1) = 0.0D0
  x(2) = 5.0D0
  y(1) = 1.0D0
  y(2) = 0.5D0
  a = -0.5D0
  b = 0.5D0
  xint = 1.0D0
  CALL DAVINT(x,y,2,a,b,ans,ierr)
  !
  !     See if test was passed.
  !
  IF ( ABS(ans-xint)>tol1 ) THEN
    Ipass = 0
    IF ( Kprint>=3 ) WRITE (Lun,99009) ierr , ans , xint
  ENDIF
  !
  !     Send message indicating passage or failure of tests.
  !
  IF ( Kprint>=2 ) THEN
    IF ( Ipass==1 ) THEN
      IF ( Kprint>=3 ) WRITE (Lun,99002)
      99002     FORMAT (/' DAVINT passed both accuracy tests.')
    ELSE
      WRITE (Lun,99003)
      99003     FORMAT (/' DAVINT failed at least one accuracy test.')
    ENDIF
  ENDIF
  !
  !     Test error returns.
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
  IF ( Kprint>=3 ) THEN
    WRITE (Lun,99004)
    99004   FORMAT (/' Test error returns from DAVINT'/' 4 error messages expected'/&
      )
  ENDIF
  DO i = 1 , 20
    x(i) = (i-1)/19.0D0 - 0.01D0
    IF ( i/=1 ) y(i) = x(i)/(EXP(x(i))-1.0)
  ENDDO
  !
  !     Test IERR = 1 error return.
  !
  y(1) = 1.0D0
  CALL DAVINT(x,y,20,0.0D0,1.0D0,ans,ierr)
  IF ( ierr/=1 ) THEN
    Ipass = 0
    fatal = .TRUE.
    IF ( Kprint>=3 ) WRITE (Lun,99010) ierr , 1
  ENDIF
  CALL XERCLR
  !
  !     Test IERR = 2 error return.
  !
  CALL DAVINT(x,y,20,1.0D0,0.0D0,ans,ierr)
  IF ( ierr/=2 ) THEN
    Ipass = 0
    fatal = .TRUE.
    IF ( Kprint>=3 ) WRITE (Lun,99010) ierr , 2
  ENDIF
  IF ( ans/=0.0D0 ) THEN
    Ipass = 0
    fatal = .TRUE.
    IF ( Kprint>=3 ) WRITE (Lun,99011)
  ENDIF
  CALL XERCLR
  !
  !     Test IERR = 5 error return.
  !
  CALL DAVINT(x,y,1,0.0D0,1.0D0,ans,ierr)
  IF ( ierr/=5 ) THEN
    Ipass = 0
    fatal = .TRUE.
    IF ( Kprint>=3 ) WRITE (Lun,99010) ierr , 5
  ENDIF
  IF ( ans/=0.0D0 ) THEN
    Ipass = 0
    fatal = .TRUE.
    IF ( Kprint>=3 ) WRITE (Lun,99011)
  ENDIF
  CALL XERCLR
  !
  !     Test IERR = 4 error return.
  !
  x(1) = 1.0D0/19.0D0
  x(2) = 0.0D0
  CALL DAVINT(x,y,20,0.0D0,1.0D0,ans,ierr)
  IF ( ierr/=4 ) THEN
    Ipass = 0
    fatal = .TRUE.
    IF ( Kprint>=3 ) WRITE (Lun,99010) ierr , 4
  ENDIF
  IF ( ans/=0.0D0 ) THEN
    Ipass = 0
    fatal = .TRUE.
    IF ( Kprint>=3 ) WRITE (Lun,99011)
  ENDIF
  CALL XERCLR
  !
  !     Test IERR = 3 error return.
  !
  x(1) = 0.0D0
  x(2) = 1.0D0/19.0D0
  CALL DAVINT(x,y,20,0.0D0,.01D0,ans,ierr)
  IF ( ierr/=3 ) THEN
    Ipass = 0
    fatal = .TRUE.
    IF ( Kprint>=3 ) WRITE (Lun,99010) ierr , 3
  ENDIF
  IF ( ans/=0.0D0 ) THEN
    Ipass = 0
    fatal = .TRUE.
    IF ( Kprint>=3 ) WRITE (Lun,99011)
  ENDIF
  CALL XERCLR
  !
  !     Reset XERMSG control variables and write summary.
  !
  CALL XSETF(kontrl)
  IF ( fatal ) THEN
    IF ( Kprint>=2 ) THEN
      WRITE (Lun,99005)
      99005     FORMAT (/' At least one incorrect argument test FAILED')
    ENDIF
  ELSEIF ( Kprint>=3 ) THEN
    WRITE (Lun,99006)
    99006   FORMAT (/' All incorrect argument tests PASSED')
  ENDIF
  !
  !     Write PASS/FAIL message.
  !
  IF ( Ipass==1.AND.Kprint>=3 ) WRITE (Lun,99007)
  99007 FORMAT (/' ***************DAVINT PASSED ALL TESTS***************')
  IF ( Ipass==0.AND.Kprint>=2 ) WRITE (Lun,99008)
  99008 FORMAT (/' ***************DAVINT FAILED SOME TESTS**************')
  RETURN
  99009 FORMAT (/' FAILED ACCURACY TEST'/' IERR=',I2,5X,'COMPUTED ANS=',&
    E20.11/14X,'CORRECT ANS=',D20.11,5X,'REQUESTED ERR=',D10.2)
  99010 FORMAT (/' IERR =',I2,' and it should =',I2/)
  99011 FORMAT (1X,'ANS .NE. 0')
END SUBROUTINE DAVNTS
