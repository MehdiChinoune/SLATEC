!*==DPNTCK.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DPNTCK
      SUBROUTINE DPNTCK(Lun,Kprint,Ipass)
      IMPLICIT NONE
!*--DPNTCK5
!***BEGIN PROLOGUE  DPNTCK
!***PURPOSE  Quick check for DPLINT, DPOLCF and DPOLVL
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (PNTCHK-S, DPNTCK-D)
!***KEYWORDS  QUICK CHECK
!***AUTHOR  Boland, W. Robert, (LANL)
!***ROUTINES CALLED  D1MACH, DPLINT, DPOLCF, DPOLVL, NUMXER, XERCLR,
!                    XGETF, XSETF
!***REVISION HISTORY  (YYMMDD)
!   920212  DATE WRITTEN
!***END PROLOGUE  DPNTCK
!     .. Scalar Arguments ..
      INTEGER Ipass , Kprint , Lun
!     .. Local Scalars ..
      DOUBLE PRECISION tol , yf
      INTEGER i , ierr , kontrl , n , nerr
      LOGICAL fatal
!     .. Local Arrays ..
      DOUBLE PRECISION c(6) , d(6) , dchk(6) , w(12) , x(6) , xchk(6) , y(6)
!     .. External Functions ..
      DOUBLE PRECISION D1MACH
      INTEGER NUMXER
      EXTERNAL D1MACH , NUMXER
!     .. External Subroutines ..
      EXTERNAL DPOLCF , DPLINT , DPOLVL , XERCLR , XGETF , XSETF
!     .. Intrinsic Functions ..
      INTRINSIC ABS , SQRT
!     .. Data statements ..
      DATA x/1.0D0 , 2.0D0 , 3.0D0 , -1.0D0 , -2.0D0 , -3.0D0/
      DATA y/0.0D0 , 9.0D0 , 64.0D0 , 0.0D0 , 9.0D0 , 64.0D0/
      DATA xchk/1.0D0 , 0.0D0 , -2.0D0 , 0.0D0 , 1.0D0 , 0.0D0/
      DATA dchk/1.0D0 , 0.0D0 , -4.0D0 , 0.0D0 , 24.0D0 , 0.0D0/
!***FIRST EXECUTABLE STATEMENT  DPNTCK
      IF ( Kprint>=2 ) WRITE (Lun,99001)
!
99001 FORMAT ('1'/' Test DPLINT, DPOLCF and DPOLVL')
!
!     Initialize variables for tests.
!
      tol = SQRT(D1MACH(4))
      Ipass = 1
      n = 6
!
!     Set up polynomial test.
!
      CALL DPLINT(n,x,y,c)
      CALL DPOLCF(0.0D0,n,x,c,d,w)
!
!     Check to see if DPOLCF test passed.
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
!     Test DPOLVL.
!
      CALL DPOLVL(5,0.0D0,yf,d,n,x,c,w,ierr)
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
      CALL DPLINT(0,x,y,c)
      IF ( NUMXER(nerr)/=2 ) THEN
        Ipass = 0
        fatal = .TRUE.
      ENDIF
      CALL XERCLR
!
      x(1) = -1.0D0
      CALL DPLINT(n,x,y,c)
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
99005 FORMAT (/' ****************DPLINT PASSED ALL TESTS**************')
      IF ( Ipass==0.AND.Kprint>=1 ) WRITE (Lun,99006)
99006 FORMAT (/' ***************DPLINT FAILED SOME TESTS**************')
      RETURN
99007 FORMAT (/'DPOLCF ',A,
     &        ' test'/' Taylor coefficients for the quintic should be'/6X,
     &        '1.000',5X,'0.000',4X,'-2.000',5X,'0.000',5X,'1.000',5X,
     &        '0.000'/' Taylor coefficients from DPOLCF are'/1X,6F10.3/)
99008 FORMAT (' Derivative test ',
     &        A/' The derivatives of the polynomial at zero as ',
     &        'computed by DPOLVL are'/1X,6F10.3/)
      END SUBROUTINE DPNTCK
