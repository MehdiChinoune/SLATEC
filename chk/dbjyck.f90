!*==DBJYCK.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DBJYCK
SUBROUTINE DBJYCK(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--DBJYCK5
  !*** Start of declarations inserted by SPAG
  INTEGER Kprint , NUMXER
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DBJYCK
  !***PURPOSE  Quick check for DBESJ and DBESY.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (BJYCK-S, DBJYCK-D)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Amos, D. E., (SNLA)
  !***DESCRIPTION
  !
  !   DBJYCK is a quick check routine for DBESJ and DBESY.  The main loops
  !   evaluate the Wronskian and test the error.  Underflow and overflow
  !   diagnostics are checked in addition to illegal arguments.
  !
  !***ROUTINES CALLED  D1MACH, DBESJ, DBESY, NUMXER, XERCLR, XGETF, XSETF
  !***REVISION HISTORY  (YYMMDD)
  !   750101  DATE WRITTEN
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891004  Removed unreachable code.  (WRB)
  !   891004  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901205  Changed usage of D1MACH(3) to D1MACH(4).  (RWC)
  !   910121  Editorial Changes.  (RWC)
  !   910501  Added TYPE record.  (WRB)
  !   910708  Code revised to test error returns for all values of
  !           KPRINT.  (WRB)
  !   910801  Editorial changes, some restructing and modifications to
  !           obtain more information when there is failure of the
  !           Wronskian.  (WRB)
  !***END PROLOGUE  DBJYCK
  INTEGER i , Ipass , ix , k , kontrl , Lun , m , n , nerr , nu , ny
  REAL(8) :: alp , del , er , fnu , fnup , rhpi , rx , tol , x
  REAL(8) :: fn(3) , w(5) , xx(5) , y(5)
  REAL(8) :: D1MACH
  LOGICAL fatal
  !***FIRST EXECUTABLE STATEMENT  DBJYCK
  IF ( Kprint>=2 ) WRITE (Lun,99001)
  99001 FORMAT (/' QUICK CHECKS FOR DBESJ AND DBESY'//)
  !
  Ipass = 1
  rhpi = 0.5D0/ATAN(1.0D0)
  xx(1) = 0.49D0
  xx(2) = 1.3D0
  xx(3) = 5.3D0
  xx(4) = 13.3D0
  xx(5) = 21.3D0
  fn(1) = 0.095D0
  fn(2) = 0.70D0
  fn(3) = 0.0D0
  tol = MAX(500.0D0*D1MACH(4),7.1D-12)
  DO m = 1 , 3
    DO n = 1 , 4
      DO nu = 1 , 4
        fnu = fn(m) + 12*(nu-1)
        DO ix = 1 , 5
          IF ( ix>=2.OR.nu<=3 ) THEN
            x = xx(ix)
            rx = rhpi/x
            CALL DBESJ(x,fnu,n,y,ny)
            IF ( ny==0 ) THEN
              CALL DBESY(x,fnu,n,w)
              fnup = fnu + n
              CALL DBESJ(x,fnup,1,y(n+1),ny)
              IF ( ny==0 ) THEN
                CALL DBESY(x,fnup,1,w(n+1))
                DO i = 1 , n
                  er = y(i+1)*w(i) - w(i+1)*y(i) - rx
                  er = ABS(er)/rx
                  IF ( er>tol ) THEN
                    Ipass = 0
                    IF ( Kprint>=2 ) WRITE (Lun,99002) m , n , nu , ix , i , &
                      x , er , tol , y(i) , y(i+1) , w(i) , w(i+1)
                    99002                   FORMAT (/' ERROR IN QUICK CHECK OF WRONSKIAN',&
                      1P/' M = ',I1,', N = ',I1,', NU = ',I1,&
                      ', IX = ',I1,', I = ',I1,/' X = ',E14.7,&
                      ', ER   = ',E14.7,', TOL = ',E14.7/' Y(I) = ',&
                      E14.7,', Y(I+1) = ',E14.7/' W(I) = ',E14.7,&
                      ', W(I+1) = ',E14.7)
                  ENDIF
                ENDDO
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !
  !     Check small values of X and order
  !
  n = 2
  fnu = 1.0D0
  x = D1MACH(4)/5.0D0
  rx = rhpi/x
  DO i = 1 , 3
    CALL DBESJ(x,fnu,n,y,ny)
    CALL DBESY(x,fnu,n,w)
    er = y(2)*w(1) - w(2)*y(1) - rx
    er = ABS(er)/rx
    IF ( er>tol ) THEN
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (Lun,99003) i , fnu , x , er , tol , y(i) , &
        y(i+1) , w(i) , w(i+1)
      99003     FORMAT (/' ERROR IN QUICK CHECK OF SMALL X AND ORDER',1P/' I = ',I1,&
        ',  FNU = ',E14.7/' X = ',E14.7,', ER = ',E14.7,', TOL = ',&
        E14.7/' Y(1) = ',E14.7,', Y(2) = ',E14.7/' W(1) = ',E14.7,&
        ', W(2) = ',E14.7)
      EXIT
    ENDIF
    fnu = D1MACH(4)/100.0D0
    x = xx(2*i-1)
    rx = rhpi/x
  ENDDO
  !
  !     Check large values of X and order
  !
  DO k = 1 , 2
    del = 30*(k-1)
    fnu = 70.0D0 + del
    DO n = 1 , 2
      x = 50.0D0 + del
      DO i = 1 , 5
        rx = rhpi/x
        CALL DBESJ(x,fnu,n,y,ny)
        IF ( ny==0 ) THEN
          CALL DBESY(x,fnu,n,w)
          IF ( n==1 ) THEN
            fnup = fnu + 1.0D0
            CALL DBESJ(x,fnup,1,y(2),ny)
            IF ( ny/=0 ) CYCLE
            CALL DBESY(x,fnup,1,w(2))
          ENDIF
          er = y(2)*w(1) - y(1)*w(2) - rx
          er = ABS(er)/rx
          IF ( er>tol ) THEN
            Ipass = 0
            IF ( Kprint>=2 ) WRITE (Lun,99004) k , n , i , x , er , tol , &
              y(1) , y(2) , w(1) , w(2)
            99004           FORMAT (/' ERROR IN QUICK CHECK OF LARGE X AND ORDER',&
              1P/' K = ',I1,', N = ',I1,', I = ',I1/' X = ',E14.7,&
              ', ER = ',E14.7,', TOL = ',E14.7/' Y(1) = ',E14.7,&
              ', Y(2) = ',E14.7/' W(1) = ',E14.7,', W(2) = ',E14.7)
            GOTO 100
          ENDIF
          x = x + 10.0D0
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  !
  !     Check underflow flags
  !
  100  x = D1MACH(1)*10.0D0
  alp = 12.3D0
  n = 3
  CALL DBESJ(x,alp,n,y,ny)
  IF ( ny/=3 ) THEN
    Ipass = 0
    IF ( Kprint>=2 ) WRITE (Lun,99005)
    99005   FORMAT (/' ERROR IN DBESJ UNDERFLOW TEST'/)
  ENDIF
  !
  !     Trigger 7 error conditions
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
  IF ( Kprint>=3 ) WRITE (Lun,99006)
  99006 FORMAT (//' TRIGGER 7 ERROR CONDITIONS'//)
  xx(1) = 1.0D0
  xx(2) = 1.0D0
  xx(3) = 1.0D0
  !
  !     Illegal arguments
  !
  DO i = 1 , 3
    xx(i) = -xx(i)
    n = INT(xx(3))
    CALL DBESJ(xx(1),xx(2),n,y,ny)
    IF ( NUMXER(nerr)/=2 ) THEN
      Ipass = 0
      fatal = .TRUE.
    ENDIF
    CALL XERCLR
    CALL DBESY(xx(1),xx(2),n,w)
    IF ( NUMXER(nerr)/=2 ) THEN
      Ipass = 0
      fatal = .TRUE.
    ENDIF
    CALL XERCLR
    xx(i) = -xx(i)
  ENDDO
  !
  !     Trigger overflow
  !
  x = D1MACH(1)*10.0D0
  n = 3
  alp = 2.3D0
  CALL DBESY(x,alp,n,w)
  IF ( NUMXER(nerr)/=6 ) THEN
    Ipass = 0
    fatal = .TRUE.
  ENDIF
  CALL XERCLR
  CALL XSETF(kontrl)
  IF ( fatal ) THEN
    IF ( Kprint>=2 ) THEN
      WRITE (Lun,99007)
      99007     FORMAT (/' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
    ENDIF
  ELSEIF ( Kprint>=3 ) THEN
    WRITE (Lun,99008)
    99008   FORMAT (/' ALL INCORRECT ARGUMENT TESTS PASSED')
  ENDIF
  !
  IF ( Ipass==1.AND.Kprint>=2 ) WRITE (Lun,99009)
  99009 FORMAT (/' *********DBESJ AND DBESY PASSED ALL TESTS*********')
  IF ( Ipass==0.AND.Kprint>=1 ) WRITE (Lun,99010)
  99010 FORMAT (/' *********DBESJ OR DBESY FAILED SOME TESTS*********')
  RETURN
END SUBROUTINE DBJYCK
