!*==DBIKCK.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DBIKCK
SUBROUTINE DBIKCK(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--DBIKCK5
  !*** Start of declarations inserted by SPAG
  INTEGER Kprint , NUMXER
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DBIKCK
  !***PURPOSE  Quick check for DBESI and DBESK.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (BIKCK-S, DBIKCK-D)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Amos, D. E., (SNLA)
  !***DESCRIPTION
  !
  !   DBIKCK is a quick check routine for DBESI and DBESK.  The main loops
  !   evaluate the Wronskian and test the error.  Underflow and overflow
  !   diagnostics are checked in addition to illegal arguments.
  !
  !***ROUTINES CALLED  D1MACH, DBESI, DBESK, NUMXER, XERCLR, XGETF, XSETF
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
  !***END PROLOGUE  DBIKCK
  INTEGER i , Ipass , ix , k , kode , kontrl , Lun , m , n , nerr , nu , &
    nw , ny
  DOUBLE PRECISION alp , del , er , fnu , fnup , rx , tol , x
  DOUBLE PRECISION fn(3) , w(5) , xx(5) , y(5)
  DOUBLE PRECISION D1MACH
  LOGICAL fatal
  !***FIRST EXECUTABLE STATEMENT  DBIKCK
  IF ( Kprint>=2 ) WRITE (Lun,99001)
  99001 FORMAT (/' QUICK CHECKS FOR DBESI AND DBESK'//)
  !
  Ipass = 1
  xx(1) = 0.49D0
  xx(2) = 1.3D0
  xx(3) = 5.3D0
  xx(4) = 13.3D0
  xx(5) = 21.3D0
  fn(1) = 0.095D0
  fn(2) = 0.70D0
  fn(3) = 0.0D0
  tol = MAX(500.0D0*D1MACH(4),7.1D-12)
  DO kode = 1 , 2
    DO m = 1 , 3
      DO n = 1 , 4
        DO nu = 1 , 4
          fnu = fn(m) + 12*(nu-1)
          DO ix = 1 , 5
            IF ( ix>=2.OR.nu<=3 ) THEN
              x = xx(ix)
              rx = 1.0D0/x
              CALL DBESI(x,fnu,kode,n,y,ny)
              IF ( ny==0 ) THEN
                CALL DBESK(x,fnu,kode,n,w,nw)
                IF ( nw==0 ) THEN
                  fnup = fnu + n
                  CALL DBESI(x,fnup,kode,1,y(n+1),ny)
                  IF ( ny==0 ) THEN
                    CALL DBESK(x,fnup,kode,1,w(n+1),nw)
                    IF ( nw==0 ) THEN
                      DO i = 1 , n
                        er = y(i+1)*w(i) + w(i+1)*y(i) - rx
                        er = ABS(er)*x
                        IF ( er>tol ) THEN
                          Ipass = 0
                          IF ( Kprint>=2 ) WRITE (Lun,99002) kode , m , n , &
                            nu , ix , i , x , er , tol , y(i) , y(i+1) , &
                            w(i) , w(i+1)
                          99002                         FORMAT (/' ERROR IN QUICK CHECK OF WRONSKIAN',&
                            1P/' KODE = ',I1,', M = ',I1,', N = ',I1,&
                            ', NU = ',I1,', IX = ',I1,', I = ',&
                            I1/' X = ',E14.7,', ER   = ',E14.7,&
                            ', TOL = ',E14.7/' Y(I) = ',E14.7,&
                            ', Y(I+1) = ',E14.7/' W(I) = ',E14.7,&
                            ', W(I+1) = ',E14.7)
                        ENDIF
                      ENDDO
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !
  !     Check small values of X and order
  !
  n = 2
  fnu = 1.0D0
  x = D1MACH(4)
  DO i = 1 , 3
    DO kode = 1 , 2
      CALL DBESI(x,fnu,kode,n,y,ny)
      CALL DBESK(x,fnu,kode,n,w,nw)
      er = y(2)*w(1) + w(2)*y(1) - 1.0D0/x
      er = ABS(er)*x
      IF ( er>tol ) THEN
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Lun,99003) i , kode , fnu , x , er , tol , &
          y(1) , y(2) , w(1) , w(2)
        99003       FORMAT (/' ERROR IN QUICK CHECK OF SMALL X AND ORDER',1P/' I = ',I1,&
          ', KODE = ',I1,', FNU = ',E14.7/' X = ',E14.7,', ER = ',&
          E14.7,', TOL = ',E14.7/' Y(1) = ',E14.7,', Y(2) = ',&
          E14.7/' W(1) = ',E14.7,', W(2) = ',E14.7)
        EXIT
      ENDIF
    ENDDO
    !
    fnu = D1MACH(4)/100.0D0
    x = xx(2*i-1)
  ENDDO
  !
  !     Check large values of X and order
  !
  kode = 2
  DO k = 1 , 2
    del = 30*(k-1)
    fnu = 45.0D0 + del
    DO n = 1 , 2
      x = 20.0D0 + del
      DO i = 1 , 5
        rx = 1.0D0/x
        CALL DBESI(x,fnu,kode,n,y,ny)
        IF ( ny==0 ) THEN
          CALL DBESK(x,fnu,kode,n,w,nw)
          IF ( nw==0 ) THEN
            IF ( n==1 ) THEN
              fnup = fnu + 1.0D0
              CALL DBESI(x,fnup,kode,1,y(2),ny)
              IF ( ny/=0 ) CYCLE
              CALL DBESK(x,fnup,kode,1,w(2),nw)
              IF ( nw/=0 ) CYCLE
            ENDIF
            er = y(2)*w(1) + y(1)*w(2) - rx
            er = ABS(er)*x
            IF ( er>tol ) THEN
              Ipass = 0
              IF ( Kprint>=2 ) WRITE (Lun,99004) k , n , i , fnup , x , er , &
                tol , y(1) , y(2) , w(1) , w(2)
              99004             FORMAT (/' ERROR IN QUICK CHECK OF LARGE X AND ORDER',&
                1P/' K = ',I1,', N = ',I1,', I = ',I1,', FNUP = ',&
                E14.7/' X = ',E14.7,', ER = ',E14.7,', TOL = ',&
                E14.7/' Y(1) = ',E14.7,', Y(2) = ',E14.7/' W(1) = ',&
                E14.7,', W(2) = ',E14.7)
              GOTO 100
            ENDIF
            x = x + 10.0D0
          ENDIF
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
  CALL DBESI(x,alp,1,n,y,ny)
  IF ( ny/=3 ) THEN
    Ipass = 0
    IF ( Kprint>=2 ) WRITE (Lun,99005)
    99005   FORMAT (/' ERROR IN DBESI UNDERFLOW TEST'/)
  ENDIF
  !
  x = LOG(D1MACH(2)/10.0D0) + 20.0D0
  alp = 1.3D0
  n = 3
  CALL DBESK(x,alp,1,n,w,nw)
  IF ( nw/=3 ) THEN
    Ipass = 0
    IF ( Kprint>=2 ) WRITE (Lun,99006)
    99006   FORMAT (/' ERROR IN DBESK UNDERFLOW TEST'/)
  ENDIF
  !
  !     Trigger 10 error conditions
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
  IF ( Kprint>=3 ) WRITE (Lun,99007)
  99007 FORMAT (//' TRIGGER 10 ERROR CONDITIONS'//)
  xx(1) = 1.0D0
  xx(2) = 1.0D0
  xx(3) = 1.0D0
  xx(4) = 1.0D0
  !
  !     Illegal arguments
  !
  DO i = 1 , 4
    xx(i) = -xx(i)
    k = INT(xx(3))
    n = INT(xx(4))
    CALL DBESI(xx(1),xx(2),k,n,y,ny)
    IF ( NUMXER(nerr)/=2 ) THEN
      Ipass = 0
      fatal = .TRUE.
    ENDIF
    CALL XERCLR
    CALL DBESK(xx(1),xx(2),k,n,w,nw)
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
  x = LOG(D1MACH(2)/10.0D0) + 20.0D0
  n = 3
  alp = 2.3D0
  CALL DBESI(x,alp,1,n,y,ny)
  IF ( NUMXER(nerr)/=6 ) THEN
    Ipass = 0
    fatal = .TRUE.
  ENDIF
  CALL XERCLR
  !
  x = D1MACH(1)*10.0D0
  CALL DBESK(x,alp,1,n,w,nw)
  IF ( NUMXER(nerr)/=6 ) THEN
    Ipass = 0
    fatal = .TRUE.
  ENDIF
  CALL XERCLR
  !
  CALL XSETF(kontrl)
  IF ( fatal ) THEN
    IF ( Kprint>=2 ) THEN
      WRITE (Lun,99008)
      99008     FORMAT (/' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
    ENDIF
  ELSEIF ( Kprint>=3 ) THEN
    WRITE (Lun,99009)
    99009   FORMAT (/' ALL INCORRECT ARGUMENT TESTS PASSED')
  ENDIF
  !
  IF ( Ipass==1.AND.Kprint>=2 ) WRITE (Lun,99010)
  99010 FORMAT (/' *********DBESI AND DBESK PASSED ALL TESTS***********')
  IF ( Ipass==0.AND.Kprint>=1 ) WRITE (Lun,99011)
  99011 FORMAT (/' *********DBESI OR DBESK FAILED SOME TESTS***********')
  RETURN
END SUBROUTINE DBIKCK
