MODULE TEST06_MOD
  IMPLICIT NONE

CONTAINS
  !** DBIKCK
  SUBROUTINE DBIKCK(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for DBESI and DBESK.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (BIKCK-S, DBIKCK-D)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Amos, D. E., (SNLA)
    !***
    ! **Description:**
    !
    !   DBIKCK is a quick check routine for DBESI and DBESK.  The main loops
    !   evaluate the Wronskian and test the error.  Underflow and overflow
    !   diagnostics are checked in addition to illegal arguments.
    !
    !***
    ! **Routines called:**  D1MACH, DBESI, DBESK, NUMXER, XERCLR, XGETF, XSETF

    !* REVISION HISTORY  (YYMMDD)
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
    
    INTEGER Kprint, NUMXER
    INTEGER i, Ipass, ix, k, kode, kontrl, Lun, m, n, nerr, nu, nw, ny
    REAL(8) :: alp, del, er, fnu, fnup, rx, tol, x
    REAL(8) :: fn(3), w(5), xx(5), y(5)
    REAL(8) :: D1MACH
    LOGICAL fatal
    !* FIRST EXECUTABLE STATEMENT  DBIKCK
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
    DO kode = 1, 2
      DO m = 1, 3
        DO n = 1, 4
          DO nu = 1, 4
            fnu = fn(m) + 12*(nu-1)
            DO ix = 1, 5
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
                        DO i = 1, n
                          er = y(i+1)*w(i) + w(i+1)*y(i) - rx
                          er = ABS(er)*x
                          IF ( er>tol ) THEN
                            Ipass = 0
                            IF ( Kprint>=2 ) WRITE (Lun,99002) kode, m, n, &
                              nu, ix, i, x, er, tol, y(i), y(i+1), &
                              w(i), w(i+1)
                            99002 FORMAT (/' ERROR IN QUICK CHECK OF WRONSKIAN',&
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
    DO i = 1, 3
      DO kode = 1, 2
        CALL DBESI(x,fnu,kode,n,y,ny)
        CALL DBESK(x,fnu,kode,n,w,nw)
        er = y(2)*w(1) + w(2)*y(1) - 1.0D0/x
        er = ABS(er)*x
        IF ( er>tol ) THEN
          Ipass = 0
          IF ( Kprint>=2 ) WRITE (Lun,99003) i, kode, fnu, x, er, tol, &
            y(1), y(2), w(1), w(2)
          99003 FORMAT (/' ERROR IN QUICK CHECK OF SMALL X AND ORDER',1P/' I = ',I1,&
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
    DO k = 1, 2
      del = 30*(k-1)
      fnu = 45.0D0 + del
      DO n = 1, 2
        x = 20.0D0 + del
        DO i = 1, 5
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
                IF ( Kprint>=2 ) WRITE (Lun,99004) k, n, i, fnup, x, er, &
                  tol, y(1), y(2), w(1), w(2)
                99004 FORMAT (/' ERROR IN QUICK CHECK OF LARGE X AND ORDER',&
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
      99005 FORMAT (/' ERROR IN DBESI UNDERFLOW TEST'/)
    ENDIF
    !
    x = LOG(D1MACH(2)/10.0D0) + 20.0D0
    alp = 1.3D0
    n = 3
    CALL DBESK(x,alp,1,n,w,nw)
    IF ( nw/=3 ) THEN
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (Lun,99006)
      99006 FORMAT (/' ERROR IN DBESK UNDERFLOW TEST'/)
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
    DO i = 1, 4
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
        99008 FORMAT (/' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
      ENDIF
    ELSEIF ( Kprint>=3 ) THEN
      WRITE (Lun,99009)
      99009 FORMAT (/' ALL INCORRECT ARGUMENT TESTS PASSED')
    ENDIF
    !
    IF ( Ipass==1.AND.Kprint>=2 ) WRITE (Lun,99010)
    99010 FORMAT (/' *********DBESI AND DBESK PASSED ALL TESTS***********')
    IF ( Ipass==0.AND.Kprint>=1 ) WRITE (Lun,99011)
    99011 FORMAT (/' *********DBESI OR DBESK FAILED SOME TESTS***********')
    RETURN
  END SUBROUTINE DBIKCK
  !** DBJYCK
  SUBROUTINE DBJYCK(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for DBESJ and DBESY.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (BJYCK-S, DBJYCK-D)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Amos, D. E., (SNLA)
    !***
    ! **Description:**
    !
    !   DBJYCK is a quick check routine for DBESJ and DBESY.  The main loops
    !   evaluate the Wronskian and test the error.  Underflow and overflow
    !   diagnostics are checked in addition to illegal arguments.
    !
    !***
    ! **Routines called:**  D1MACH, DBESJ, DBESY, NUMXER, XERCLR, XGETF, XSETF

    !* REVISION HISTORY  (YYMMDD)
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
    
    INTEGER Kprint, NUMXER
    INTEGER i, Ipass, ix, k, kontrl, Lun, m, n, nerr, nu, ny
    REAL(8) :: alp, del, er, fnu, fnup, rhpi, rx, tol, x
    REAL(8) :: fn(3), w(5), xx(5), y(5)
    REAL(8) :: D1MACH
    LOGICAL fatal
    !* FIRST EXECUTABLE STATEMENT  DBJYCK
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
    DO m = 1, 3
      DO n = 1, 4
        DO nu = 1, 4
          fnu = fn(m) + 12*(nu-1)
          DO ix = 1, 5
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
                  DO i = 1, n
                    er = y(i+1)*w(i) - w(i+1)*y(i) - rx
                    er = ABS(er)/rx
                    IF ( er>tol ) THEN
                      Ipass = 0
                      IF ( Kprint>=2 ) WRITE (Lun,99002) m, n, nu, ix, i, &
                        x, er, tol, y(i), y(i+1), w(i), w(i+1)
                      99002 FORMAT (/' ERROR IN QUICK CHECK OF WRONSKIAN',&
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
    DO i = 1, 3
      CALL DBESJ(x,fnu,n,y,ny)
      CALL DBESY(x,fnu,n,w)
      er = y(2)*w(1) - w(2)*y(1) - rx
      er = ABS(er)/rx
      IF ( er>tol ) THEN
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Lun,99003) i, fnu, x, er, tol, y(i), &
          y(i+1), w(i), w(i+1)
        99003 FORMAT (/' ERROR IN QUICK CHECK OF SMALL X AND ORDER',1P/' I = ',I1,&
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
    DO k = 1, 2
      del = 30*(k-1)
      fnu = 70.0D0 + del
      DO n = 1, 2
        x = 50.0D0 + del
        DO i = 1, 5
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
              IF ( Kprint>=2 ) WRITE (Lun,99004) k, n, i, x, er, tol, &
                y(1), y(2), w(1), w(2)
              99004 FORMAT (/' ERROR IN QUICK CHECK OF LARGE X AND ORDER',&
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
      99005 FORMAT (/' ERROR IN DBESJ UNDERFLOW TEST'/)
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
    DO i = 1, 3
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
        99007 FORMAT (/' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
      ENDIF
    ELSEIF ( Kprint>=3 ) THEN
      WRITE (Lun,99008)
      99008 FORMAT (/' ALL INCORRECT ARGUMENT TESTS PASSED')
    ENDIF
    !
    IF ( Ipass==1.AND.Kprint>=2 ) WRITE (Lun,99009)
    99009 FORMAT (/' *********DBESJ AND DBESY PASSED ALL TESTS*********')
    IF ( Ipass==0.AND.Kprint>=1 ) WRITE (Lun,99010)
    99010 FORMAT (/' *********DBESJ OR DBESY FAILED SOME TESTS*********')
    RETURN
  END SUBROUTINE DBJYCK
  !** DEG8CK
  SUBROUTINE DEG8CK(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for DEXINT and DGAUS8.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (EG8CK-S, DEG8CK-D)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Amos, D. E., (SNLA)
    !***
    ! **Description:**
    !
    !   DEG8CK is a quick check routine for DEXINT and DGAUS8.  Exponential
    !   integrals from DEXINT are checked against quadratures from DGAUS8.
    !
    !***
    ! **Routines called:**  D1MACH, DEXINT, DFEIN, DGAUS8
    !***
    ! COMMON BLOCKS    DFEINX

    !* REVISION HISTORY  (YYMMDD)
    !   800501  DATE WRITTEN
    !   890718  Added check when testing error conditions.  (WRB)
    !   890718  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   910708  Code revised to test error returns for all values of
    !           KPRINT.  (WRB)
    !   920206  Corrected argument list in CALL to DEXINT.  (WRB)
    
    INTEGER Kprint
    COMMON /DFEINX/ X, A, FKM
    INTEGER i, icase, ie, ierr, ii, ik, Ipass, ix, iy, k, ke, kk, &
      kode, kx, Lun, m, n, nm, nz
    REAL(8) :: A, ans, atol, bb, en, er, ex, FKM, sig, sum, &
      tol, t1, t2, X, xx, y
    REAL(8) :: D1MACH
    DIMENSION en(4), y(4), xx(5)
    LOGICAL fatal
    !* FIRST EXECUTABLE STATEMENT  DEG8CK
    IF ( Kprint>=2 ) WRITE (Lun,99001)
    !
    99001 FORMAT ('1'/' QUICK CHECK FOR DEXINT AND DGAUS8'/)
    Ipass = 1
    tol = SQRT(MAX(D1MACH(4),1.0D-18))
    DO kode = 1, 2
      ik = kode - 1
      FKM = ik
      DO n = 1, 25, 8
        DO m = 1, 4
          nm = n + m - 1
          DO ix = 1, 25, 8
            X = ix - 0.20D0
            CALL DEXINT(X,n,kode,m,tol,en,nz,ierr)
            kx = INT( X + 0.5D0 )
            IF ( kx==0 ) kx = 1
            icase = 1
            A = n
            IF ( kx>n ) THEN
              icase = 2
              A = nm
              IF ( kx<nm ) THEN
                icase = 3
                A = kx
              ENDIF
            ENDIF
            sig = 3.0D0/X
            t2 = 1.0D0
            sum = 0.0D0
            DO
              t1 = t2
              t2 = t2 + sig
              atol = tol
              CALL DGAUS8(DFEIN,t1,t2,atol,ans,ierr)
              sum = sum + ans
              IF ( ABS(ans)<ABS(sum)*tol ) THEN
                ex = 1.0D0
                IF ( kode==1 ) ex = EXP(-X)
                bb = A
                IF ( icase==3 ) THEN
                  iy = kx - n + 1
                  y(iy) = sum
                  ke = m - iy
                  ie = iy - 1
                  kk = iy
                  ii = iy
                ELSEIF ( icase/=2 ) THEN
                  y(1) = sum
                  IF ( m==1 ) GOTO 5
                  ke = m - 1
                  kk = 1
                ELSE
                  y(m) = sum
                  IF ( m==1 ) GOTO 5
                  ie = m - 1
                  ii = m
                  EXIT
                ENDIF
                !
                !             Forward recur
                !
                DO k = 1, ke
                  y(kk+1) = (ex-X*y(kk))/bb
                  bb = bb + 1.0D0
                  kk = kk + 1
                ENDDO
                IF ( icase==3 ) EXIT
                GOTO 5
              ENDIF
            ENDDO
            bb = A - 1.0D0
            !
            !             Backward recur
            !
            DO i = 1, ie
              y(ii-1) = (ex-bb*y(ii))/X
              bb = bb - 1.0D0
              ii = ii - 1
            ENDDO
            5 CONTINUE
            DO i = 1, m
              er = ABS((y(i)-en(i))/y(i))
              IF ( er>tol ) THEN
                WRITE (Lun,99002)
                99002 FORMAT (//' ERROR IN DEG8CK COMPARISON TEST'/)
                Ipass = 0
                GOTO 100
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    !     Trigger 6 error conditions.
    !
    100  fatal = .FALSE.
    !
    IF ( Kprint>=3 ) WRITE (Lun,99003)
    99003 FORMAT (/' TRIGGER 6 ERROR CONDITIONS'/)
    xx(1) = 1.0D0
    xx(2) = 1.0D0
    xx(3) = 1.0D0
    xx(4) = 1.0D0
    xx(5) = 0.01D0
    DO i = 1, 5
      xx(i) = -xx(i)
      k = INT( xx(2) )
      n = INT( xx(3) )
      m = INT( xx(4) )
      CALL DEXINT(xx(i),n,k,m,xx(5),en,nz,ierr)
      IF ( ierr/=1 ) THEN
        Ipass = 0
        fatal = .TRUE.
        WRITE (Lun,99004) i
        99004 FORMAT (' Error occurred with DO index I =',I2)
      ENDIF
      xx(i) = -xx(i)
    ENDDO
    X = 0.0D0
    tol = 1.0D-2
    CALL DEXINT(X,1,1,1,tol,en,nz,ierr)
    IF ( ierr/=1 ) THEN
      Ipass = 0
      fatal = .TRUE.
      WRITE (Lun,99005)
      99005 FORMAT (' Error occurred with X = 0.0')
    ENDIF
    IF ( fatal ) THEN
      IF ( Kprint>=2 ) THEN
        WRITE (Lun,99006)
        99006 FORMAT (/' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
      ENDIF
    ELSEIF ( Kprint>=3 ) THEN
      WRITE (Lun,99007)
      99007 FORMAT (/' ALL INCORRECT ARGUMENT TESTS PASSED')
    ENDIF
    !
    IF ( Ipass==1.AND.Kprint>=2 ) WRITE (Lun,99008)
    99008 FORMAT (/' *********DEXINT AND DGAUS8 PASSED ALL TESTS*********')
    IF ( Ipass==0.AND.Kprint>=1 ) WRITE (Lun,99009)
    99009 FORMAT (/' *********DEXINT OR DGAUS8 FAILED SOME TESTS*********')
    RETURN
  END SUBROUTINE DEG8CK
  !** DFEIN
  REAL(8) FUNCTION DFEIN(T)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to DEG8CK.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)
    !***
    ! COMMON BLOCKS    DFEINX

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    COMMON /DFEINX/ X, A, FKM
    REAL(8) :: X, A, FKM, T, aln
    !* FIRST EXECUTABLE STATEMENT  DFEIN
    aln = (FKM-T)*X - A*LOG(T)
    DFEIN = EXP(aln)
  END FUNCTION DFEIN
END MODULE TEST06_MOD
!** TEST06
PROGRAM TEST06
  USE TEST06_MOD
  IMPLICIT NONE
  !>
  !***
  !  Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C
  !***
  ! **Type:**      DOUBLE PRECISION (TEST05-S, TEST06-D)
  !***
  ! **Keywords:**  QUICK CHECK DRIVER
  !***
  ! **Author:**  SLATEC Common Mathematical Library Committee
  !***
  ! **Description:**
  !
  !- Usage:
  !     One input data record is required
  !         READ (LIN, '(I1)') KPRINT
  !
  !- Arguments:
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
  !- Description:
  !     Driver for testing SLATEC subprograms
  !        DEXINT   DGAUS8
  !        DBESI    DBESK
  !        DBESJ    DBESY
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  DBIKCK, DBJYCK, DEG8CK, I1MACH, XERMAX, XSETF,
  !                    XSETUN

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  
  INTEGER I1MACH
  INTEGER ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST06
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
  !     Test DEXINT and DQAUS8
  !
  CALL DEG8CK(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test DBESI and DBESK
  !
  CALL DBIKCK(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test DBESJ and DBESY
  !
  CALL DBJYCK(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST06 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST06  *************')
  ENDIF
  STOP
END PROGRAM TEST06
