!** CUOIK
PURE SUBROUTINE CUOIK(Z,Fnu,Kode,Ikflg,N,Y,Nuf,Tol,Elim,Alim)
  !> Subsidiary to CBESH, CBESI and CBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CUOIK-A, ZUOIK-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     CUOIK COMPUTES THE LEADING TERMS OF THE UNIFORM ASYMPTOTIC
  !     EXPANSIONS FOR THE I AND K FUNCTIONS AND COMPARES THEM
  !     (IN LOGARITHMIC FORM) TO ALIM AND ELIM FOR OVER AND UNDERFLOW
  !     WHERE ALIM<ELIM. IF THE MAGNITUDE, BASED ON THE LEADING
  !     EXPONENTIAL, IS LESS THAN ALIM OR GREATER THAN -ALIM, THEN
  !     THE RESULT IS ON SCALE. IF NOT, THEN A REFINED TEST USING OTHER
  !     MULTIPLIERS (IN LOGARITHMIC FORM) IS MADE BASED ON ELIM. HERE
  !     EXP(-ELIM)=SMALLEST MACHINE NUMBER*1.0E+3 AND EXP(-ALIM)=
  !     EXP(-ELIM)/TOL
  !
  !     IKFLG=1 MEANS THE I SEQUENCE IS TESTED
  !          =2 MEANS THE K SEQUENCE IS TESTED
  !     NUF = 0 MEANS THE LAST MEMBER OF THE SEQUENCE IS ON SCALE
  !         =-1 MEANS AN OVERFLOW WOULD OCCUR
  !     IKFLG=1 AND NUF>0 MEANS THE LAST NUF Y VALUES WERE SET TO ZERO
  !             THE FIRST N-NUF VALUES MUST BE SET BY ANOTHER ROUTINE
  !     IKFLG=2 AND NUF=N MEANS ALL Y VALUES WERE SET TO ZERO
  !     IKFLG=2 AND 0<NUF<N NOT CONSIDERED. Y MUST BE SET BY
  !             ANOTHER ROUTINE
  !
  !***
  ! **See also:**  CBESH, CBESI, CBESK
  !***
  ! **Routines called:**  CUCHK, CUNHJ, CUNIK, R1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : tiny_sp
  !
  INTEGER, INTENT(IN) :: Ikflg, Kode, N
  INTEGER, INTENT(OUT) :: Nuf
  REAL(SP), INTENT(IN) :: Alim, Elim, Fnu, Tol
  COMPLEX(SP), INTENT(IN) :: Z
  COMPLEX(SP), INTENT(OUT) :: Y(N)
  !
  INTEGER :: i, iform, init, nn, nw
  COMPLEX(SP) :: arg, asum, bsum, cwrk(16), cz, phi, summ, zb, zeta1, zeta2, zn, zr
  REAL(SP) :: aarg, aphi, ascle, ax, ay, fnn, gnn, gnu, rcz, x, yy
  REAL(SP), PARAMETER :: aic = 1.265512123484645396_SP
  !* FIRST EXECUTABLE STATEMENT  CUOIK
  Nuf = 0
  nn = N
  x = REAL(Z)
  zr = Z
  IF( x<0._SP ) zr = -Z
  zb = zr
  yy = AIMAG(zr)
  ax = ABS(x)*1.7321_SP
  ay = ABS(yy)
  iform = 1
  IF( ay>ax ) iform = 2
  gnu = MAX(Fnu,1._SP)
  IF( Ikflg/=1 ) THEN
    fnn = nn
    gnn = Fnu + fnn - 1._SP
    gnu = MAX(gnn,fnn)
  END IF
  !-----------------------------------------------------------------------
  !     ONLY THE MAGNITUDE OF ARG AND PHI ARE NEEDED ALONG WITH THE
  !     REAL PARTS OF ZETA1, ZETA2 AND ZB. NO ATTEMPT IS MADE TO GET
  !     THE SIGN OF THE IMAGINARY PART CORRECT.
  !-----------------------------------------------------------------------
  IF( iform==2 ) THEN
    zn = -zr*CMPLX(0._SP,1._SP,SP)
    IF( yy<=0._SP ) zn = CONJG(-zn)
    CALL CUNHJ(zn,gnu,1,Tol,phi,arg,zeta1,zeta2,asum,bsum)
    cz = -zeta1 + zeta2
    aarg = ABS(arg)
  ELSE
    init = 0
    CALL CUNIK(zr,gnu,Ikflg,1,Tol,init,phi,zeta1,zeta2,summ,cwrk)
    cz = -zeta1 + zeta2
  END IF
  IF( Kode==2 ) cz = cz - zb
  IF( Ikflg==2 ) cz = -cz
  aphi = ABS(phi)
  rcz = REAL(cz)
  !-----------------------------------------------------------------------
  !     OVERFLOW TEST
  !-----------------------------------------------------------------------
  IF( rcz>Elim ) THEN
    Nuf = -1
    RETURN
  ELSE
    IF( rcz<Alim ) THEN
      !-----------------------------------------------------------------------
      !     UNDERFLOW TEST
      !-----------------------------------------------------------------------
      IF( rcz>=(-Elim) ) THEN
        IF( rcz>(-Alim) ) GOTO 50
        rcz = rcz + LOG(aphi)
        IF( iform==2 ) rcz = rcz - 0.25_SP*LOG(aarg) - aic
        IF( rcz>(-Elim) ) THEN
          ascle = 1.E+3_SP*tiny_sp/Tol
          cz = cz + LOG(phi)
          IF( iform/=1 ) cz = cz - CMPLX(0.25_SP,0._SP,SP)*LOG(arg)- CMPLX(aic,0._SP,SP)
          ax = EXP(rcz)/Tol
          ay = AIMAG(cz)
          cz = CMPLX(ax,0._SP,SP)*CMPLX(COS(ay),SIN(ay),SP)
          CALL CUCHK(cz,nw,ascle,Tol)
          IF( nw/=1 ) GOTO 50
        END IF
      END IF
      DO i = 1, nn
        Y(i) = (0._SP,0._SP)
      END DO
      Nuf = nn
      RETURN
    ELSE
      rcz = rcz + LOG(aphi)
      IF( iform==2 ) rcz = rcz - 0.25_SP*LOG(aarg) - aic
      IF( rcz>Elim ) THEN
        Nuf = -1
        RETURN
      END IF
    END IF
    50  IF( Ikflg==2 ) RETURN
    IF( N==1 ) RETURN
  END IF
  !-----------------------------------------------------------------------
  !     SET UNDERFLOWS ON I SEQUENCE
  !-----------------------------------------------------------------------
  100  gnu = Fnu + (nn-1)
  IF( iform==2 ) THEN
    CALL CUNHJ(zn,gnu,1,Tol,phi,arg,zeta1,zeta2,asum,bsum)
    cz = -zeta1 + zeta2
    aarg = ABS(arg)
  ELSE
    init = 0
    CALL CUNIK(zr,gnu,Ikflg,1,Tol,init,phi,zeta1,zeta2,summ,cwrk)
    cz = -zeta1 + zeta2
  END IF
  IF( Kode==2 ) cz = cz - zb
  aphi = ABS(phi)
  rcz = REAL(cz)
  IF( rcz>=(-Elim) ) THEN
    IF( rcz>(-Alim) ) RETURN
    rcz = rcz + LOG(aphi)
    IF( iform==2 ) rcz = rcz - 0.25_SP*LOG(aarg) - aic
    IF( rcz>(-Elim) ) THEN
      ascle = 1.E+3_SP*tiny_sp/Tol
      cz = cz + LOG(phi)
      IF( iform/=1 ) cz = cz - CMPLX(0.25_SP,0._SP,SP)*LOG(arg)- CMPLX(aic,0._SP,SP)
      ax = EXP(rcz)/Tol
      ay = AIMAG(cz)
      cz = CMPLX(ax,0._SP,SP)*CMPLX(COS(ay),SIN(ay),SP)
      CALL CUCHK(cz,nw,ascle,Tol)
      IF( nw/=1 ) RETURN
    END IF
  END IF
  Y(nn) = (0._SP,0._SP)
  nn = nn - 1
  Nuf = Nuf + 1
  IF( nn==0 ) RETURN
  GOTO 100

  RETURN
END SUBROUTINE CUOIK