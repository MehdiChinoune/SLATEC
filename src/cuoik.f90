!DECK CUOIK
SUBROUTINE CUOIK(Z,Fnu,Kode,Ikflg,N,Y,Nuf,Tol,Elim,Alim)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CUOIK
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to CBESH, CBESI and CBESK
  !***LIBRARY   SLATEC
  !***TYPE      ALL (CUOIK-A, ZUOIK-A)
  !***AUTHOR  Amos, D. E., (SNL)
  !***DESCRIPTION
  !
  !     CUOIK COMPUTES THE LEADING TERMS OF THE UNIFORM ASYMPTOTIC
  !     EXPANSIONS FOR THE I AND K FUNCTIONS AND COMPARES THEM
  !     (IN LOGARITHMIC FORM) TO ALIM AND ELIM FOR OVER AND UNDERFLOW
  !     WHERE ALIM.LT.ELIM. IF THE MAGNITUDE, BASED ON THE LEADING
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
  !     IKFLG=1 AND NUF.GT.0 MEANS THE LAST NUF Y VALUES WERE SET TO ZERO
  !             THE FIRST N-NUF VALUES MUST BE SET BY ANOTHER ROUTINE
  !     IKFLG=2 AND NUF.EQ.N MEANS ALL Y VALUES WERE SET TO ZERO
  !     IKFLG=2 AND 0.LT.NUF.LT.N NOT CONSIDERED. Y MUST BE SET BY
  !             ANOTHER ROUTINE
  !
  !***SEE ALSO  CBESH, CBESI, CBESK
  !***ROUTINES CALLED  CUCHK, CUNHJ, CUNIK, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  CUOIK
  COMPLEX arg, asum, bsum, cwrk, cz, czero, phi, sum, Y, Z, zb, &
    zeta1, zeta2, zn, zr
  REAL aarg, aic, Alim, aphi, ascle, ax, ay, Elim, fnn, Fnu, gnn, &
    gnu, rcz, Tol, x, yy, R1MACH
  INTEGER i, iform, Ikflg, init, Kode, N, nn, Nuf, nw
  DIMENSION Y(N), cwrk(16)
  DATA czero/(0.0E0,0.0E0)/
  DATA aic/1.265512123484645396E+00/
  !***FIRST EXECUTABLE STATEMENT  CUOIK
  Nuf = 0
  nn = N
  x = REAL(Z)
  zr = Z
  IF ( x<0.0E0 ) zr = -Z
  zb = zr
  yy = AIMAG(zr)
  ax = ABS(x)*1.7321E0
  ay = ABS(yy)
  iform = 1
  IF ( ay>ax ) iform = 2
  gnu = MAX(Fnu,1.0E0)
  IF ( Ikflg/=1 ) THEN
    fnn = nn
    gnn = Fnu + fnn - 1.0E0
    gnu = MAX(gnn,fnn)
  ENDIF
  !-----------------------------------------------------------------------
  !     ONLY THE MAGNITUDE OF ARG AND PHI ARE NEEDED ALONG WITH THE
  !     REAL PARTS OF ZETA1, ZETA2 AND ZB. NO ATTEMPT IS MADE TO GET
  !     THE SIGN OF THE IMAGINARY PART CORRECT.
  !-----------------------------------------------------------------------
  IF ( iform==2 ) THEN
    zn = -zr*CMPLX(0.0E0,1.0E0)
    IF ( yy<=0.0E0 ) zn = CONJG(-zn)
    CALL CUNHJ(zn,gnu,1,Tol,phi,arg,zeta1,zeta2,asum,bsum)
    cz = -zeta1 + zeta2
    aarg = ABS(arg)
  ELSE
    init = 0
    CALL CUNIK(zr,gnu,Ikflg,1,Tol,init,phi,zeta1,zeta2,sum,cwrk)
    cz = -zeta1 + zeta2
  ENDIF
  IF ( Kode==2 ) cz = cz - zb
  IF ( Ikflg==2 ) cz = -cz
  aphi = ABS(phi)
  rcz = REAL(cz)
  !-----------------------------------------------------------------------
  !     OVERFLOW TEST
  !-----------------------------------------------------------------------
  IF ( rcz>Elim ) THEN
    Nuf = -1
    GOTO 99999
  ELSE
    IF ( rcz<Alim ) THEN
      !-----------------------------------------------------------------------
      !     UNDERFLOW TEST
      !-----------------------------------------------------------------------
      IF ( rcz>=(-Elim) ) THEN
        IF ( rcz>(-Alim) ) GOTO 50
        rcz = rcz + ALOG(aphi)
        IF ( iform==2 ) rcz = rcz - 0.25E0*ALOG(aarg) - aic
        IF ( rcz>(-Elim) ) THEN
          ascle = 1.0E+3*R1MACH(1)/Tol
          cz = cz + CLOG(phi)
          IF ( iform/=1 ) cz = cz - CMPLX(0.25E0,0.0E0)*CLOG(arg)&
            - CMPLX(aic,0.0E0)
          ax = EXP(rcz)/Tol
          ay = AIMAG(cz)
          cz = CMPLX(ax,0.0E0)*CMPLX(COS(ay),SIN(ay))
          CALL CUCHK(cz,nw,ascle,Tol)
          IF ( nw/=1 ) GOTO 50
        ENDIF
      ENDIF
      DO i = 1, nn
        Y(i) = czero
      ENDDO
      Nuf = nn
      RETURN
    ELSE
      rcz = rcz + ALOG(aphi)
      IF ( iform==2 ) rcz = rcz - 0.25E0*ALOG(aarg) - aic
      IF ( rcz>Elim ) THEN
        Nuf = -1
        GOTO 99999
      ENDIF
    ENDIF
    50     IF ( Ikflg==2 ) RETURN
    IF ( N==1 ) RETURN
  ENDIF
  !-----------------------------------------------------------------------
  !     SET UNDERFLOWS ON I SEQUENCE
  !-----------------------------------------------------------------------
  100  gnu = Fnu + (nn-1)
  IF ( iform==2 ) THEN
    CALL CUNHJ(zn,gnu,1,Tol,phi,arg,zeta1,zeta2,asum,bsum)
    cz = -zeta1 + zeta2
    aarg = ABS(arg)
  ELSE
    init = 0
    CALL CUNIK(zr,gnu,Ikflg,1,Tol,init,phi,zeta1,zeta2,sum,cwrk)
    cz = -zeta1 + zeta2
  ENDIF
  IF ( Kode==2 ) cz = cz - zb
  aphi = ABS(phi)
  rcz = REAL(cz)
  IF ( rcz>=(-Elim) ) THEN
    IF ( rcz>(-Alim) ) RETURN
    rcz = rcz + ALOG(aphi)
    IF ( iform==2 ) rcz = rcz - 0.25E0*ALOG(aarg) - aic
    IF ( rcz>(-Elim) ) THEN
      ascle = 1.0E+3*R1MACH(1)/Tol
      cz = cz + CLOG(phi)
      IF ( iform/=1 ) cz = cz - CMPLX(0.25E0,0.0E0)*CLOG(arg)&
        - CMPLX(aic,0.0E0)
      ax = EXP(rcz)/Tol
      ay = AIMAG(cz)
      cz = CMPLX(ax,0.0E0)*CMPLX(COS(ay),SIN(ay))
      CALL CUCHK(cz,nw,ascle,Tol)
      IF ( nw/=1 ) RETURN
    ENDIF
  ENDIF
  Y(nn) = czero
  nn = nn - 1
  Nuf = Nuf + 1
  IF ( nn==0 ) RETURN
  GOTO 100
  99999 CONTINUE
  END SUBROUTINE CUOIK
