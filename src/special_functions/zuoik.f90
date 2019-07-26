!** ZUOIK
SUBROUTINE ZUOIK(Zr,Zi,Fnu,Kode,Ikflg,N,Yr,Yi,Nuf,Tol,Elim,Alim)
  !> Subsidiary to ZBESH, ZBESI and ZBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CUOIK-A, ZUOIK-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     ZUOIK COMPUTES THE LEADING TERMS OF THE UNIFORM ASYMPTOTIC
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
  ! **See also:**  ZBESH, ZBESI, ZBESK
  !***
  ! **Routines called:**  D1MACH, ZABS, ZLOG, ZUCHK, ZUNHJ, ZUNIK

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !   930122  Added ZLOG to EXTERNAL statement.  (RWC)
  USE service, ONLY : tiny_dp
  !     COMPLEX ARG,ASUM,BSUM,CWRK,CZ,CZERO,PHI,SUM,Y,Z,ZB,ZETA1,ZETA2,ZN,
  !    *ZR
  INTEGER :: i, iform, Ikflg, init, Kode, N, nn, Nuf, nw
  REAL(DP) :: aarg, Alim, aphi, argi, argr, asumi, asumr, ascle, ax, ay, bsumi, &
    bsumr, cwrki(16), cwrkr(16), czi, czr, Elim, fnn, Fnu, gnn, gnu, phii, phir, &
    rcz, str, sti, sumi, sumr, Tol, Yi(N), Yr(N), zbi, zbr, zeta1i, zeta1r, &
    zeta2i, zeta2r, Zi, zni, znr, Zr, zri, zrr
  REAL(DP), PARAMETER :: zeror = 0._DP, zeroi = 0._DP
  REAL(DP), PARAMETER :: aic = 1.265512123484645396_DP
  !* FIRST EXECUTABLE STATEMENT  ZUOIK
  Nuf = 0
  nn = N
  zrr = Zr
  zri = Zi
  IF( Zr<0._DP ) THEN
    zrr = -Zr
    zri = -Zi
  END IF
  zbr = zrr
  zbi = zri
  ax = ABS(Zr)*1.7321_DP
  ay = ABS(Zi)
  iform = 1
  IF( ay>ax ) iform = 2
  gnu = MAX(Fnu,1._DP)
  IF( Ikflg/=1 ) THEN
    fnn = nn
    gnn = Fnu + fnn - 1._DP
    gnu = MAX(gnn,fnn)
  END IF
  !-----------------------------------------------------------------------
  !     ONLY THE MAGNITUDE OF ARG AND PHI ARE NEEDED ALONG WITH THE
  !     REAL PARTS OF ZETA1, ZETA2 AND ZB. NO ATTEMPT IS MADE TO GET
  !     THE SIGN OF THE IMAGINARY PART CORRECT.
  !-----------------------------------------------------------------------
  IF( iform==2 ) THEN
    znr = zri
    zni = -zrr
    IF( Zi<=0._DP ) znr = -znr
    CALL ZUNHJ(znr,zni,gnu,1,Tol,phir,phii,argr,argi,zeta1r,zeta1i,zeta2r,&
      zeta2i,asumr,asumi,bsumr,bsumi)
    czr = -zeta1r + zeta2r
    czi = -zeta1i + zeta2i
    aarg = ZABS(argr,argi)
  ELSE
    init = 0
    CALL ZUNIK(zrr,zri,gnu,Ikflg,1,Tol,init,phir,phii,zeta1r,zeta1i,zeta2r,&
      zeta2i,sumr,sumi,cwrkr,cwrki)
    czr = -zeta1r + zeta2r
    czi = -zeta1i + zeta2i
  END IF
  IF( Kode/=1 ) THEN
    czr = czr - zbr
    czi = czi - zbi
  END IF
  IF( Ikflg/=1 ) THEN
    czr = -czr
    czi = -czi
  END IF
  aphi = ZABS(phir,phii)
  rcz = czr
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
        IF( iform==2 ) rcz = rcz - 0.25_DP*LOG(aarg) - aic
        IF( rcz>(-Elim) ) THEN
          ascle = 1.E3_DP*tiny_dp/Tol
          CALL ZLOG(phir,phii,str,sti)
          czr = czr + str
          czi = czi + sti
          IF( iform/=1 ) THEN
            CALL ZLOG(argr,argi,str,sti)
            czr = czr - 0.25_DP*str - aic
            czi = czi - 0.25_DP*sti
          END IF
          ax = EXP(rcz)/Tol
          ay = czi
          czr = ax*COS(ay)
          czi = ax*SIN(ay)
          CALL ZUCHK(czr,czi,nw,ascle,Tol)
          IF( nw==0 ) GOTO 50
        END IF
      END IF
      DO i = 1, nn
        Yr(i) = zeror
        Yi(i) = zeroi
      END DO
      Nuf = nn
      RETURN
    ELSE
      rcz = rcz + LOG(aphi)
      IF( iform==2 ) rcz = rcz - 0.25_DP*LOG(aarg) - aic
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
    CALL ZUNHJ(znr,zni,gnu,1,Tol,phir,phii,argr,argi,zeta1r,zeta1i,zeta2r,&
      zeta2i,asumr,asumi,bsumr,bsumi)
    czr = -zeta1r + zeta2r
    czi = -zeta1i + zeta2i
    aarg = ZABS(argr,argi)
  ELSE
    init = 0
    CALL ZUNIK(zrr,zri,gnu,Ikflg,1,Tol,init,phir,phii,zeta1r,zeta1i,zeta2r,&
      zeta2i,sumr,sumi,cwrkr,cwrki)
    czr = -zeta1r + zeta2r
    czi = -zeta1i + zeta2i
  END IF
  IF( Kode/=1 ) THEN
    czr = czr - zbr
    czi = czi - zbi
  END IF
  aphi = ZABS(phir,phii)
  rcz = czr
  IF( rcz>=(-Elim) ) THEN
    IF( rcz>(-Alim) ) RETURN
    rcz = rcz + LOG(aphi)
    IF( iform==2 ) rcz = rcz - 0.25_DP*LOG(aarg) - aic
    IF( rcz>(-Elim) ) THEN
      ascle = 1.E3_DP*tiny_dp/Tol
      CALL ZLOG(phir,phii,str,sti)
      czr = czr + str
      czi = czi + sti
      IF( iform/=1 ) THEN
        CALL ZLOG(argr,argi,str,sti)
        czr = czr - 0.25_DP*str - aic
        czi = czi - 0.25_DP*sti
      END IF
      ax = EXP(rcz)/Tol
      ay = czi
      czr = ax*COS(ay)
      czi = ax*SIN(ay)
      CALL ZUCHK(czr,czi,nw,ascle,Tol)
      IF( nw==0 ) RETURN
    END IF
  END IF
  Yr(nn) = zeror
  Yi(nn) = zeroi
  nn = nn - 1
  Nuf = Nuf + 1
  IF( nn==0 ) RETURN
  GOTO 100
  RETURN
END SUBROUTINE ZUOIK
