!*==ZUOIK.f90  processed by SPAG 6.72Dc at 11:02 on  6 Feb 2019
!DECK ZUOIK
SUBROUTINE ZUOIK(Zr,Zi,Fnu,Kode,Ikflg,N,Yr,Yi,Nuf,Tol,Elim,Alim)
  IMPLICIT NONE
  !*--ZUOIK5
  !***BEGIN PROLOGUE  ZUOIK
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to ZBESH, ZBESI and ZBESK
  !***LIBRARY   SLATEC
  !***TYPE      ALL (CUOIK-A, ZUOIK-A)
  !***AUTHOR  Amos, D. E., (SNL)
  !***DESCRIPTION
  !
  !     ZUOIK COMPUTES THE LEADING TERMS OF THE UNIFORM ASYMPTOTIC
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
  !***SEE ALSO  ZBESH, ZBESI, ZBESK
  !***ROUTINES CALLED  D1MACH, ZABS, ZLOG, ZUCHK, ZUNHJ, ZUNIK
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !   930122  Added ZLOG to EXTERNAL statement.  (RWC)
  !***END PROLOGUE  ZUOIK
  !     COMPLEX ARG,ASUM,BSUM,CWRK,CZ,CZERO,PHI,SUM,Y,Z,ZB,ZETA1,ZETA2,ZN,
  !    *ZR
  DOUBLE PRECISION aarg , aic , Alim , aphi , argi , argr , asumi , asumr , &
    ascle , ax , ay , bsumi , bsumr , cwrki , cwrkr , czi , &
    czr , Elim , fnn , Fnu , gnn , gnu , phii , phir , rcz , &
    str , sti , sumi , sumr , Tol , Yi , Yr , zbi , zbr , &
    zeroi , zeror , zeta1i , zeta1r , zeta2i , zeta2r , Zi , &
    zni , znr , Zr , zri , zrr , D1MACH , ZABS
  INTEGER i , idum , iform , Ikflg , init , Kode , N , nn , Nuf , nw
  DIMENSION Yr(N) , Yi(N) , cwrkr(16) , cwrki(16)
  EXTERNAL ZABS , ZLOG
  DATA zeror , zeroi/0.0D0 , 0.0D0/
  DATA aic/1.265512123484645396D+00/
  !***FIRST EXECUTABLE STATEMENT  ZUOIK
  Nuf = 0
  nn = N
  zrr = Zr
  zri = Zi
  IF ( Zr<0.0D0 ) THEN
    zrr = -Zr
    zri = -Zi
  ENDIF
  zbr = zrr
  zbi = zri
  ax = ABS(Zr)*1.7321D0
  ay = ABS(Zi)
  iform = 1
  IF ( ay>ax ) iform = 2
  gnu = MAX(Fnu,1.0D0)
  IF ( Ikflg/=1 ) THEN
    fnn = nn
    gnn = Fnu + fnn - 1.0D0
    gnu = MAX(gnn,fnn)
  ENDIF
  !-----------------------------------------------------------------------
  !     ONLY THE MAGNITUDE OF ARG AND PHI ARE NEEDED ALONG WITH THE
  !     REAL PARTS OF ZETA1, ZETA2 AND ZB. NO ATTEMPT IS MADE TO GET
  !     THE SIGN OF THE IMAGINARY PART CORRECT.
  !-----------------------------------------------------------------------
  IF ( iform==2 ) THEN
    znr = zri
    zni = -zrr
    IF ( Zi<=0.0D0 ) znr = -znr
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
  ENDIF
  IF ( Kode/=1 ) THEN
    czr = czr - zbr
    czi = czi - zbi
  ENDIF
  IF ( Ikflg/=1 ) THEN
    czr = -czr
    czi = -czi
  ENDIF
  aphi = ZABS(phir,phii)
  rcz = czr
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
        rcz = rcz + LOG(aphi)
        IF ( iform==2 ) rcz = rcz - 0.25D0*LOG(aarg) - aic
        IF ( rcz>(-Elim) ) THEN
          ascle = 1.0D+3*D1MACH(1)/Tol
          CALL ZLOG(phir,phii,str,sti,idum)
          czr = czr + str
          czi = czi + sti
          IF ( iform/=1 ) THEN
            CALL ZLOG(argr,argi,str,sti,idum)
            czr = czr - 0.25D0*str - aic
            czi = czi - 0.25D0*sti
          ENDIF
          ax = EXP(rcz)/Tol
          ay = czi
          czr = ax*COS(ay)
          czi = ax*SIN(ay)
          CALL ZUCHK(czr,czi,nw,ascle,Tol)
          IF ( nw==0 ) GOTO 50
        ENDIF
      ENDIF
      DO i = 1 , nn
        Yr(i) = zeror
        Yi(i) = zeroi
      ENDDO
      Nuf = nn
      RETURN
    ELSE
      rcz = rcz + LOG(aphi)
      IF ( iform==2 ) rcz = rcz - 0.25D0*LOG(aarg) - aic
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
  ENDIF
  IF ( Kode/=1 ) THEN
    czr = czr - zbr
    czi = czi - zbi
  ENDIF
  aphi = ZABS(phir,phii)
  rcz = czr
  IF ( rcz>=(-Elim) ) THEN
    IF ( rcz>(-Alim) ) RETURN
    rcz = rcz + LOG(aphi)
    IF ( iform==2 ) rcz = rcz - 0.25D0*LOG(aarg) - aic
    IF ( rcz>(-Elim) ) THEN
      ascle = 1.0D+3*D1MACH(1)/Tol
      CALL ZLOG(phir,phii,str,sti,idum)
      czr = czr + str
      czi = czi + sti
      IF ( iform/=1 ) THEN
        CALL ZLOG(argr,argi,str,sti,idum)
        czr = czr - 0.25D0*str - aic
        czi = czi - 0.25D0*sti
      ENDIF
      ax = EXP(rcz)/Tol
      ay = czi
      czr = ax*COS(ay)
      czi = ax*SIN(ay)
      CALL ZUCHK(czr,czi,nw,ascle,Tol)
      IF ( nw==0 ) RETURN
    ENDIF
  ENDIF
  Yr(nn) = zeror
  Yi(nn) = zeroi
  nn = nn - 1
  Nuf = Nuf + 1
  IF ( nn==0 ) RETURN
  GOTO 100
  99999 END SUBROUTINE ZUOIK
