!** ZKSCL
SUBROUTINE ZKSCL(Zrr,Zri,Fnu,N,Yr,Yi,Nz,Rzr,Rzi,Ascle,Tol,Elim)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to ZBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CKSCL-A, ZKSCL-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     SET K FUNCTIONS TO ZERO ON UNDERFLOW, CONTINUE RECURRENCE
  !     ON SCALED FUNCTIONS UNTIL TWO MEMBERS COME ON SCALE, THEN
  !     RETURN WITH MIN(NZ+2,N) VALUES SCALED BY 1/TOL.
  !
  !***
  ! **See also:**  ZBESK
  !***
  ! **Routines called:**  ZABS, ZLOG, ZUCHK

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !   930122  Added ZLOG to EXTERNAL statement.  (RWC)

  !     COMPLEX CK,CS,CY,CZERO,RZ,S1,S2,Y,ZR,ZD,CELM
  INTEGER i, ic, idum, kk, N, nn, nw, Nz
  REAL(8) :: acs, as, Ascle, cki, ckr, csi, csr, cyi(2), cyr(2), &
    Elim, fn, Fnu, Rzi, Rzr, str, s1i, s1r, s2i, &
    s2r, Tol, Yi(N), Yr(N), zeroi, zeror, Zri, Zrr, zdr, zdi, celmr, elm, helim, alas
  REAL(8), EXTERNAL :: ZABS
  EXTERNAL :: ZLOG
  DATA zeror, zeroi/0.0D0, 0.0D0/
  !* FIRST EXECUTABLE STATEMENT  ZKSCL
  Nz = 0
  ic = 0
  nn = MIN(2,N)
  DO i = 1, nn
    s1r = Yr(i)
    s1i = Yi(i)
    cyr(i) = s1r
    cyi(i) = s1i
    as = ZABS(s1r,s1i)
    acs = -Zrr + LOG(as)
    Nz = Nz + 1
    Yr(i) = zeror
    Yi(i) = zeroi
    IF ( acs>=(-Elim) ) THEN
      CALL ZLOG(s1r,s1i,csr,csi,idum)
      csr = csr - Zrr
      csi = csi - Zri
      str = EXP(csr)/Tol
      csr = str*COS(csi)
      csi = str*SIN(csi)
      CALL ZUCHK(csr,csi,nw,Ascle,Tol)
      IF ( nw==0 ) THEN
        Yr(i) = csr
        Yi(i) = csi
        ic = i
        Nz = Nz - 1
      ENDIF
    ENDIF
  ENDDO
  IF ( N==1 ) RETURN
  IF ( ic<=1 ) THEN
    Yr(1) = zeror
    Yi(1) = zeroi
    Nz = 2
  ENDIF
  IF ( N==2 ) RETURN
  IF ( Nz==0 ) RETURN
  fn = Fnu + 1.0D0
  ckr = fn*Rzr
  cki = fn*Rzi
  s1r = cyr(1)
  s1i = cyi(1)
  s2r = cyr(2)
  s2i = cyi(2)
  helim = 0.5D0*Elim
  elm = EXP(-Elim)
  celmr = elm
  zdr = Zrr
  zdi = Zri
  !
  !     FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE IF
  !     S2 GETS LARGER THAN EXP(ELIM/2)
  !
  DO i = 3, N
    kk = i
    csr = s2r
    csi = s2i
    s2r = ckr*csr - cki*csi + s1r
    s2i = cki*csr + ckr*csi + s1i
    s1r = csr
    s1i = csi
    ckr = ckr + Rzr
    cki = cki + Rzi
    as = ZABS(s2r,s2i)
    alas = LOG(as)
    acs = -zdr + alas
    Nz = Nz + 1
    Yr(i) = zeror
    Yi(i) = zeroi
    IF ( acs>=(-Elim) ) THEN
      CALL ZLOG(s2r,s2i,csr,csi,idum)
      csr = csr - zdr
      csi = csi - zdi
      str = EXP(csr)/Tol
      csr = str*COS(csi)
      csi = str*SIN(csi)
      CALL ZUCHK(csr,csi,nw,Ascle,Tol)
      IF ( nw==0 ) THEN
        Yr(i) = csr
        Yi(i) = csi
        Nz = Nz - 1
        IF ( ic==kk-1 ) GOTO 100
        ic = kk
        CYCLE
      ENDIF
    ENDIF
    IF ( alas>=helim ) THEN
      zdr = zdr - Elim
      s1r = s1r*celmr
      s1i = s1i*celmr
      s2r = s2r*celmr
      s2i = s2i*celmr
    ENDIF
  ENDDO
  Nz = N
  IF ( ic==N ) Nz = N - 1
  GOTO 200
  100  Nz = kk - 2
  200 CONTINUE
  DO i = 1, Nz
    Yr(i) = zeror
    Yi(i) = zeroi
  ENDDO
END SUBROUTINE ZKSCL
