!** ZACON
SUBROUTINE ZACON(Zr,Zi,Fnu,Kode,Mr,N,Yr,Yi,Nz,Rl,Fnul,Tol,Elim,Alim)
  !>
  !  Subsidiary to ZBESH and ZBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CACON-A, ZACON-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     ZACON APPLIES THE ANALYTIC CONTINUATION FORMULA
  !
  !         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
  !                 MP=PI*MR*CMPLX(0.0,1.0)
  !
  !     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
  !     HALF Z PLANE
  !
  !***
  ! **See also:**  ZBESH, ZBESK
  !***
  ! **Routines called:**  D1MACH, ZABS, ZBINU, ZBKNU, ZMLT, ZS1S2

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : XERMSG, D1MACH
  !     COMPLEX CK,CONE,CSCL,CSCR,CSGN,CSPN,CY,CZERO,C1,C2,RZ,SC1,SC2,ST,
  !    *S1,S2,Y,Z,ZN
  INTEGER i, inu, iuf, kflag, Kode, Mr, N, nn, nw, Nz
  REAL(8) :: Alim, arg, ascle, as2, azn, bry(3), bscle, cki, ckr, cpn, cscl, &
    cscr, csgni, csgnr, cspni, cspnr, csr, csrr(3), cssr(3), cyi(2), cyr(2), &
    c1i, c1m, c1r, c2i, c2r, Elim, fmr, fn, Fnu, Fnul, pti, ptr, razn, Rl, rzi, &
    rzr, sc1i, sc1r, sc2i, sc2r, sgn, spn, sti, str, s1i, s1r, s2i, s2r, Tol, &
    Yi(N), Yr(N), yy, Zi, zni, znr, Zr
  REAL(8), PARAMETER :: pi = 3.14159265358979324D0
  REAL(8), PARAMETER :: zeror = 0.0D0, coner = 1.0D0
  !* FIRST EXECUTABLE STATEMENT  ZACON
  Nz = 0
  znr = -Zr
  zni = -Zi
  nn = N
  CALL ZBINU(znr,zni,Fnu,Kode,nn,Yr,Yi,nw,Rl,Fnul,Tol,Elim,Alim)
  IF ( nw>=0 ) THEN
    !-----------------------------------------------------------------------
    !     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
    !-----------------------------------------------------------------------
    nn = MIN(2,N)
    CALL ZBKNU(znr,zni,Fnu,Kode,nn,cyr,cyi,nw,Tol,Elim,Alim)
    IF ( nw==0 ) THEN
      s1r = cyr(1)
      s1i = cyi(1)
      fmr = Mr
      sgn = -SIGN(pi,fmr)
      csgnr = zeror
      csgni = sgn
      IF ( Kode/=1 ) THEN
        yy = -zni
        cpn = COS(yy)
        spn = SIN(yy)
        CALL ZMLT(csgnr,csgni,cpn,spn,csgnr,csgni)
      END IF
      !-----------------------------------------------------------------------
      !     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
      !     WHEN FNU IS LARGE
      !-----------------------------------------------------------------------
      inu = INT( Fnu )
      arg = (Fnu-inu)*sgn
      cpn = COS(arg)
      spn = SIN(arg)
      cspnr = cpn
      cspni = spn
      IF ( MOD(inu,2)/=0 ) THEN
        cspnr = -cspnr
        cspni = -cspni
      END IF
      iuf = 0
      c1r = s1r
      c1i = s1i
      c2r = Yr(1)
      c2i = Yi(1)
      ascle = 1.0D+3*D1MACH(1)/Tol
      IF ( Kode/=1 ) THEN
        CALL ZS1S2(znr,zni,c1r,c1i,c2r,c2i,nw,ascle,Alim,iuf)
        Nz = Nz + nw
        sc1r = c1r
        sc1i = c1i
      END IF
      CALL ZMLT(cspnr,cspni,c1r,c1i,str,sti)
      CALL ZMLT(csgnr,csgni,c2r,c2i,ptr,pti)
      Yr(1) = str + ptr
      Yi(1) = sti + pti
      IF ( N==1 ) RETURN
      cspnr = -cspnr
      cspni = -cspni
      s2r = cyr(2)
      s2i = cyi(2)
      c1r = s2r
      c1i = s2i
      c2r = Yr(2)
      c2i = Yi(2)
      IF ( Kode/=1 ) THEN
        CALL ZS1S2(znr,zni,c1r,c1i,c2r,c2i,nw,ascle,Alim,iuf)
        Nz = Nz + nw
        sc2r = c1r
        sc2i = c1i
      END IF
      CALL ZMLT(cspnr,cspni,c1r,c1i,str,sti)
      CALL ZMLT(csgnr,csgni,c2r,c2i,ptr,pti)
      Yr(2) = str + ptr
      Yi(2) = sti + pti
      IF ( N==2 ) RETURN
      cspnr = -cspnr
      cspni = -cspni
      azn = ZABS(znr,zni)
      razn = 1.0D0/azn
      str = znr*razn
      sti = -zni*razn
      rzr = (str+str)*razn
      rzi = (sti+sti)*razn
      fn = Fnu + 1.0D0
      ckr = fn*rzr
      cki = fn*rzi
      !-----------------------------------------------------------------------
      !     SCALE NEAR EXPONENT EXTREMES DURING RECURRENCE ON K FUNCTIONS
      !-----------------------------------------------------------------------
      cscl = 1.0D0/Tol
      cscr = Tol
      cssr(1) = cscl
      cssr(2) = coner
      cssr(3) = cscr
      csrr(1) = cscr
      csrr(2) = coner
      csrr(3) = cscl
      bry(1) = ascle
      bry(2) = 1.0D0/ascle
      bry(3) = D1MACH(2)
      as2 = ZABS(s2r,s2i)
      kflag = 2
      IF ( as2<=bry(1) ) THEN
        kflag = 1
      ELSEIF ( as2>=bry(2) ) THEN
        kflag = 3
      END IF
      bscle = bry(kflag)
      s1r = s1r*cssr(kflag)
      s1i = s1i*cssr(kflag)
      s2r = s2r*cssr(kflag)
      s2i = s2i*cssr(kflag)
      csr = csrr(kflag)
      DO i = 3, N
        str = s2r
        sti = s2i
        s2r = ckr*str - cki*sti + s1r
        s2i = ckr*sti + cki*str + s1i
        s1r = str
        s1i = sti
        c1r = s2r*csr
        c1i = s2i*csr
        str = c1r
        sti = c1i
        c2r = Yr(i)
        c2i = Yi(i)
        IF ( Kode/=1 ) THEN
          IF ( iuf>=0 ) THEN
            CALL ZS1S2(znr,zni,c1r,c1i,c2r,c2i,nw,ascle,Alim,iuf)
            Nz = Nz + nw
            sc1r = sc2r
            sc1i = sc2i
            sc2r = c1r
            sc2i = c1i
            IF ( iuf==3 ) THEN
              iuf = -4
              s1r = sc1r*cssr(kflag)
              s1i = sc1i*cssr(kflag)
              s2r = sc2r*cssr(kflag)
              s2i = sc2i*cssr(kflag)
              str = sc2r
              sti = sc2i
            END IF
          END IF
        END IF
        ptr = cspnr*c1r - cspni*c1i
        pti = cspnr*c1i + cspni*c1r
        Yr(i) = ptr + csgnr*c2r - csgni*c2i
        Yi(i) = pti + csgnr*c2i + csgni*c2r
        ckr = ckr + rzr
        cki = cki + rzi
        cspnr = -cspnr
        cspni = -cspni
        IF ( kflag<3 ) THEN
          ptr = ABS(c1r)
          pti = ABS(c1i)
          c1m = MAX(ptr,pti)
          IF ( c1m>bscle ) THEN
            kflag = kflag + 1
            bscle = bry(kflag)
            s1r = s1r*csr
            s1i = s1i*csr
            s2r = str
            s2i = sti
            s1r = s1r*cssr(kflag)
            s1i = s1i*cssr(kflag)
            s2r = s2r*cssr(kflag)
            s2i = s2i*cssr(kflag)
            csr = csrr(kflag)
          END IF
        END IF
      END DO
      RETURN
    END IF
  END IF
  Nz = -1
  IF ( nw==(-2) ) Nz = -2
END SUBROUTINE ZACON
