!DECK ZMLRI
SUBROUTINE ZMLRI(Zr,Zi,Fnu,Kode,N,Yr,Yi,Nz,Tol)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  ZMLRI
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to ZBESI and ZBESK
  !***LIBRARY   SLATEC
  !***TYPE      ALL (CMLRI-A, ZMLRI-A)
  !***AUTHOR  Amos, D. E., (SNL)
  !***DESCRIPTION
  !
  !     ZMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY THE
  !     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES.
  !
  !***SEE ALSO  ZBESI, ZBESK
  !***ROUTINES CALLED  D1MACH, DGAMLN, ZABS, ZEXP, ZLOG, ZMLT
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !   930122  Added ZEXP and ZLOG to EXTERNAL statement.  (RWC)
  !***END PROLOGUE  ZMLRI
  !     COMPLEX CK,CNORM,CONE,CTWO,CZERO,PT,P1,P2,RZ,SUM,Y,Z
  REAL(8) :: ack, ak, ap, at, az, bk, cki, ckr, cnormi, &
    cnormr, conei, coner, fkap, fkk, flam, fnf, Fnu, &
    pti, ptr, p1i, p1r, p2i, p2r, raz, rho, rho2, &
    rzi, rzr, scle, sti, str, sumi, sumr, tfnf, Tol, &
    tst, Yi, Yr, zeroi, zeror, Zi, Zr, DGAMLN, &
    D1MACH, ZABS
  INTEGER i, iaz, idum, ifnu, inu, itime, k, kk, km, Kode, m, N, &
    Nz
  DIMENSION Yr(N), Yi(N)
  EXTERNAL ZABS, ZEXP, ZLOG
  DATA zeror, zeroi, coner, conei/0.0D0, 0.0D0, 1.0D0, 0.0D0/
  !***FIRST EXECUTABLE STATEMENT  ZMLRI
  scle = D1MACH(1)/Tol
  Nz = 0
  az = ZABS(Zr,Zi)
  iaz = INT( az )
  ifnu = INT( Fnu )
  inu = ifnu + N - 1
  at = iaz + 1.0D0
  raz = 1.0D0/az
  str = Zr*raz
  sti = -Zi*raz
  ckr = str*at*raz
  cki = sti*at*raz
  rzr = (str+str)*raz
  rzi = (sti+sti)*raz
  p1r = zeror
  p1i = zeroi
  p2r = coner
  p2i = conei
  ack = (at+1.0D0)*raz
  rho = ack + SQRT(ack*ack-1.0D0)
  rho2 = rho*rho
  tst = (rho2+rho2)/((rho2-1.0D0)*(rho-1.0D0))
  tst = tst/Tol
  !-----------------------------------------------------------------------
  !     COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES
  !-----------------------------------------------------------------------
  ak = at
  DO i = 1, 80
    ptr = p2r
    pti = p2i
    p2r = p1r - (ckr*ptr-cki*pti)
    p2i = p1i - (cki*ptr+ckr*pti)
    p1r = ptr
    p1i = pti
    ckr = ckr + rzr
    cki = cki + rzi
    ap = ZABS(p2r,p2i)
    IF ( ap>tst*ak*ak ) GOTO 100
    ak = ak + 1.0D0
  ENDDO
  Nz = -2
  RETURN
  100  i = i + 1
  k = 0
  IF ( inu>=iaz ) THEN
    !-----------------------------------------------------------------------
    !     COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS
    !-----------------------------------------------------------------------
    p1r = zeror
    p1i = zeroi
    p2r = coner
    p2i = conei
    at = inu + 1.0D0
    str = Zr*raz
    sti = -Zi*raz
    ckr = str*at*raz
    cki = sti*at*raz
    ack = at*raz
    tst = SQRT(ack/Tol)
    itime = 1
    DO k = 1, 80
      ptr = p2r
      pti = p2i
      p2r = p1r - (ckr*ptr-cki*pti)
      p2i = p1i - (ckr*pti+cki*ptr)
      p1r = ptr
      p1i = pti
      ckr = ckr + rzr
      cki = cki + rzi
      ap = ZABS(p2r,p2i)
      IF ( ap>=tst ) THEN
        IF ( itime==2 ) GOTO 200
        ack = ZABS(ckr,cki)
        flam = ack + SQRT(ack*ack-1.0D0)
        fkap = ap/ZABS(p1r,p1i)
        rho = MIN(flam,fkap)
        tst = tst*SQRT(rho/(rho*rho-1.0D0))
        itime = 2
      ENDIF
    ENDDO
    Nz = -2
    RETURN
  ENDIF
  !-----------------------------------------------------------------------
  !     BACKWARD RECURRENCE AND SUM NORMALIZING RELATION
  !-----------------------------------------------------------------------
  200  k = k + 1
  kk = MAX(i+iaz,k+inu)
  fkk = kk
  p1r = zeror
  p1i = zeroi
  !-----------------------------------------------------------------------
  !     SCALE P2 AND SUM BY SCLE
  !-----------------------------------------------------------------------
  p2r = scle
  p2i = zeroi
  fnf = Fnu - ifnu
  tfnf = fnf + fnf
  bk = DGAMLN(fkk+tfnf+1.0D0,idum) - DGAMLN(fkk+1.0D0,idum)&
    - DGAMLN(tfnf+1.0D0,idum)
  bk = EXP(bk)
  sumr = zeror
  sumi = zeroi
  km = kk - inu
  DO i = 1, km
    ptr = p2r
    pti = p2i
    p2r = p1r + (fkk+fnf)*(rzr*ptr-rzi*pti)
    p2i = p1i + (fkk+fnf)*(rzi*ptr+rzr*pti)
    p1r = ptr
    p1i = pti
    ak = 1.0D0 - tfnf/(fkk+tfnf)
    ack = bk*ak
    sumr = sumr + (ack+bk)*p1r
    sumi = sumi + (ack+bk)*p1i
    bk = ack
    fkk = fkk - 1.0D0
  ENDDO
  Yr(N) = p2r
  Yi(N) = p2i
  IF ( N/=1 ) THEN
    DO i = 2, N
      ptr = p2r
      pti = p2i
      p2r = p1r + (fkk+fnf)*(rzr*ptr-rzi*pti)
      p2i = p1i + (fkk+fnf)*(rzi*ptr+rzr*pti)
      p1r = ptr
      p1i = pti
      ak = 1.0D0 - tfnf/(fkk+tfnf)
      ack = bk*ak
      sumr = sumr + (ack+bk)*p1r
      sumi = sumi + (ack+bk)*p1i
      bk = ack
      fkk = fkk - 1.0D0
      m = N - i + 1
      Yr(m) = p2r
      Yi(m) = p2i
    ENDDO
  ENDIF
  IF ( ifnu>0 ) THEN
    DO i = 1, ifnu
      ptr = p2r
      pti = p2i
      p2r = p1r + (fkk+fnf)*(rzr*ptr-rzi*pti)
      p2i = p1i + (fkk+fnf)*(rzr*pti+rzi*ptr)
      p1r = ptr
      p1i = pti
      ak = 1.0D0 - tfnf/(fkk+tfnf)
      ack = bk*ak
      sumr = sumr + (ack+bk)*p1r
      sumi = sumi + (ack+bk)*p1i
      bk = ack
      fkk = fkk - 1.0D0
    ENDDO
  ENDIF
  ptr = Zr
  pti = Zi
  IF ( Kode==2 ) ptr = zeror
  CALL ZLOG(rzr,rzi,str,sti,idum)
  p1r = -fnf*str + ptr
  p1i = -fnf*sti + pti
  ap = DGAMLN(1.0D0+fnf,idum)
  ptr = p1r - ap
  pti = p1i
  !-----------------------------------------------------------------------
  !     THE DIVISION CEXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW
  !     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES
  !-----------------------------------------------------------------------
  p2r = p2r + sumr
  p2i = p2i + sumi
  ap = ZABS(p2r,p2i)
  p1r = 1.0D0/ap
  CALL ZEXP(ptr,pti,str,sti)
  ckr = str*p1r
  cki = sti*p1r
  ptr = p2r*p1r
  pti = -p2i*p1r
  CALL ZMLT(ckr,cki,ptr,pti,cnormr,cnormi)
  DO i = 1, N
    str = Yr(i)*cnormr - Yi(i)*cnormi
    Yi(i) = Yr(i)*cnormi + Yi(i)*cnormr
    Yr(i) = str
  ENDDO
  RETURN
END SUBROUTINE ZMLRI
