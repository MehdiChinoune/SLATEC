!DECK PPADD
SUBROUTINE PPADD(N,Ierror,A,C,Cbp,Bp,Bh)
  IMPLICIT NONE
  REAL A, Bh, Bp, BSRH, C, CNV, db, EPS, PPSGF, PPSPF, psg, &
    PSGF, scnv, sgn, xl, xm, xr
  INTEGER i3, icv, Ierror, if, ig, IK, is, it, iz, j, K, modiz, &
    N, NCMplx, nhalf, NM, NPP, nt
  !***BEGIN PROLOGUE  PPADD
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to BLKTRI
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (PPADD-S)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !   PPADD computes the eigenvalues of the periodic tridiagonal matrix
  !   with coefficients AN,BN,CN.
  !
  !   N    is the order of the BH and BP polynomials.
  !   BP   contains the eigenvalues on output.
  !   CBP  is the same as BP except type complex.
  !   BH   is used to temporarily store the roots of the B HAT polynomial
  !        which enters through BP.
  !
  !***SEE ALSO  BLKTRI
  !***ROUTINES CALLED  BSRH, PPSGF, PPSPF, PSGF
  !***COMMON BLOCKS    CBLKT
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  PPADD
  !
  COMPLEX cx, fsg, hsg, dd, f, fp, fpp, cdis, r1, r2, r3, Cbp
  DIMENSION A(*), C(*), Bp(*), Bh(*), Cbp(*)
  COMMON /CBLKT / NPP, K, EPS, CNV, NM, NCMplx, IK
  EXTERNAL PSGF, PPSPF, PPSGF
  !***FIRST EXECUTABLE STATEMENT  PPADD
  scnv = SQRT(CNV)
  iz = N
  IF ( Bp(N)<Bp(1) ) THEN
    DO j = 1, N
      nt = N - j
      Bh(j) = Bp(nt+1)
    ENDDO
  ELSEIF ( Bp(N)==Bp(1) ) THEN
    GOTO 300
  ELSE
    DO j = 1, N
      Bh(j) = Bp(j)
    ENDDO
  ENDIF
  NCMplx = 0
  modiz = MOD(iz,2)
  is = 1
  IF ( modiz/=0 ) THEN
    IF ( A(1)<0 ) GOTO 100
    IF ( A(1)==0 ) GOTO 300
  ENDIF
  xl = Bh(1)
  db = Bh(3) - Bh(1)
  DO
    xl = xl - db
    IF ( PSGF(xl,iz,C,A,Bh)>0 ) THEN
      sgn = -1.
      Cbp(1) = CMPLX(BSRH(xl,Bh(1),iz,C,A,Bh,PSGF,sgn),0.)
      is = 2
      EXIT
    ENDIF
  ENDDO
  100 CONTINUE
  IF = iz - 1
  IF ( modiz/=0 ) THEN
    IF ( A(1)<0 ) THEN
    ELSEIF ( A(1)==0 ) THEN
      GOTO 300
    ELSE
      GOTO 200
    ENDIF
  ENDIF
  xr = Bh(iz)
  db = Bh(iz) - Bh(iz-2)
  DO
    xr = xr + db
    IF ( PSGF(xr,iz,C,A,Bh)>=0 ) THEN
      sgn = 1.
      Cbp(iz) = CMPLX(BSRH(Bh(iz),xr,iz,C,A,Bh,PSGF,sgn),0.)
      if = iz - 2
      EXIT
    ENDIF
  ENDDO
  200 CONTINUE
  DO ig = is, if, 2
    xl = Bh(ig)
    xr = Bh(ig+1)
    sgn = -1.
    xm = BSRH(xl,xr,iz,C,A,Bh,PPSPF,sgn)
    psg = PSGF(xm,iz,C,A,Bh)
    IF ( ABS(psg)>EPS ) THEN
      IF ( psg*PPSGF(xm,iz,C,A,Bh)<0 ) THEN
        !
        !     CASE OF A REAL ZERO
        !
        sgn = 1.
        Cbp(ig) = CMPLX(BSRH(Bh(ig),xm,iz,C,A,Bh,PSGF,sgn),0.)
        sgn = -1.
        Cbp(ig+1) = CMPLX(BSRH(xm,Bh(ig+1),iz,C,A,Bh,PSGF,sgn),0.)
        CYCLE
      ELSEIF ( psg*PPSGF(xm,iz,C,A,Bh)/=0 ) THEN
        !
        !     CASE OF A COMPLEX ZERO
        !
        it = 0
        icv = 0
        cx = CMPLX(xm,0.)
        GOTO 250
      ENDIF
    ENDIF
    !
    !     CASE OF A MULTIPLE ZERO
    !
    Cbp(ig) = CMPLX(xm,0.)
    Cbp(ig+1) = CMPLX(xm,0.)
    CYCLE
    250    fsg = (1.,0.)
    hsg = (1.,0.)
    fp = (0.,0.)
    fpp = (0.,0.)
    DO j = 1, iz
      dd = 1./(cx-Bh(j))
      fsg = fsg*A(j)*dd
      hsg = hsg*C(j)*dd
      fp = fp + dd
      fpp = fpp - dd*dd
    ENDDO
    IF ( modiz/=0 ) THEN
      f = (1.,0.) + fsg + hsg
    ELSE
      f = (1.,0.) - fsg - hsg
    ENDIF
    i3 = 0
    IF ( ABS(fp)>0 ) THEN
      i3 = 1
      r3 = -f/fp
    ENDIF
    IF ( ABS(fpp)<=0 ) THEN
      r1 = r3
    ELSE
      cdis = SQRT(fp**2-2.*f*fpp)
      r1 = cdis - fp
      r2 = -fp - cdis
      IF ( ABS(r1)<=ABS(r2) ) THEN
        r1 = r2/fpp
      ELSE
        r1 = r1/fpp
      ENDIF
      r2 = 2.*f/fpp/r1
      IF ( ABS(r2)<ABS(r1) ) r1 = r2
      IF ( i3>0 ) THEN
        IF ( ABS(r3)<ABS(r1) ) r1 = r3
      ENDIF
    ENDIF
    cx = cx + r1
    it = it + 1
    IF ( it>50 ) GOTO 300
    IF ( ABS(r1)>scnv ) GOTO 250
    IF ( icv<=0 ) THEN
      icv = 1
      GOTO 250
    ELSE
      Cbp(ig) = cx
      Cbp(ig+1) = CONJG(cx)
    ENDIF
  ENDDO
  IF ( ABS(Cbp(N))<ABS(Cbp(1)) ) THEN
    nhalf = N/2
    DO j = 1, nhalf
      nt = N - j
      cx = Cbp(j)
      Cbp(j) = Cbp(nt+1)
      Cbp(nt+1) = cx
    ENDDO
  ELSEIF ( ABS(Cbp(N))==ABS(Cbp(1)) ) THEN
    GOTO 300
  ENDIF
  NCMplx = 1
  DO j = 2, iz
    IF ( AIMAG(Cbp(j))/=0 ) GOTO 99999
  ENDDO
  NCMplx = 0
  DO j = 2, iz
    Bp(j) = REAL(Cbp(j))
  ENDDO
  GOTO 99999
  300  Ierror = 4
  99999 CONTINUE
  END SUBROUTINE PPADD
