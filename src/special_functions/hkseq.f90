!** HKSEQ
SUBROUTINE HKSEQ(X,M,H,Ierr)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to BSKIN
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (HKSEQ-S, DHKSEQ-D)
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !   HKSEQ is an adaptation of subroutine PSIFN described in the
  !   reference below.  HKSEQ generates the sequence
  !   H(K,X) = (-X)**(K+1)*(PSI(K,X) PSI(K,X+0.5))/GAMMA(K+1), for
  !            K=0,...,M.
  !
  !***
  ! **See also:**  BSKIN
  !***
  ! **References:**  D. E. Amos, A portable Fortran subroutine for
  !                 derivatives of the Psi function, Algorithm 610, ACM
  !                 Transactions on Mathematical Software 9, 4 (1983),
  !                 pp. 494-502.
  !***
  ! **Routines called:**  I1MACH, R1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   820601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
  
  INTEGER i, Ierr, j, k, M, mx, nx
  INTEGER I1MACH
  REAL b, fk, fln, fn, fnp, H, hrx, rln, rxsq, r1m5, s, slope, &
    t, tk, trm, trmh, trmr, tst, u, v, wdtol, X, xdmy, xh, &
    xinc, xm, xmin, yint
  REAL R1MACH
  DIMENSION b(22), trm(22), trmr(25), trmh(25), u(25), v(25), H(*)
  SAVE b
  !-----------------------------------------------------------------------
  !             SCALED BERNOULLI NUMBERS 2.0*B(2K)*(1-2**(-2K))
  !-----------------------------------------------------------------------
  DATA b(1), b(2), b(3), b(4), b(5), b(6), b(7), b(8), b(9), &
    b(10), b(11), b(12), b(13), b(14), b(15), b(16), b(17), &
    b(18), b(19), b(20), b(21), b(22)/1.00000000000000000E+00, &
    -5.00000000000000000E-01, 2.50000000000000000E-01, &
    -6.25000000000000000E-02, 4.68750000000000000E-02, &
    -6.64062500000000000E-02, 1.51367187500000000E-01, &
    -5.06103515625000000E-01, 2.33319091796875000E+00, &
    -1.41840972900390625E+01, 1.09941936492919922E+02, &
    -1.05824747562408447E+03, 1.23842434241771698E+04, &
    -1.73160495905935764E+05, 2.85103429084961116E+06, &
    -5.45964619322445132E+07, 1.20316174668075304E+09, &
    -3.02326315271452307E+10, 8.59229286072319606E+11, &
    -2.74233104097776039E+13, 9.76664637943633248E+14, &
    -3.85931586838450360E+16/
  !
  !* FIRST EXECUTABLE STATEMENT  HKSEQ
  Ierr = 0
  wdtol = MAX(R1MACH(4),1.0E-18)
  fn = M - 1
  fnp = fn + 1.0E0
  !-----------------------------------------------------------------------
  !     COMPUTE XMIN
  !-----------------------------------------------------------------------
  r1m5 = R1MACH(5)
  rln = r1m5*I1MACH(11)
  rln = MIN(rln,18.06E0)
  fln = MAX(rln,3.0E0) - 3.0E0
  yint = 3.50E0 + 0.40E0*fln
  slope = 0.21E0 + fln*(0.0006038E0*fln+0.008677E0)
  xm = yint + slope*fn
  mx = INT(xm) + 1
  xmin = mx
  !-----------------------------------------------------------------------
  !     GENERATE H(M-1,XDMY)*XDMY**(M) BY THE ASYMPTOTIC EXPANSION
  !-----------------------------------------------------------------------
  xdmy = X
  xinc = 0.0E0
  IF ( X<xmin ) THEN
    nx = INT(X)
    xinc = xmin - nx
    xdmy = X + xinc
  ENDIF
  rxsq = 1.0E0/(xdmy*xdmy)
  hrx = 0.5E0/xdmy
  tst = 0.5E0*wdtol
  t = fnp*hrx
  !-----------------------------------------------------------------------
  !     INITIALIZE COEFFICIENT ARRAY
  !-----------------------------------------------------------------------
  s = t*b(3)
  IF ( ABS(s)>=tst ) THEN
    tk = 2.0E0
    DO k = 4, 22
      t = t*((tk+fn+1.0E0)/(tk+1.0E0))*((tk+fn)/(tk+2.0E0))*rxsq
      trm(k) = t*b(k)
      IF ( ABS(trm(k))<tst ) GOTO 100
      s = s + trm(k)
      tk = tk + 2.0E0
    ENDDO
    GOTO 200
  ENDIF
  100  H(M) = s + 0.5E0
  IF ( M/=1 ) THEN
    !-----------------------------------------------------------------------
    !     GENERATE LOWER DERIVATIVES, I.LT.M-1
    !-----------------------------------------------------------------------
    DO i = 2, M
      fnp = fn
      fn = fn - 1.0E0
      s = fnp*hrx*b(3)
      IF ( ABS(s)>=tst ) THEN
        fk = fnp + 3.0E0
        DO k = 4, 22
          trm(k) = trm(k)*fnp/fk
          IF ( ABS(trm(k))<tst ) GOTO 120
          s = s + trm(k)
          fk = fk + 2.0E0
        ENDDO
        GOTO 200
      ENDIF
      120      mx = M - i + 1
      H(mx) = s + 0.5E0
    ENDDO
  ENDIF
  IF ( xinc==0.0E0 ) RETURN
  !-----------------------------------------------------------------------
  !     RECUR BACKWARD FROM XDMY TO X
  !-----------------------------------------------------------------------
  xh = X + 0.5E0
  s = 0.0E0
  nx = INT(xinc)
  DO i = 1, nx
    trmr(i) = X/(X+nx-i)
    u(i) = trmr(i)
    trmh(i) = X/(xh+nx-i)
    v(i) = trmh(i)
    s = s + u(i) - v(i)
  ENDDO
  mx = nx + 1
  trmr(mx) = X/xdmy
  u(mx) = trmr(mx)
  H(1) = H(1)*trmr(mx) + s
  IF ( M==1 ) RETURN
  DO j = 2, M
    s = 0.0E0
    DO i = 1, nx
      trmr(i) = trmr(i)*u(i)
      trmh(i) = trmh(i)*v(i)
      s = s + trmr(i) - trmh(i)
    ENDDO
    trmr(mx) = trmr(mx)*u(mx)
    H(j) = H(j)*trmr(mx) + s
  ENDDO
  RETURN
  200  Ierr = 2
END SUBROUTINE HKSEQ
