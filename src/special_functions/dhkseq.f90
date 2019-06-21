!** DHKSEQ
SUBROUTINE DHKSEQ(X,M,H,Ierr)
  !> Subsidiary to DBSKIN
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (HKSEQ-S, DHKSEQ-D)
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !   DHKSEQ is an adaptation of subroutine DPSIFN described in the
  !   reference below.  DHKSEQ generates the sequence
  !   H(K,X) = (-X)**(K+1)*(PSI(K,X) PSI(K,X+0.5))/GAMMA(K+1), for
  !            K=0,...,M.
  !
  !***
  ! **See also:**  DBSKIN
  !***
  ! **References:**  D. E. Amos, A portable Fortran subroutine for
  !                 derivatives of the Psi function, Algorithm 610, ACM
  !                 Transactions on Mathematical Software 9, 4 (1983),
  !                 pp. 494-502.
  !***
  ! **Routines called:**  D1MACH, I1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   820601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
  USE service, ONLY : D1MACH, I1MACH
  INTEGER :: Ierr, M
  REAL(DP) :: H(M), X
  INTEGER :: i, j, k, mx, nx
  REAL(DP) :: fk, fln, fn, fnp, hrx, rln, rxsq, r1m5, s, slope, t, tk, trm(22), &
    trmh(25), trmr(25), tst, u(25), v(25), wdtol, xdmy, xh, xinc, xm, xmin, yint
  !-----------------------------------------------------------------------
  !             SCALED BERNOULLI NUMBERS 2.0*B(2K)*(1-2**(-2K))
  !-----------------------------------------------------------------------
  REAL(DP), PARAMETER :: b(22) = [ 1.00000000000000000E+00_DP, -5.00000000000000000E-01_DP, &
    2.50000000000000000E-01_DP, -6.25000000000000000E-02_DP, 4.68750000000000000E-02_DP, &
    -6.64062500000000000E-02_DP, 1.51367187500000000E-01_DP, -5.06103515625000000E-01_DP, &
    2.33319091796875000E+00_DP, -1.41840972900390625E+01_DP, 1.09941936492919922E+02_DP, &
    -1.05824747562408447E+03_DP, 1.23842434241771698E+04_DP, -1.73160495905935764E+05_DP, &
    2.85103429084961116E+06_DP, -5.45964619322445132E+07_DP, 1.20316174668075304E+09_DP, &
    -3.02326315271452307E+10_DP, 8.59229286072319606E+11_DP, -2.74233104097776039E+13_DP, &
    9.76664637943633248E+14_DP, -3.85931586838450360E+16_DP ]
  !
  !* FIRST EXECUTABLE STATEMENT  DHKSEQ
  Ierr = 0
  wdtol = MAX(D1MACH(4),1.E-18_DP)
  fn = M - 1
  fnp = fn + 1._DP
  !-----------------------------------------------------------------------
  !     COMPUTE XMIN
  !-----------------------------------------------------------------------
  r1m5 = D1MACH(5)
  rln = r1m5*I1MACH(14)
  rln = MIN(rln,18.06_DP)
  fln = MAX(rln,3._DP) - 3._DP
  yint = 3.50_DP + 0.40_DP*fln
  slope = 0.21_DP + fln*(0.0006038_DP*fln+0.008677_DP)
  xm = yint + slope*fn
  mx = INT(xm) + 1
  xmin = mx
  !-----------------------------------------------------------------------
  !     GENERATE H(M-1,XDMY)*XDMY**(M) BY THE ASYMPTOTIC EXPANSION
  !-----------------------------------------------------------------------
  xdmy = X
  xinc = 0._DP
  IF( X<xmin ) THEN
    nx = INT(X)
    xinc = xmin - nx
    xdmy = X + xinc
  END IF
  rxsq = 1._DP/(xdmy*xdmy)
  hrx = 0.5_DP/xdmy
  tst = 0.5_DP*wdtol
  t = fnp*hrx
  !-----------------------------------------------------------------------
  !     INITIALIZE COEFFICIENT ARRAY
  !-----------------------------------------------------------------------
  s = t*b(3)
  IF( ABS(s)>=tst ) THEN
    tk = 2._DP
    DO k = 4, 22
      t = t*((tk+fn+1._DP)/(tk+1._DP))*((tk+fn)/(tk+2._DP))*rxsq
      trm(k) = t*b(k)
      IF( ABS(trm(k))<tst ) GOTO 100
      s = s + trm(k)
      tk = tk + 2._DP
    END DO
    GOTO 200
  END IF
  100  H(M) = s + 0.5_DP
  IF( M/=1 ) THEN
    !-----------------------------------------------------------------------
    !     GENERATE LOWER DERIVATIVES, I<M-1
    !-----------------------------------------------------------------------
    DO i = 2, M
      fnp = fn
      fn = fn - 1._DP
      s = fnp*hrx*b(3)
      IF( ABS(s)>=tst ) THEN
        fk = fnp + 3._DP
        DO k = 4, 22
          trm(k) = trm(k)*fnp/fk
          IF( ABS(trm(k))<tst ) GOTO 120
          s = s + trm(k)
          fk = fk + 2._DP
        END DO
        GOTO 200
      END IF
      120  mx = M - i + 1
      H(mx) = s + 0.5_DP
    END DO
  END IF
  IF( xinc==0._DP ) RETURN
  !-----------------------------------------------------------------------
  !     RECUR BACKWARD FROM XDMY TO X
  !-----------------------------------------------------------------------
  xh = X + 0.5_DP
  s = 0._DP
  nx = INT(xinc)
  DO i = 1, nx
    trmr(i) = X/(X+nx-i)
    u(i) = trmr(i)
    trmh(i) = X/(xh+nx-i)
    v(i) = trmh(i)
    s = s + u(i) - v(i)
  END DO
  mx = nx + 1
  trmr(mx) = X/xdmy
  u(mx) = trmr(mx)
  H(1) = H(1)*trmr(mx) + s
  IF( M==1 ) RETURN
  DO j = 2, M
    s = 0._DP
    DO i = 1, nx
      trmr(i) = trmr(i)*u(i)
      trmh(i) = trmh(i)*v(i)
      s = s + trmr(i) - trmh(i)
    END DO
    trmr(mx) = trmr(mx)*u(mx)
    H(j) = H(j)*trmr(mx) + s
  END DO
  RETURN
  200  Ierr = 2
END SUBROUTINE DHKSEQ
