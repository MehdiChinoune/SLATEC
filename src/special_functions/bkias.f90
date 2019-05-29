!** BKIAS
SUBROUTINE BKIAS(X,N,Ktrms,T,Ans,Ind,Ms,Gmrn,H,Ierr)
  !>
  !  Subsidiary to BSKIN
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (BKIAS-S, DBKIAS-D)
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !     BKIAS computes repeated integrals of the K0 Bessel function
  !     by the asymptotic expansion
  !
  !***
  ! **See also:**  BSKIN
  !***
  ! **Routines called:**  BDIFF, GAMRN, HKSEQ, R1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   820601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  USE service, ONLY : R1MACH
  INTEGER :: Ind, Ktrms, Ms, N, Ierr
  REAL :: Ans, Gmrn, H(30), T(50), X
  INTEGER i, ii, j, jmi, jn, k, kk, km, mm, mp
  REAL :: den1, den2, den3, er, err, fj, fk, fln, fm1, g1, gs, hn, rat, rg1, rxp, &
    rz, rzx, s(31), ss, sumi, sumj, tol, v(52), w(52), xp(16), z
  !-----------------------------------------------------------------------
  !             COEFFICIENTS OF POLYNOMIAL P(J-1,X), J=1,15
  !-----------------------------------------------------------------------
  REAL, PARAMETER :: b(120) = [ 1.00000000000000000E+00, 1.00000000000000000E+00, &
    -2.00000000000000000E+00, 1.00000000000000000E+00, -8.00000000000000000E+00, &
    6.00000000000000000E+00, 1.00000000000000000E+00, -2.20000000000000000E+01, &
    5.80000000000000000E+01, -2.40000000000000000E+01, 1.00000000000000000E+00, &
    -5.20000000000000000E+01, 3.28000000000000000E+02, -4.44000000000000000E+02, &
    1.20000000000000000E+02, 1.00000000000000000E+00, -1.14000000000000000E+02, &
    1.45200000000000000E+03, -4.40000000000000000E+03, 3.70800000000000000E+03, &
    -7.20000000000000000E+02, 1.00000000000000000E+00, -2.40000000000000000E+02, &
    5.61000000000000000E+03, -3.21200000000000000E+04, 5.81400000000000000E+04, &
    -3.39840000000000000E+04, 5.04000000000000000E+03, 1.00000000000000000E+00, &
    -4.94000000000000000E+02, 1.99500000000000000E+04, -1.95800000000000000E+05, &
    6.44020000000000000E+05, -7.85304000000000000E+05, 3.41136000000000000E+05, &
    -4.03200000000000000E+04, 1.00000000000000000E+00, -1.00400000000000000E+03, &
    6.72600000000000000E+04, -1.06250000000000000E+06, 5.76550000000000000E+06, &
    -1.24400640000000000E+07, 1.10262960000000000E+07, -3.73392000000000000E+06, &
    3.62880000000000000E+05, 1.00000000000000000E+00, -2.02600000000000000E+03, &
    2.18848000000000000E+05, -5.32616000000000000E+06, 4.47650000000000000E+07, &
    -1.55357384000000000E+08, 2.38904904000000000E+08, -1.62186912000000000E+08, &
    4.43390400000000000E+07, -3.62880000000000000E+06, 1.00000000000000000E+00, &
    -4.07200000000000000E+03, 6.95038000000000000E+05, -2.52439040000000000E+07, &
    3.14369720000000000E+08, -1.64838430400000000E+09, 4.00269508800000000E+09, &
    -4.64216395200000000E+09, 2.50748121600000000E+09, -5.68356480000000000E+08, &
    3.99168000000000000E+07, 1.00000000000000000E+00, -8.16600000000000000E+03, &
    2.17062600000000000E+06, -1.14876376000000000E+08, 2.05148277600000000E+09, &
    -1.55489607840000000E+10, 5.60413987840000000E+10, -1.01180433024000000E+11, &
    9.21997902240000000E+10, -4.07883018240000000E+10, 7.82771904000000000E+09, &
    -4.79001600000000000E+08, 1.00000000000000000E+00, -1.63560000000000000E+04, &
    6.69969600000000000E+06, -5.07259276000000000E+08, 1.26698177760000000E+10, &
    -1.34323420224000000E+11, 6.87720046384000000E+11, -1.81818864230400000E+12, &
    2.54986547342400000E+12, -1.88307966182400000E+12, 6.97929436800000000E+11, &
    -1.15336085760000000E+11, 6.22702080000000000E+09, 1.00000000000000000E+00, &
    -3.27380000000000000E+04, 2.05079880000000000E+07, -2.18982980800000000E+09, &
    7.50160522280000000E+10, -1.08467651241600000E+12, 7.63483214939200000E+12, &
    -2.82999100661120000E+13, 5.74943734645920000E+13, -6.47283751398720000E+13, &
    3.96895780558080000E+13, -1.25509040179200000E+13, 1.81099255680000000E+12, &
    -8.71782912000000000E+10, 1.00000000000000000E+00, -6.55040000000000000E+04, &
    6.24078900000000000E+07, -9.29252692000000000E+09, 4.29826006340000000E+11, &
    -8.30844432796800000E+12, 7.83913848313120000E+13, -3.94365587815520000E+14, &
    1.11174747256968000E+15, -1.79717122069056000E+15, 1.66642448627145600E+15, &
    -8.65023253219584000E+14, 2.36908271543040000E+14, -3.01963769856000000E+13, &
    1.30767436800000000E+12 ]
  !-----------------------------------------------------------------------
  !             BOUNDS B(M,K), K=M-3
  !-----------------------------------------------------------------------
  REAL, PARAMETER :: bnd(15) = [ 1.0E0, 1.0E0, 1.0E0, 1.0E0, 3.10E0, 5.18E0, 11.7E0, &
    29.8E0, 90.4E0, 297.0E0, 1070.0E0, 4290.0E0, 18100.0E0, 84700.0E0, 408000.0E0 ]
  REAL, PARAMETER :: hrtpi = 8.86226925452758014E-01
  !
  !* FIRST EXECUTABLE STATEMENT  BKIAS
  Ierr = 0
  tol = MAX(R1MACH(4),1.0E-18)
  fln = N
  rz = 1.0E0/(X+fln)
  rzx = X*rz
  z = 0.5E0*(X+fln)
  IF ( Ind<=1 ) Gmrn = GAMRN(z)
  gs = hrtpi*Gmrn
  g1 = gs + gs
  rg1 = 1.0E0/g1
  Gmrn = (rz+rz)/Gmrn
  IF ( Ind>1 ) GOTO 200
  !-----------------------------------------------------------------------
  !     EVALUATE ERROR FOR M=MS
  !-----------------------------------------------------------------------
  hn = 0.5E0*fln
  den2 = Ktrms + Ktrms + N
  den3 = den2 - 2.0E0
  den1 = X + den2
  err = rg1*(X+X)/(den1-1.0E0)
  IF ( N/=0 ) rat = 1.0E0/(fln*fln)
  IF ( Ktrms/=0 ) THEN
    fj = Ktrms
    rat = 0.25E0/(hrtpi*den3*SQRT(fj))
  END IF
  err = err*rat
  fj = -3.0E0
  DO j = 1, 15
    IF ( j<=5 ) err = err/den1
    fm1 = MAX(1.0E0,fj)
    fj = fj + 1.0E0
    er = bnd(j)*err
    IF ( Ktrms==0 ) THEN
      er = er*(1.0E0+hn/fm1)
      IF ( er<tol ) GOTO 100
      IF ( j>=5 ) err = err/fln
    ELSE
      er = er/fm1
      IF ( er<tol ) GOTO 100
      IF ( j>=5 ) err = err/den3
    END IF
  END DO
  Ierr = 2
  RETURN
  100  Ms = j
  200  mm = Ms + Ms
  mp = mm + 1
  !-----------------------------------------------------------------------
  !     H(K)=(-Z)**(K)*(PSI(K-1,Z)-PSI(K-1,Z+0.5))/GAMMA(K), K=1,2,...,MM
  !-----------------------------------------------------------------------
  IF ( Ind>1 ) THEN
    rat = z/(z-0.5E0)
    rxp = rat
    DO i = 1, mm
      H(i) = rxp*(1.0E0-H(i))
      rxp = rxp*rat
    END DO
  ELSE
    CALL HKSEQ(z,mm,H,Ierr)
  END IF
  !-----------------------------------------------------------------------
  !     SCALED S SEQUENCE
  !-----------------------------------------------------------------------
  s(1) = 1.0E0
  fk = 1.0E0
  DO k = 2, mp
    ss = 0.0E0
    km = k - 1
    i = km
    DO j = 1, km
      ss = ss + s(j)*H(i)
      i = i - 1
    END DO
    s(k) = ss/fk
    fk = fk + 1.0E0
  END DO
  !-----------------------------------------------------------------------
  !     SCALED S-TILDA SEQUENCE
  !-----------------------------------------------------------------------
  IF ( Ktrms/=0 ) THEN
    fk = 0.0E0
    ss = 0.0E0
    rg1 = rg1/z
    DO k = 1, Ktrms
      v(k) = z/(z+fk)
      w(k) = T(k)*v(k)
      ss = ss + w(k)
      fk = fk + 1.0E0
    END DO
    s(1) = s(1) - ss*rg1
    DO i = 2, mp
      ss = 0.0E0
      DO k = 1, Ktrms
        w(k) = w(k)*v(k)
        ss = ss + w(k)
      END DO
      s(i) = s(i) - ss*rg1
    END DO
  END IF
  !-----------------------------------------------------------------------
  !     SUM ON J
  !-----------------------------------------------------------------------
  sumj = 0.0E0
  jn = 1
  rxp = 1.0E0
  xp(1) = 1.0E0
  DO j = 1, Ms
    jn = jn + j - 1
    xp(j+1) = xp(j)*rzx
    rxp = rxp*rz
    !-----------------------------------------------------------------------
    !     SUM ON I
    !-----------------------------------------------------------------------
    sumi = 0.0E0
    ii = jn
    DO i = 1, j
      jmi = j - i + 1
      kk = j + i + 1
      DO k = 1, jmi
        v(k) = s(kk)*xp(k)
        kk = kk + 1
      END DO
      CALL BDIFF(jmi,v)
      sumi = sumi + b(ii)*v(jmi)*xp(i+1)
      ii = ii + 1
    END DO
    sumj = sumj + sumi*rxp
  END DO
  Ans = gs*(s(1)-sumj)
  RETURN
END SUBROUTINE BKIAS
