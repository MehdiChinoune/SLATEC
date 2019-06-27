!** DGAMRN
REAL(DP) ELEMENTAL FUNCTION DGAMRN(X)
  !> Subsidiary to DBSKIN
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (GAMRN-S, DGAMRN-D)
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract   * A Double Precision Routine *
  !         DGAMRN computes the GAMMA function ratio GAMMA(X)/GAMMA(X+0.5)
  !         for real X>0. If X>=XMIN, an asymptotic expansion is
  !         evaluated. If X<XMIN, an integer is added to X to form a
  !         new value of X>=XMIN and the asymptotic expansion is eval-
  !         uated for this new value of X. Successive application of the
  !         recurrence relation
  !
  !                      W(X)=W(X+1)*(1+0.5/X)
  !
  !         reduces the argument to its original value. XMIN and comp-
  !         utational tolerances are computed as a function of the number
  !         of digits carried in a word by calls to I1MACH and D1MACH.
  !         However, the computational accuracy is limited to the max-
  !         imum of unit roundoff (=D1MACH(4)) and 1.0D-18 since critical
  !         constants are given to only 18 digits.
  !
  !         Input      X is Double Precision
  !           X      - Argument, X>0.0D0
  !
  !         Output      DGAMRN is DOUBLE PRECISION
  !           DGAMRN  - Ratio  GAMMA(X)/GAMMA(X+0.5)
  !
  !***
  ! **See also:**  DBSKIN
  !***
  ! **References:**  Y. L. Luke, The Special Functions and Their
  !                 Approximations, Vol. 1, Math In Sci. And
  !                 Eng. Series 53, Academic Press, New York, 1969,
  !                 pp. 34-35.
  !***
  ! **Routines called:**  D1MACH, I1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   820601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !   920520  Added REFERENCES section.  (WRB)
  USE service, ONLY : D1MACH, I1MACH
  REAL(DP), INTENT(IN) :: X
  INTEGER :: i, i1m11, k, mx, nx
  REAL(DP) :: fln, rln, s, tol, trm, xdmy, xinc, xm, xmin, xp, xsq
  !
  REAL(DP), PARAMETER :: gr(12) = [ 1.00000000000000000E+00_DP, -1.56250000000000000E-02_DP, &
    2.56347656250000000E-03_DP, -1.27983093261718750E-03_DP, 1.34351104497909546E-03_DP, &
    -2.43289663922041655E-03_DP, 6.75423753364157164E-03_DP, -2.66369606131178216E-02_DP, &
    1.41527455519564332E-01_DP, -9.74384543032201613E-01_DP, 8.43686251229783675E+00_DP, &
    -8.97258321640552515E+01_DP ]
  !
  !* FIRST EXECUTABLE STATEMENT  DGAMRN
  nx = INT(X)
  tol = MAX(D1MACH(4),1.E-18_DP)
  i1m11 = I1MACH(14)
  rln = D1MACH(5)*i1m11
  fln = MIN(rln,20._DP)
  fln = MAX(fln,3._DP)
  fln = fln - 3._DP
  xm = 2._DP + fln*(0.2366_DP+0.01723_DP*fln)
  mx = INT(xm) + 1
  xmin = mx
  xdmy = X - 0.25_DP
  xinc = 0._DP
  IF( X<xmin ) THEN
    xinc = xmin - nx
    xdmy = xdmy + xinc
  END IF
  s = 1._DP
  IF( xdmy*tol<=1._DP ) THEN
    xsq = 1._DP/(xdmy*xdmy)
    xp = xsq
    DO k = 2, 12
      trm = gr(k)*xp
      IF( ABS(trm)<tol ) EXIT
      s = s + trm
      xp = xp*xsq
    END DO
  END IF
  s = s/SQRT(xdmy)
  IF( xinc/=0._DP ) THEN
    nx = INT(xinc)
    xp = 0._DP
    DO i = 1, nx
      s = s*(1._DP+0.5_DP/(X+xp))
      xp = xp + 1._DP
    END DO
    DGAMRN = s
    RETURN
  END IF
  DGAMRN = s

  RETURN
END FUNCTION DGAMRN