!** GAMRN
REAL(SP) FUNCTION GAMRN(X)
  !>
  !  Subsidiary to BSKIN
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (GAMRN-S, DGAMRN-D)
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract
  !         GAMRN computes the GAMMA function ratio GAMMA(X)/GAMMA(X+0.5)
  !         for real X.gt.0. If X.ge.XMIN, an asymptotic expansion is
  !         evaluated. If X.lt.XMIN, an integer is added to X to form a
  !         new value of X.ge.XMIN and the asymptotic expansion is eval-
  !         uated for this new value of X. Successive application of the
  !         recurrence relation
  !
  !                      W(X)=W(X+1)*(1+0.5/X)
  !
  !         reduces the argument to its original value. XMIN and comp-
  !         utational tolerances are computed as a function of the number
  !         of digits carried in a word by calls to I1MACH and R1MACH.
  !         However, the computational accuracy is limited to the max-
  !         imum of unit roundoff (=R1MACH(4)) and 1.0E-18 since critical
  !         constants are given to only 18 digits.
  !
  !         Input
  !           X      - Argument, X.gt.0.0
  !
  !         OUTPUT
  !           GAMRN  - Ratio  GAMMA(X)/GAMMA(X+0.5)
  !
  !***
  ! **See also:**  BSKIN
  !***
  ! **References:**  Y. L. Luke, The Special Functions and Their
  !                 Approximations, Vol. 1, Math In Sci. And
  !                 Eng. Series 53, Academic Press, New York, 1969,
  !                 pp. 34-35.
  !***
  ! **Routines called:**  I1MACH, R1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   820601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !   920520  Added REFERENCES section.  (WRB)
  USE service, ONLY : R1MACH, I1MACH
  INTEGER i, i1m11, k, mx, nx
  REAL(SP) fln, rln, s, tol, trm, X, xdmy, xinc, xm, xmin, xp, xsq
  !
  REAL(SP), PARAMETER :: gr(12) = [ 1.00000000000000000E+00,-1.56250000000000000E-02, &
    2.56347656250000000E-03,-1.27983093261718750E-03,  1.34351104497909546E-03, &
    -2.43289663922041655E-03, 6.75423753364157164E-03,-2.66369606131178216E-02, &
    1.41527455519564332E-01,-9.74384543032201613E-01,  8.43686251229783675E+00, &
    -8.97258321640552515E+01 ]
  !
  !* FIRST EXECUTABLE STATEMENT  GAMRN
  nx = INT(X)
  tol = MAX(R1MACH(4),1.0E-18)
  i1m11 = I1MACH(11)
  rln = R1MACH(5)*i1m11
  fln = MIN(rln,20.0E0)
  fln = MAX(fln,3.0E0)
  fln = fln - 3.0E0
  xm = 2.0E0 + fln*(0.2366E0+0.01723E0*fln)
  mx = INT(xm) + 1
  xmin = mx
  xdmy = X - 0.25E0
  xinc = 0.0E0
  IF ( X<xmin ) THEN
    xinc = xmin - nx
    xdmy = xdmy + xinc
  END IF
  s = 1.0E0
  IF ( xdmy*tol<=1.0E0 ) THEN
    xsq = 1.0E0/(xdmy*xdmy)
    xp = xsq
    DO k = 2, 12
      trm = gr(k)*xp
      IF ( ABS(trm)<tol ) EXIT
      s = s + trm
      xp = xp*xsq
    END DO
  END IF
  s = s/SQRT(xdmy)
  IF ( xinc/=0.0E0 ) THEN
    nx = INT(xinc)
    xp = 0.0E0
    DO i = 1, nx
      s = s*(1.0E0+0.5E0/(X+xp))
      xp = xp + 1.0E0
    END DO
    GAMRN = s
    RETURN
  END IF
  GAMRN = s
  RETURN
END FUNCTION GAMRN
