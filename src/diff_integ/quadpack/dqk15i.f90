!** DQK15I
PURE SUBROUTINE DQK15I(F,Boun,Inf,A,B,Result,Abserr,Resabs,Resasc)
  !> The original (infinite integration range is mapped onto the
  !  interval (0,1) and (A,B) is a part of (0,1). it is the purpose to compute
  !  I = Integral of transformed integrand over (A,B),
  !  J = Integral of ABS(Transformed Integrand) over (A,B).
  !***
  ! **Library:**   SLATEC (QUADPACK)
  !***
  ! **Category:**  H2A3A2, H2A4A2
  !***
  ! **Type:**      DOUBLE PRECISION (QK15I-S, DQK15I-D)
  !***
  ! **Keywords:**  15-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
  !***
  ! **Author:**  Piessens, Robert
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !           de Doncker, Elise
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !***
  ! **Description:**
  !
  !           Integration Rule
  !           Standard Fortran subroutine
  !           Double precision version
  !
  !           PARAMETERS
  !            ON ENTRY
  !              F      - Double precision
  !                       Function subprogram defining the integrand
  !                       FUNCTION F(X). The actual name for F needs to be
  !                       Declared E X T E R N A L in the calling program.
  !
  !              BOUN   - Double precision
  !                       Finite bound of original integration
  !                       Range (SET TO ZERO IF INF = +2)
  !
  !              INF    - Integer
  !                       If INF = -1, the original interval is
  !                                   (-INFINITY,BOUND),
  !                       If INF = +1, the original interval is
  !                                   (BOUND,+INFINITY),
  !                       If INF = +2, the original interval is
  !                                   (-INFINITY,+INFINITY) AND
  !                       The integral is computed as the sum of two
  !                       integrals, one over (-INFINITY,0) and one over
  !                       (0,+INFINITY).
  !
  !              A      - Double precision
  !                       Lower limit for integration over subrange
  !                       of (0,1)
  !
  !              B      - Double precision
  !                       Upper limit for integration over subrange
  !                       of (0,1)
  !
  !            ON RETURN
  !              RESULT - Double precision
  !                       Approximation to the integral I
  !                       Result is computed by applying the 15-POINT
  !                       KRONROD RULE(RESK) obtained by optimal addition
  !                       of abscissae to the 7-POINT GAUSS RULE(RESG).
  !
  !              ABSERR - Double precision
  !                       Estimate of the modulus of the absolute error,
  !                       WHICH SHOULD EQUAL or EXCEED ABS(I-RESULT)
  !
  !              RESABS - Double precision
  !                       Approximation to the integral J
  !
  !              RESASC - Double precision
  !                       Approximation to the integral of
  !                       ABS((TRANSFORMED INTEGRAND)-I/(B-A)) over (A,B)
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   800101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : D1MACH
  !
  INTERFACE
    REAL(DP) PURE FUNCTION F(X)
      IMPORT DP
      REAL(DP), INTENT(IN) :: X
    END FUNCTION F
  END INTERFACE
  INTEGER, INTENT(IN) :: Inf
  REAL(DP), INTENT(IN) :: A, B, Boun
  REAL(DP), INTENT(OUT) :: Abserr, Resabs, Resasc, Result
  !
  INTEGER :: j
  REAL(DP) :: absc, absc1, absc2, centr, dinf, epmach, fc, fsum, fval1, fval2, &
    fv1(7), fv2(7), hlgth, resg, resk, reskh, tabsc1, tabsc2, uflow
  !
  !           THE ABSCISSAE AND WEIGHTS ARE SUPPLIED FOR THE INTERVAL
  !           (-1,1).  BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND
  !           THEIR CORRESPONDING WEIGHTS ARE GIVEN.
  !
  !           XGK    - ABSCISSAE OF THE 15-POINT KRONROD RULE
  !                    XGK(2), XGK(4), ... ABSCISSAE OF THE 7-POINT
  !                    GAUSS RULE
  !                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
  !                    ADDED TO THE 7-POINT GAUSS RULE
  !
  !           WGK    - WEIGHTS OF THE 15-POINT KRONROD RULE
  !
  !           WG     - WEIGHTS OF THE 7-POINT GAUSS RULE, CORRESPONDING
  !                    TO THE ABSCISSAE XGK(2), XGK(4), ...
  !                    WG(1), WG(3), ... ARE SET TO ZERO.
  !
  REAL(DP), PARAMETER :: xgk(8) = [ 0.9914553711208126E+00_DP, 0.9491079123427585E+00_DP, &
    0.8648644233597691E+00_DP, 0.7415311855993944E+00_DP, 0.5860872354676911E+00_DP, &
    0.4058451513773972E+00_DP, 0.2077849550078985E+00_DP, 0.0000000000000000E+00_DP ]
  REAL(DP), PARAMETER :: wgk(8) = [ 0.2293532201052922E-01_DP, 0.6309209262997855E-01_DP, &
    0.1047900103222502E+00_DP, 0.1406532597155259E+00_DP, 0.1690047266392679E+00_DP, &
    0.1903505780647854E+00_DP, 0.2044329400752989E+00_DP, 0.2094821410847278E+00_DP ]
  REAL(DP), PARAMETER :: wg(8) = [ 0.0000000000000000E+00_DP, 0.1294849661688697E+00_DP, &
    0.0000000000000000E+00_DP, 0.2797053914892767E+00_DP, 0.0000000000000000E+00_DP, &
    0.3818300505051189E+00_DP, 0.0000000000000000E+00_DP, 0.4179591836734694E+00_DP ]
  !
  !           LIST OF MAJOR VARIABLES
  !           -----------------------
  !
  !           CENTR  - MID POINT OF THE INTERVAL
  !           HLGTH  - HALF-LENGTH OF THE INTERVAL
  !           ABSC*  - ABSCISSA
  !           TABSC* - TRANSFORMED ABSCISSA
  !           FVAL*  - FUNCTION VALUE
  !           RESG   - RESULT OF THE 7-POINT GAUSS FORMULA
  !           RESK   - RESULT OF THE 15-POINT KRONROD FORMULA
  !           RESKH  - APPROXIMATION TO THE MEAN VALUE OF THE TRANSFORMED
  !                    INTEGRAND OVER (A,B), I.E. TO I/(B-A)
  !
  !           MACHINE DEPENDENT CONSTANTS
  !           ---------------------------
  !
  !           EPMACH IS THE LARGEST RELATIVE SPACING.
  !           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
  !
  !* FIRST EXECUTABLE STATEMENT  DQK15I
  epmach = D1MACH(4)
  uflow = D1MACH(1)
  dinf = MIN(1,Inf)
  !
  centr = 0.5_DP*(A+B)
  hlgth = 0.5_DP*(B-A)
  tabsc1 = Boun + dinf*(1._DP-centr)/centr
  fval1 = F(tabsc1)
  IF( Inf==2 ) fval1 = fval1 + F(-tabsc1)
  fc = (fval1/centr)/centr
  !
  !           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO
  !           THE INTEGRAL, AND ESTIMATE THE ERROR.
  !
  resg = wg(8)*fc
  resk = wgk(8)*fc
  Resabs = ABS(resk)
  DO j = 1, 7
    absc = hlgth*xgk(j)
    absc1 = centr - absc
    absc2 = centr + absc
    tabsc1 = Boun + dinf*(1._DP-absc1)/absc1
    tabsc2 = Boun + dinf*(1._DP-absc2)/absc2
    fval1 = F(tabsc1)
    fval2 = F(tabsc2)
    IF( Inf==2 ) fval1 = fval1 + F(-tabsc1)
    IF( Inf==2 ) fval2 = fval2 + F(-tabsc2)
    fval1 = (fval1/absc1)/absc1
    fval2 = (fval2/absc2)/absc2
    fv1(j) = fval1
    fv2(j) = fval2
    fsum = fval1 + fval2
    resg = resg + wg(j)*fsum
    resk = resk + wgk(j)*fsum
    Resabs = Resabs + wgk(j)*(ABS(fval1)+ABS(fval2))
  END DO
  reskh = resk*0.5_DP
  Resasc = wgk(8)*ABS(fc-reskh)
  DO j = 1, 7
    Resasc = Resasc + wgk(j)*(ABS(fv1(j)-reskh)+ABS(fv2(j)-reskh))
  END DO
  Result = resk*hlgth
  Resasc = Resasc*hlgth
  Resabs = Resabs*hlgth
  Abserr = ABS((resk-resg)*hlgth)
  IF( Resasc/=0._DP .AND. Abserr/=0._DP )&
    Abserr = Resasc*MIN(1._DP,(0.2E+03_DP*Abserr/Resasc)**1.5_DP)
  IF( Resabs>uflow/(0.5E+02_DP*epmach) ) Abserr = MAX((epmach*0.5E+02_DP)*Resabs,Abserr)
  !
END SUBROUTINE DQK15I