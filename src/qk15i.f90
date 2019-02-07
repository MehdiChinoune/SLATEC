!*==QK15I.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK QK15I
SUBROUTINE QK15I(F,Boun,Inf,A,B,Result,Abserr,Resabs,Resasc)
  IMPLICIT NONE
  !*--QK15I5
  !***BEGIN PROLOGUE  QK15I
  !***PURPOSE  The original (infinite integration range is mapped
  !            onto the interval (0,1) and (A,B) is a part of (0,1).
  !            it is the purpose to compute
  !            I = Integral of transformed integrand over (A,B),
  !            J = Integral of ABS(Transformed Integrand) over (A,B).
  !***LIBRARY   SLATEC (QUADPACK)
  !***CATEGORY  H2A3A2, H2A4A2
  !***TYPE      SINGLE PRECISION (QK15I-S, DQK15I-D)
  !***KEYWORDS  15-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
  !***AUTHOR  Piessens, Robert
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !           de Doncker, Elise
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !***DESCRIPTION
  !
  !           Integration Rule
  !           Standard Fortran subroutine
  !           Real version
  !
  !           PARAMETERS
  !            ON ENTRY
  !              F      - Real
  !                       Function subprogram defining the integrand
  !                       FUNCTION F(X). The actual name for F needs to be
  !                       Declared E X T E R N A L in the calling program.
  !
  !              BOUN   - Real
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
  !              A      - Real
  !                       Lower limit for integration over subrange
  !                       of (0,1)
  !
  !              B      - Real
  !                       Upper limit for integration over subrange
  !                       of (0,1)
  !
  !            ON RETURN
  !              RESULT - Real
  !                       Approximation to the integral I
  !                       Result is computed by applying the 15-POINT
  !                       KRONROD RULE(RESK) obtained by optimal addition
  !                       of abscissae to the 7-POINT GAUSS RULE(RESG).
  !
  !              ABSERR - Real
  !                       Estimate of the modulus of the absolute error,
  !                       WHICH SHOULD EQUAL or EXCEED ABS(I-RESULT)
  !
  !              RESABS - Real
  !                       Approximation to the integral J
  !
  !              RESASC - Real
  !                       Approximation to the integral of
  !                       ABS((TRANSFORMED INTEGRAND)-I/(B-A)) over (A,B)
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   800101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  QK15I
  !
  REAL A , absc , absc1 , absc2 , Abserr , B , Boun , centr , dinf , &
    R1MACH , epmach , F , fc , fsum , fval1 , fval2 , fv1 , fv2 , hlgth , &
    Resabs , Resasc , resg , resk , reskh , Result , tabsc1 , tabsc2 , &
    uflow , wg , wgk , xgk
  INTEGER Inf , j
  EXTERNAL F
  !
  DIMENSION fv1(7) , fv2(7) , xgk(8) , wgk(8) , wg(8)
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
  SAVE xgk , wgk , wg
  DATA xgk(1) , xgk(2) , xgk(3) , xgk(4) , xgk(5) , xgk(6) , xgk(7) , xgk(8)&
    /0.9914553711208126E+00 , 0.9491079123427585E+00 , &
    0.8648644233597691E+00 , 0.7415311855993944E+00 , &
    0.5860872354676911E+00 , 0.4058451513773972E+00 , &
    0.2077849550078985E+00 , 0.0000000000000000E+00/
  !
  DATA wgk(1) , wgk(2) , wgk(3) , wgk(4) , wgk(5) , wgk(6) , wgk(7) , wgk(8)&
    /0.2293532201052922E-01 , 0.6309209262997855E-01 , &
    0.1047900103222502E+00 , 0.1406532597155259E+00 , &
    0.1690047266392679E+00 , 0.1903505780647854E+00 , &
    0.2044329400752989E+00 , 0.2094821410847278E+00/
  !
  DATA wg(1) , wg(2) , wg(3) , wg(4) , wg(5) , wg(6) , wg(7) , &
    wg(8)/0.0000000000000000E+00 , 0.1294849661688697E+00 , &
    0.0000000000000000E+00 , 0.2797053914892767E+00 , &
    0.0000000000000000E+00 , 0.3818300505051189E+00 , &
    0.0000000000000000E+00 , 0.4179591836734694E+00/
  !
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
  !***FIRST EXECUTABLE STATEMENT  QK15I
  epmach = R1MACH(4)
  uflow = R1MACH(1)
  dinf = MIN(1,Inf)
  !
  centr = 0.5E+00*(A+B)
  hlgth = 0.5E+00*(B-A)
  tabsc1 = Boun + dinf*(0.1E+01-centr)/centr
  fval1 = F(tabsc1)
  IF ( Inf==2 ) fval1 = fval1 + F(-tabsc1)
  fc = (fval1/centr)/centr
  !
  !           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO
  !           THE INTEGRAL, AND ESTIMATE THE ERROR.
  !
  resg = wg(8)*fc
  resk = wgk(8)*fc
  Resabs = ABS(resk)
  DO j = 1 , 7
    absc = hlgth*xgk(j)
    absc1 = centr - absc
    absc2 = centr + absc
    tabsc1 = Boun + dinf*(0.1E+01-absc1)/absc1
    tabsc2 = Boun + dinf*(0.1E+01-absc2)/absc2
    fval1 = F(tabsc1)
    fval2 = F(tabsc2)
    IF ( Inf==2 ) fval1 = fval1 + F(-tabsc1)
    IF ( Inf==2 ) fval2 = fval2 + F(-tabsc2)
    fval1 = (fval1/absc1)/absc1
    fval2 = (fval2/absc2)/absc2
    fv1(j) = fval1
    fv2(j) = fval2
    fsum = fval1 + fval2
    resg = resg + wg(j)*fsum
    resk = resk + wgk(j)*fsum
    Resabs = Resabs + wgk(j)*(ABS(fval1)+ABS(fval2))
  ENDDO
  reskh = resk*0.5E+00
  Resasc = wgk(8)*ABS(fc-reskh)
  DO j = 1 , 7
    Resasc = Resasc + wgk(j)*(ABS(fv1(j)-reskh)+ABS(fv2(j)-reskh))
  ENDDO
  Result = resk*hlgth
  Resasc = Resasc*hlgth
  Resabs = Resabs*hlgth
  Abserr = ABS((resk-resg)*hlgth)
  IF ( Resasc/=0.0E+00.AND.Abserr/=0.E0 )&
    Abserr = Resasc*MIN(0.1E+01,(0.2E+03*Abserr/Resasc)**1.5E+00)
  IF ( Resabs>uflow/(0.5E+02*epmach) ) Abserr = MAX((epmach*0.5E+02)*Resabs,&
    Abserr)
END SUBROUTINE QK15I
