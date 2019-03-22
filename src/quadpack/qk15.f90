!** QK15
SUBROUTINE QK15(F,A,B,Result,Abserr,Resabs,Resasc)
  IMPLICIT NONE
  !>
  !***
  !  To compute I = Integral of F over (A,B), with error
  !                           estimate
  !                       J = integral of ABS(F) over (A,B)
  !***
  ! **Library:**   SLATEC (QUADPACK)
  !***
  ! **Category:**  H2A1A2
  !***
  ! **Type:**      SINGLE PRECISION (QK15-S, DQK15-D)
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
  !           Integration rules
  !           Standard fortran subroutine
  !           Real version
  !
  !           PARAMETERS
  !            ON ENTRY
  !              F      - Real
  !                       Function subprogram defining the integrand
  !                       FUNCTION F(X). The actual name for F needs to be
  !                       Declared E X T E R N A L in the calling program.
  !
  !              A      - Real
  !                       Lower limit of integration
  !
  !              B      - Real
  !                       Upper limit of integration
  !
  !            ON RETURN
  !              RESULT - Real
  !                       Approximation to the integral I
  !                       Result is computed by applying the 15-POINT
  !                       KRONROD RULE (RESK) obtained by optimal addition
  !                       of abscissae to the 7-POINT GAUSS RULE(RESG).
  !
  !              ABSERR - Real
  !                       Estimate of the modulus of the absolute error,
  !                       which should not exceed ABS(I-RESULT)
  !
  !              RESABS - Real
  !                       Approximation to the integral J
  !
  !              RESASC - Real
  !                       Approximation to the integral of ABS(F-I/(B-A))
  !                       over (A,B)
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  R1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   800101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)

  !
  REAL A, absc, Abserr, B, centr, dhlgth, epmach, fc, fsum, &
    fval1, fval2, fv1, fv2, hlgth, Resabs, Resasc, resg, resk, &
    reskh, Result, uflow, wg, wgk, xgk
  INTEGER j, jtw, jtwm1
  REAL, EXTERNAL :: F
  REAL, EXTERNAL :: R1MACH
  !
  DIMENSION fv1(7), fv2(7), wg(4), wgk(8), xgk(8)
  !
  !           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
  !           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
  !           CORRESPONDING WEIGHTS ARE GIVEN.
  !
  !           XGK    - ABSCISSAE OF THE 15-POINT KRONROD RULE
  !                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 7-POINT
  !                    GAUSS RULE
  !                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
  !                    ADDED TO THE 7-POINT GAUSS RULE
  !
  !           WGK    - WEIGHTS OF THE 15-POINT KRONROD RULE
  !
  !           WG     - WEIGHTS OF THE 7-POINT GAUSS RULE
  !
  SAVE xgk, wgk, wg
  DATA xgk(1), xgk(2), xgk(3), xgk(4), xgk(5), xgk(6), xgk(7), xgk(8)&
    /0.9914553711208126E+00, 0.9491079123427585E+00, &
    0.8648644233597691E+00, 0.7415311855993944E+00, &
    0.5860872354676911E+00, 0.4058451513773972E+00, &
    0.2077849550078985E+00, 0.0E+00/
  DATA wgk(1), wgk(2), wgk(3), wgk(4), wgk(5), wgk(6), wgk(7), wgk(8)&
    /0.2293532201052922E-01, 0.6309209262997855E-01, &
    0.1047900103222502E+00, 0.1406532597155259E+00, &
    0.1690047266392679E+00, 0.1903505780647854E+00, &
    0.2044329400752989E+00, 0.2094821410847278E+00/
  DATA wg(1), wg(2), wg(3), wg(4)/0.1294849661688697E+00, &
    0.2797053914892767E+00, 0.3818300505051189E+00, &
    0.4179591836734694E+00/
  !
  !
  !           LIST OF MAJOR VARIABLES
  !           -----------------------
  !
  !           CENTR  - MID POINT OF THE INTERVAL
  !           HLGTH  - HALF-LENGTH OF THE INTERVAL
  !           ABSC   - ABSCISSA
  !           FVAL*  - FUNCTION VALUE
  !           RESG   - RESULT OF THE 7-POINT GAUSS FORMULA
  !           RESK   - RESULT OF THE 15-POINT KRONROD FORMULA
  !           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
  !                    I.E. TO I/(B-A)
  !
  !           MACHINE DEPENDENT CONSTANTS
  !           ---------------------------
  !
  !           EPMACH IS THE LARGEST RELATIVE SPACING.
  !           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
  !
  !* FIRST EXECUTABLE STATEMENT  QK15
  epmach = R1MACH(4)
  uflow = R1MACH(1)
  !
  centr = 0.5E+00*(A+B)
  hlgth = 0.5E+00*(B-A)
  dhlgth = ABS(hlgth)
  !
  !           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO
  !           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
  !
  fc = F(centr)
  resg = fc*wg(4)
  resk = fc*wgk(8)
  Resabs = ABS(resk)
  DO j = 1, 3
    jtw = j*2
    absc = hlgth*xgk(jtw)
    fval1 = F(centr-absc)
    fval2 = F(centr+absc)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1 + fval2
    resg = resg + wg(j)*fsum
    resk = resk + wgk(jtw)*fsum
    Resabs = Resabs + wgk(jtw)*(ABS(fval1)+ABS(fval2))
  ENDDO
  DO j = 1, 4
    jtwm1 = j*2 - 1
    absc = hlgth*xgk(jtwm1)
    fval1 = F(centr-absc)
    fval2 = F(centr+absc)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1 + fval2
    resk = resk + wgk(jtwm1)*fsum
    Resabs = Resabs + wgk(jtwm1)*(ABS(fval1)+ABS(fval2))
  ENDDO
  reskh = resk*0.5E+00
  Resasc = wgk(8)*ABS(fc-reskh)
  DO j = 1, 7
    Resasc = Resasc + wgk(j)*(ABS(fv1(j)-reskh)+ABS(fv2(j)-reskh))
  ENDDO
  Result = resk*hlgth
  Resabs = Resabs*dhlgth
  Resasc = Resasc*dhlgth
  Abserr = ABS((resk-resg)*hlgth)
  IF ( Resasc/=0.0E+00.AND.Abserr/=0.0E+00 )&
    Abserr = Resasc*MIN(0.1E+01,(0.2E+03*Abserr/Resasc)**1.5E+00)
  IF ( Resabs>uflow/(0.5E+02*epmach) ) Abserr = MAX((epmach*0.5E+02)*Resabs,&
    Abserr)
END SUBROUTINE QK15
