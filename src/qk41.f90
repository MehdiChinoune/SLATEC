!*==QK41.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK QK41
SUBROUTINE QK41(F,A,B,Result,Abserr,Resabs,Resasc)
  IMPLICIT NONE
  !*--QK415
  !***BEGIN PROLOGUE  QK41
  !***PURPOSE  To compute I = Integral of F over (A,B), with error
  !                           estimate
  !                       J = Integral of ABS(F) over (A,B)
  !***LIBRARY   SLATEC (QUADPACK)
  !***CATEGORY  H2A1A2
  !***TYPE      SINGLE PRECISION (QK41-S, DQK41-D)
  !***KEYWORDS  41-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
  !***AUTHOR  Piessens, Robert
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !           de Doncker, Elise
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !***DESCRIPTION
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
  !                       declared E X T E R N A L in the calling program.
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
  !                       RESULT is computed by applying the 41-POINT
  !                       GAUSS-KRONROD RULE (RESK) obtained by optimal
  !                       addition of abscissae to the 20-POINT GAUSS
  !                       RULE (RESG).
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
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   800101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  QK41
  !
  REAL A, absc, Abserr, B, centr, dhlgth, epmach, F, fc, fsum, &
    fval1, fval2, fv1, fv2, hlgth, Resabs, Resasc, resg, resk, &
    reskh, Result, R1MACH, uflow, wg, wgk, xgk
  INTEGER j, jtw, jtwm1
  EXTERNAL F
  !
  DIMENSION fv1(20), fv2(20), xgk(21), wgk(21), wg(10)
  !
  !           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
  !           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
  !           CORRESPONDING WEIGHTS ARE GIVEN.
  !
  !           XGK    - ABSCISSAE OF THE 41-POINT GAUSS-KRONROD RULE
  !                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 20-POINT
  !                    GAUSS RULE
  !                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
  !                    ADDED TO THE 20-POINT GAUSS RULE
  !
  !           WGK    - WEIGHTS OF THE 41-POINT GAUSS-KRONROD RULE
  !
  !           WG     - WEIGHTS OF THE 20-POINT GAUSS RULE
  !
  SAVE xgk, wgk, wg
  DATA xgk(1), xgk(2), xgk(3), xgk(4), xgk(5), xgk(6), xgk(7), xgk(8)&
    , xgk(9), xgk(10), xgk(11), xgk(12), xgk(13), xgk(14), xgk(15)&
    , xgk(16), xgk(17), xgk(18), xgk(19), xgk(20), xgk(21)&
    /0.9988590315882777E+00, 0.9931285991850949E+00, &
    0.9815078774502503E+00, 0.9639719272779138E+00, &
    0.9408226338317548E+00, 0.9122344282513259E+00, &
    0.8782768112522820E+00, 0.8391169718222188E+00, &
    0.7950414288375512E+00, 0.7463319064601508E+00, &
    0.6932376563347514E+00, 0.6360536807265150E+00, &
    0.5751404468197103E+00, 0.5108670019508271E+00, &
    0.4435931752387251E+00, 0.3737060887154196E+00, &
    0.3016278681149130E+00, 0.2277858511416451E+00, &
    0.1526054652409227E+00, 0.7652652113349733E-01, 0.0E+00/
  DATA wgk(1), wgk(2), wgk(3), wgk(4), wgk(5), wgk(6), wgk(7), wgk(8)&
    , wgk(9), wgk(10), wgk(11), wgk(12), wgk(13), wgk(14), wgk(15)&
    , wgk(16), wgk(17), wgk(18), wgk(19), wgk(20), wgk(21)&
    /0.3073583718520532E-02, 0.8600269855642942E-02, &
    0.1462616925697125E-01, 0.2038837346126652E-01, &
    0.2588213360495116E-01, 0.3128730677703280E-01, &
    0.3660016975820080E-01, 0.4166887332797369E-01, &
    0.4643482186749767E-01, 0.5094457392372869E-01, &
    0.5519510534828599E-01, 0.5911140088063957E-01, &
    0.6265323755478117E-01, 0.6583459713361842E-01, &
    0.6864867292852162E-01, 0.7105442355344407E-01, &
    0.7303069033278667E-01, 0.7458287540049919E-01, &
    0.7570449768455667E-01, 0.7637786767208074E-01, &
    0.7660071191799966E-01/
  DATA wg(1), wg(2), wg(3), wg(4), wg(5), wg(6), wg(7), wg(8), &
    wg(9), wg(10)/0.1761400713915212E-01, 0.4060142980038694E-01, &
    0.6267204833410906E-01, 0.8327674157670475E-01, &
    0.1019301198172404E+00, 0.1181945319615184E+00, &
    0.1316886384491766E+00, 0.1420961093183821E+00, &
    0.1491729864726037E+00, 0.1527533871307259E+00/
  !
  !
  !           LIST OF MAJOR VARIABLES
  !           -----------------------
  !
  !           CENTR  - MID POINT OF THE INTERVAL
  !           HLGTH  - HALF-LENGTH OF THE INTERVAL
  !           ABSC   - ABSCISSA
  !           FVAL*  - FUNCTION VALUE
  !           RESG   - RESULT OF THE 20-POINT GAUSS FORMULA
  !           RESK   - RESULT OF THE 41-POINT KRONROD FORMULA
  !           RESKH  - APPROXIMATION TO MEAN VALUE OF F OVER (A,B), I.E.
  !                    TO I/(B-A)
  !
  !           MACHINE DEPENDENT CONSTANTS
  !           ---------------------------
  !
  !           EPMACH IS THE LARGEST RELATIVE SPACING.
  !           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
  !
  !***FIRST EXECUTABLE STATEMENT  QK41
  epmach = R1MACH(4)
  uflow = R1MACH(1)
  !
  centr = 0.5E+00*(A+B)
  hlgth = 0.5E+00*(B-A)
  dhlgth = ABS(hlgth)
  !
  !           COMPUTE THE 41-POINT GAUSS-KRONROD APPROXIMATION TO
  !           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
  !
  resg = 0.0E+00
  fc = F(centr)
  resk = wgk(21)*fc
  Resabs = ABS(resk)
  DO j = 1, 10
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
  DO j = 1, 10
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
  Resasc = wgk(21)*ABS(fc-reskh)
  DO j = 1, 20
    Resasc = Resasc + wgk(j)*(ABS(fv1(j)-reskh)+ABS(fv2(j)-reskh))
  ENDDO
  Result = resk*hlgth
  Resabs = Resabs*dhlgth
  Resasc = Resasc*dhlgth
  Abserr = ABS((resk-resg)*hlgth)
  IF ( Resasc/=0.0E+00.AND.Abserr/=0.E+00 )&
    Abserr = Resasc*MIN(0.1E+01,(0.2E+03*Abserr/Resasc)**1.5E+00)
  IF ( Resabs>uflow/(0.5E+02*epmach) ) Abserr = MAX((epmach*0.5E+02)*Resabs,&
    Abserr)
END SUBROUTINE QK41
