!*==DQK15.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DQK15
SUBROUTINE DQK15(F,A,B,Result,Abserr,Resabs,Resasc)
  IMPLICIT NONE
  !*--DQK155
  !***BEGIN PROLOGUE  DQK15
  !***PURPOSE  To compute I = Integral of F over (A,B), with error
  !                           estimate
  !                       J = integral of ABS(F) over (A,B)
  !***LIBRARY   SLATEC (QUADPACK)
  !***CATEGORY  H2A1A2
  !***TYPE      DOUBLE PRECISION (QK15-S, DQK15-D)
  !***KEYWORDS  15-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
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
  !           Double precision version
  !
  !           PARAMETERS
  !            ON ENTRY
  !              F      - Double precision
  !                       Function subprogram defining the integrand
  !                       FUNCTION F(X). The actual name for F needs to be
  !                       Declared E X T E R N A L in the calling program.
  !
  !              A      - Double precision
  !                       Lower limit of integration
  !
  !              B      - Double precision
  !                       Upper limit of integration
  !
  !            ON RETURN
  !              RESULT - Double precision
  !                       Approximation to the integral I
  !                       Result is computed by applying the 15-POINT
  !                       KRONROD RULE (RESK) obtained by optimal addition
  !                       of abscissae to the 7-POINT GAUSS RULE(RESG).
  !
  !              ABSERR - Double precision
  !                       Estimate of the modulus of the absolute error,
  !                       which should not exceed ABS(I-RESULT)
  !
  !              RESABS - Double precision
  !                       Approximation to the integral J
  !
  !              RESASC - Double precision
  !                       Approximation to the integral of ABS(F-I/(B-A))
  !                       over (A,B)
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   800101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DQK15
  !
  REAL(8) :: A , absc , Abserr , B , centr , dhlgth , D1MACH , &
    epmach , F , fc , fsum , fval1 , fval2 , fv1 , fv2 , &
    hlgth , Resabs , Resasc , resg , resk , reskh , Result , &
    uflow , wg , wgk , xgk
  INTEGER j , jtw , jtwm1
  EXTERNAL F
  !
  DIMENSION fv1(7) , fv2(7) , wg(4) , wgk(8) , xgk(8)
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
  !
  ! GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS
  ! AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
  ! BELL LABS, NOV. 1981.
  !
  SAVE wg , xgk , wgk
  DATA wg(1)/0.129484966168869693270611432679082D0/
  DATA wg(2)/0.279705391489276667901467771423780D0/
  DATA wg(3)/0.381830050505118944950369775488975D0/
  DATA wg(4)/0.417959183673469387755102040816327D0/
  !
  DATA xgk(1)/0.991455371120812639206854697526329D0/
  DATA xgk(2)/0.949107912342758524526189684047851D0/
  DATA xgk(3)/0.864864423359769072789712788640926D0/
  DATA xgk(4)/0.741531185599394439863864773280788D0/
  DATA xgk(5)/0.586087235467691130294144838258730D0/
  DATA xgk(6)/0.405845151377397166906606412076961D0/
  DATA xgk(7)/0.207784955007898467600689403773245D0/
  DATA xgk(8)/0.000000000000000000000000000000000D0/
  !
  DATA wgk(1)/0.022935322010529224963732008058970D0/
  DATA wgk(2)/0.063092092629978553290700663189204D0/
  DATA wgk(3)/0.104790010322250183839876322541518D0/
  DATA wgk(4)/0.140653259715525918745189590510238D0/
  DATA wgk(5)/0.169004726639267902826583426598550D0/
  DATA wgk(6)/0.190350578064785409913256402421014D0/
  DATA wgk(7)/0.204432940075298892414161999234649D0/
  DATA wgk(8)/0.209482141084727828012999174891714D0/
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
  !***FIRST EXECUTABLE STATEMENT  DQK15
  epmach = D1MACH(4)
  uflow = D1MACH(1)
  !
  centr = 0.5D+00*(A+B)
  hlgth = 0.5D+00*(B-A)
  dhlgth = ABS(hlgth)
  !
  !           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO
  !           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
  !
  fc = F(centr)
  resg = fc*wg(4)
  resk = fc*wgk(8)
  Resabs = ABS(resk)
  DO j = 1 , 3
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
  DO j = 1 , 4
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
  reskh = resk*0.5D+00
  Resasc = wgk(8)*ABS(fc-reskh)
  DO j = 1 , 7
    Resasc = Resasc + wgk(j)*(ABS(fv1(j)-reskh)+ABS(fv2(j)-reskh))
  ENDDO
  Result = resk*hlgth
  Resabs = Resabs*dhlgth
  Resasc = Resasc*dhlgth
  Abserr = ABS((resk-resg)*hlgth)
  IF ( Resasc/=0.0D+00.AND.Abserr/=0.0D+00 )&
    Abserr = Resasc*MIN(0.1D+01,(0.2D+03*Abserr/Resasc)**1.5D+00)
  IF ( Resabs>uflow/(0.5D+02*epmach) ) Abserr = MAX((epmach*0.5D+02)*Resabs,&
    Abserr)
END SUBROUTINE DQK15
