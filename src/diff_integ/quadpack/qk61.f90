!** QK61
SUBROUTINE QK61(F,A,B,Result,Abserr,Resabs,Resasc)
  !> To compute I = Integral of F over (A,B) with error
  !                           estimate
  !                       J = Integral of ABS(F) over (A,B)
  !***
  ! **Library:**   SLATEC (QUADPACK)
  !***
  ! **Category:**  H2A1A2
  !***
  ! **Type:**      SINGLE PRECISION (QK61-S, DQK61-D)
  !***
  ! **Keywords:**  61-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
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
  !        Integration rule
  !        Standard fortran subroutine
  !        Real version
  !
  !
  !        PARAMETERS
  !         ON ENTRY
  !           F      - Real
  !                    Function subprogram defining the integrand
  !                    function F(X). The actual name for F needs to be
  !                    declared E X T E R N A L in the calling program.
  !
  !           A      - Real
  !                    Lower limit of integration
  !
  !           B      - Real
  !                    Upper limit of integration
  !
  !         ON RETURN
  !           RESULT - Real
  !                    Approximation to the integral I
  !                    RESULT is computed by applying the 61-point
  !                    Kronrod rule (RESK) obtained by optimal addition of
  !                    abscissae to the 30-point Gauss rule (RESG).
  !
  !           ABSERR - Real
  !                    Estimate of the modulus of the absolute error,
  !                    which should equal or exceed ABS(I-RESULT)
  !
  !           RESABS - Real
  !                    Approximation to the integral J
  !
  !           RESASC - Real
  !                    Approximation to the integral of ABS(F-I/(B-A))
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
  USE service, ONLY : R1MACH
  !
  INTERFACE
    REAL(SP) FUNCTION F(X)
      IMPORT SP
      REAL(SP) :: X
    END FUNCTION F
  END INTERFACE
  REAL(SP) :: Abserr, B, A, Resabs, Resasc, Result
  INTEGER :: j, jtw, jtwm1
  REAL(SP) :: absc, centr, dhlgth, epmach, fc, fsum, fval1, fval2, fv1(30), &
    fv2(30), hlgth, resg, resk, reskh, uflow
  !
  !           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE
  !           INTERVAL (-1,1). BECAUSE OF SYMMETRY ONLY THE POSITIVE
  !           ABSCISSAE AND THEIR CORRESPONDING WEIGHTS ARE GIVEN.
  !
  !           XGK   - ABSCISSAE OF THE 61-POINT KRONROD RULE
  !                   XGK(2), XGK(4)  ... ABSCISSAE OF THE 30-POINT
  !                   GAUSS RULE
  !                   XGK(1), XGK(3)  ... OPTIMALLY ADDED ABSCISSAE
  !                   TO THE 30-POINT GAUSS RULE
  !
  !           WGK   - WEIGHTS OF THE 61-POINT KRONROD RULE
  !
  !           WG    - WEIGHTS OF THE 30-POINT GAUSS RULE
  !
  REAL(SP), PARAMETER :: xgk(31) = [ 0.9994844100504906_SP, 0.9968934840746495_SP, &
    0.9916309968704046_SP, 0.9836681232797472_SP, 0.9731163225011263_SP, &
    0.9600218649683075_SP, 0.9443744447485600_SP, 0.9262000474292743_SP, &
    0.9055733076999078_SP, 0.8825605357920527_SP, 0.8572052335460611_SP, &
    0.8295657623827684_SP, 0.7997278358218391_SP, 0.7677774321048262_SP, &
    0.7337900624532268_SP, 0.6978504947933158_SP, 0.6600610641266270_SP, &
    0.6205261829892429_SP, 0.5793452358263617_SP, 0.5366241481420199_SP, &
    0.4924804678617786_SP, 0.4470337695380892_SP, 0.4004012548303944_SP, &
    0.3527047255308781_SP, 0.3040732022736251_SP, 0.2546369261678898_SP, &
    0.2045251166823099_SP, 0.1538699136085835_SP, 0.1028069379667370_SP, &
    0.5147184255531770E-01_SP, 0.0_SP ]
  REAL(SP), PARAMETER :: wgk(31) = [ 0.1389013698677008E-02_SP, 0.3890461127099884E-02_SP, &
    0.6630703915931292E-02_SP, 0.9273279659517763E-02_SP, 0.1182301525349634E-01_SP, &
    0.1436972950704580E-01_SP, 0.1692088918905327E-01_SP, 0.1941414119394238E-01_SP, &
    0.2182803582160919E-01_SP, 0.2419116207808060E-01_SP, 0.2650995488233310E-01_SP, &
    0.2875404876504129E-01_SP, 0.3090725756238776E-01_SP, 0.3298144705748373E-01_SP, &
    0.3497933802806002E-01_SP, 0.3688236465182123E-01_SP, 0.3867894562472759E-01_SP, &
    0.4037453895153596E-01_SP, 0.4196981021516425E-01_SP, 0.4345253970135607E-01_SP, &
    0.4481480013316266E-01_SP, 0.4605923827100699E-01_SP, 0.4718554656929915E-01_SP, &
    0.4818586175708713E-01_SP, 0.4905543455502978E-01_SP, 0.4979568342707421E-01_SP, &
    0.5040592140278235E-01_SP, 0.5088179589874961E-01_SP, 0.5122154784925877E-01_SP, &
    0.5142612853745903E-01_SP, 0.5149472942945157E-01_SP ]
  REAL(SP), PARAMETER :: wg(15) = [ 0.7968192496166606E-02_SP, 0.1846646831109096E-01_SP, &
    0.2878470788332337E-01_SP, 0.3879919256962705E-01_SP, 0.4840267283059405E-01_SP, &
    0.5749315621761907E-01_SP, 0.6597422988218050E-01_SP, 0.7375597473770521E-01_SP, &
    0.8075589522942022E-01_SP, 0.8689978720108298E-01_SP, 0.9212252223778613E-01_SP, &
    0.9636873717464426E-01_SP, 0.9959342058679527E-01_SP, 0.1017623897484055E+00_SP, &
    0.1028526528935588E+00_SP ]
  !
  !           LIST OF MAJOR VARIABLES
  !           -----------------------
  !
  !           CENTR  - MID POINT OF THE INTERVAL
  !           HLGTH  - HALF-LENGTH OF THE INTERVAL
  !           ABSC   - ABSCISSA
  !           FVAL*  - FUNCTION VALUE
  !           RESG   - RESULT OF THE 30-POINT GAUSS RULE
  !           RESK   - RESULT OF THE 61-POINT KRONROD RULE
  !           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F
  !                    OVER (A,B), I.E. TO I/(B-A)
  !
  !           MACHINE DEPENDENT CONSTANTS
  !           ---------------------------
  !
  !           EPMACH IS THE LARGEST RELATIVE SPACING.
  !           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
  !
  !* FIRST EXECUTABLE STATEMENT  QK61
  epmach = R1MACH(4)
  uflow = R1MACH(1)
  !
  centr = 0.5_SP*(B+A)
  hlgth = 0.5_SP*(B-A)
  dhlgth = ABS(hlgth)
  !
  !           COMPUTE THE 61-POINT KRONROD APPROXIMATION TO THE
  !           INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
  !
  resg = 0._SP
  fc = F(centr)
  resk = wgk(31)*fc
  Resabs = ABS(resk)
  DO j = 1, 15
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
  END DO
  DO j = 1, 15
    jtwm1 = j*2 - 1
    absc = hlgth*xgk(jtwm1)
    fval1 = F(centr-absc)
    fval2 = F(centr+absc)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1 + fval2
    resk = resk + wgk(jtwm1)*fsum
    Resabs = Resabs + wgk(jtwm1)*(ABS(fval1)+ABS(fval2))
  END DO
  reskh = resk*0.5_SP
  Resasc = wgk(31)*ABS(fc-reskh)
  DO j = 1, 30
    Resasc = Resasc + wgk(j)*(ABS(fv1(j)-reskh)+ABS(fv2(j)-reskh))
  END DO
  Result = resk*hlgth
  Resabs = Resabs*dhlgth
  Resasc = Resasc*dhlgth
  Abserr = ABS((resk-resg)*hlgth)
  IF( Resasc/=0._SP .AND. Abserr/=0._SP )&
    Abserr = Resasc*MIN(1._SP,(200._SP*Abserr/Resasc)**1.5_SP)
  IF( Resabs>uflow/(50._SP*epmach) ) Abserr = MAX((epmach*50._SP)*Resabs,&
    Abserr)
END SUBROUTINE QK61
