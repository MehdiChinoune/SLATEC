!** DQK41
PURE SUBROUTINE DQK41(F,A,B,Result,Abserr,Resabs,Resasc)
  !> To compute I = Integral of F over (A,B) with error estimate
  !  J = Integral of ABS(F) over (A,B)
  !***
  ! **Library:**   SLATEC (QUADPACK)
  !***
  ! **Category:**  H2A1A2
  !***
  ! **Type:**      DOUBLE PRECISION (QK41-S, DQK41-D)
  !***
  ! **Keywords:**  41-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
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
  !           Double precision version
  !
  !           PARAMETERS
  !            ON ENTRY
  !              F      - Double precision
  !                       Function subprogram defining the integrand
  !                       FUNCTION F(X). The actual name for F needs to be
  !                       declared E X T E R N A L in the calling program.
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
  !                       RESULT is computed by applying the 41-POINT
  !                       GAUSS-KRONROD RULE (RESK) obtained by optimal
  !                       addition of abscissae to the 20-POINT GAUSS
  !                       RULE (RESG).
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
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   800101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : tiny_dp, eps_dp
  !
  INTERFACE
    REAL(DP) PURE FUNCTION F(X)
      IMPORT DP
      REAL(DP), INTENT(IN) :: X
    END FUNCTION F
  END INTERFACE
  REAL(DP), INTENT(IN) :: A, B
  REAL(DP), INTENT(OUT) :: Abserr, Resabs, Resasc, Result
  !
  INTEGER :: j, jtw, jtwm1
  REAL(DP) :: absc, centr, dhlgth, epmach, fc, fsum, fval1, fval2, fv1(20), &
    fv2(20), hlgth, resg, resk, reskh, uflow
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
  !
  ! GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS
  ! AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
  ! BELL LABS, NOV. 1981.
  !
  REAL(DP), PARAMETER :: wg(10) = [ 0.017614007139152118311861962351853_DP, &
    0.040601429800386941331039952274932_DP, 0.062672048334109063569506535187042_DP, &
    0.083276741576704748724758143222046_DP, 0.101930119817240435036750135480350_DP, &
    0.118194531961518417312377377711382_DP, 0.131688638449176626898494499748163_DP, &
    0.142096109318382051329298325067165_DP, 0.149172986472603746787828737001969_DP, &
    0.152753387130725850698084331955098_DP ]
  REAL(DP), PARAMETER :: xgk(21) = [ 0.998859031588277663838315576545863_DP, &
    0.993128599185094924786122388471320_DP, 0.981507877450250259193342994720217_DP, &
    0.963971927277913791267666131197277_DP, 0.940822633831754753519982722212443_DP, &
    0.912234428251325905867752441203298_DP, 0.878276811252281976077442995113078_DP, &
    0.839116971822218823394529061701521_DP, 0.795041428837551198350638833272788_DP, &
    0.746331906460150792614305070355642_DP, 0.693237656334751384805490711845932_DP, &
    0.636053680726515025452836696226286_DP, 0.575140446819710315342946036586425_DP, &
    0.510867001950827098004364050955251_DP, 0.443593175238725103199992213492640_DP, &
    0.373706088715419560672548177024927_DP, 0.301627868114913004320555356858592_DP, &
    0.227785851141645078080496195368575_DP, 0.152605465240922675505220241022678_DP, &
    0.076526521133497333754640409398838_DP, 0.000000000000000000000000000000000_DP ]
  REAL(DP), PARAMETER :: wgk(21) = [ 0.003073583718520531501218293246031_DP, &
    0.008600269855642942198661787950102_DP, 0.014626169256971252983787960308868_DP, &
    0.020388373461266523598010231432755_DP, 0.025882133604951158834505067096153_DP, &
    0.031287306777032798958543119323801_DP, 0.036600169758200798030557240707211_DP, &
    0.041668873327973686263788305936895_DP, 0.046434821867497674720231880926108_DP, &
    0.050944573923728691932707670050345_DP, 0.055195105348285994744832372419777_DP, &
    0.059111400880639572374967220648594_DP, 0.062653237554781168025870122174255_DP, &
    0.065834597133618422111563556969398_DP, 0.068648672928521619345623411885368_DP, &
    0.071054423553444068305790361723210_DP, 0.073030690332786667495189417658913_DP, &
    0.074582875400499188986581418362488_DP, 0.075704497684556674659542775376617_DP, &
    0.076377867672080736705502835038061_DP, 0.076600711917999656445049901530102_DP ]
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
  !* FIRST EXECUTABLE STATEMENT  DQK41
  epmach = eps_dp
  uflow = tiny_dp
  !
  centr = 0.5_DP*(A+B)
  hlgth = 0.5_DP*(B-A)
  dhlgth = ABS(hlgth)
  !
  !           COMPUTE THE 41-POINT GAUSS-KRONROD APPROXIMATION TO
  !           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
  !
  resg = 0._DP
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
  END DO
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
  END DO
  reskh = resk*0.5_DP
  Resasc = wgk(21)*ABS(fc-reskh)
  DO j = 1, 20
    Resasc = Resasc + wgk(j)*(ABS(fv1(j)-reskh)+ABS(fv2(j)-reskh))
  END DO
  Result = resk*hlgth
  Resabs = Resabs*dhlgth
  Resasc = Resasc*dhlgth
  Abserr = ABS((resk-resg)*hlgth)
  IF( Resasc/=0._DP .AND. Abserr/=0.D+00 )&
    Abserr = Resasc*MIN(1._DP,(0.2E+03_DP*Abserr/Resasc)**1.5_DP)
  IF( Resabs>uflow/(0.5E+02_DP*epmach) ) Abserr = MAX((epmach*0.5E+02_DP)*Resabs,Abserr)
  !
END SUBROUTINE DQK41