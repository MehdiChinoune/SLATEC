!DECK DQK51
SUBROUTINE DQK51(F,A,B,Result,Abserr,Resabs,Resasc)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DQK51
  !***PURPOSE  To compute I = Integral of F over (A,B) with error
  !                           estimate
  !                       J = Integral of ABS(F) over (A,B)
  !***LIBRARY   SLATEC (QUADPACK)
  !***CATEGORY  H2A1A2
  !***TYPE      DOUBLE PRECISION (QK51-S, DQK51-D)
  !***KEYWORDS  51-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
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
  !                       Function subroutine defining the integrand
  !                       function F(X). The actual name for F needs to be
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
  !                       RESULT is computed by applying the 51-point
  !                       Kronrod rule (RESK) obtained by optimal addition
  !                       of abscissae to the 25-point Gauss rule (RESG).
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
  !   910819  Added WGK(26) to code.  (WRB)
  !***END PROLOGUE  DQK51
  !
  REAL(8) :: A, absc, Abserr, B, centr, dhlgth, D1MACH, &
    epmach, F, fc, fsum, fval1, fval2, fv1, fv2, &
    hlgth, Resabs, Resasc, resg, resk, reskh, Result, &
    uflow, wg, wgk, xgk
  INTEGER j, jtw, jtwm1
  EXTERNAL F
  !
  DIMENSION fv1(25), fv2(25), xgk(26), wgk(26), wg(13)
  !
  !           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
  !           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
  !           CORRESPONDING WEIGHTS ARE GIVEN.
  !
  !           XGK    - ABSCISSAE OF THE 51-POINT KRONROD RULE
  !                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 25-POINT
  !                    GAUSS RULE
  !                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
  !                    ADDED TO THE 25-POINT GAUSS RULE
  !
  !           WGK    - WEIGHTS OF THE 51-POINT KRONROD RULE
  !
  !           WG     - WEIGHTS OF THE 25-POINT GAUSS RULE
  !
  !
  ! GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS
  ! AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
  ! BELL LABS, NOV. 1981.
  !
  SAVE wg, xgk, wgk
  DATA wg(1)/0.011393798501026287947902964113235D0/
  DATA wg(2)/0.026354986615032137261901815295299D0/
  DATA wg(3)/0.040939156701306312655623487711646D0/
  DATA wg(4)/0.054904695975835191925936891540473D0/
  DATA wg(5)/0.068038333812356917207187185656708D0/
  DATA wg(6)/0.080140700335001018013234959669111D0/
  DATA wg(7)/0.091028261982963649811497220702892D0/
  DATA wg(8)/0.100535949067050644202206890392686D0/
  DATA wg(9)/0.108519624474263653116093957050117D0/
  DATA wg(10)/0.114858259145711648339325545869556D0/
  DATA wg(11)/0.119455763535784772228178126512901D0/
  DATA wg(12)/0.122242442990310041688959518945852D0/
  DATA wg(13)/0.123176053726715451203902873079050D0/
  !
  DATA xgk(1)/0.999262104992609834193457486540341D0/
  DATA xgk(2)/0.995556969790498097908784946893902D0/
  DATA xgk(3)/0.988035794534077247637331014577406D0/
  DATA xgk(4)/0.976663921459517511498315386479594D0/
  DATA xgk(5)/0.961614986425842512418130033660167D0/
  DATA xgk(6)/0.942974571228974339414011169658471D0/
  DATA xgk(7)/0.920747115281701561746346084546331D0/
  DATA xgk(8)/0.894991997878275368851042006782805D0/
  DATA xgk(9)/0.865847065293275595448996969588340D0/
  DATA xgk(10)/0.833442628760834001421021108693570D0/
  DATA xgk(11)/0.797873797998500059410410904994307D0/
  DATA xgk(12)/0.759259263037357630577282865204361D0/
  DATA xgk(13)/0.717766406813084388186654079773298D0/
  DATA xgk(14)/0.673566368473468364485120633247622D0/
  DATA xgk(15)/0.626810099010317412788122681624518D0/
  DATA xgk(16)/0.577662930241222967723689841612654D0/
  DATA xgk(17)/0.526325284334719182599623778158010D0/
  DATA xgk(18)/0.473002731445714960522182115009192D0/
  DATA xgk(19)/0.417885382193037748851814394594572D0/
  DATA xgk(20)/0.361172305809387837735821730127641D0/
  DATA xgk(21)/0.303089538931107830167478909980339D0/
  DATA xgk(22)/0.243866883720988432045190362797452D0/
  DATA xgk(23)/0.183718939421048892015969888759528D0/
  DATA xgk(24)/0.122864692610710396387359818808037D0/
  DATA xgk(25)/0.061544483005685078886546392366797D0/
  DATA xgk(26)/0.000000000000000000000000000000000D0/
  !
  DATA wgk(1)/0.001987383892330315926507851882843D0/
  DATA wgk(2)/0.005561932135356713758040236901066D0/
  DATA wgk(3)/0.009473973386174151607207710523655D0/
  DATA wgk(4)/0.013236229195571674813656405846976D0/
  DATA wgk(5)/0.016847817709128298231516667536336D0/
  DATA wgk(6)/0.020435371145882835456568292235939D0/
  DATA wgk(7)/0.024009945606953216220092489164881D0/
  DATA wgk(8)/0.027475317587851737802948455517811D0/
  DATA wgk(9)/0.030792300167387488891109020215229D0/
  DATA wgk(10)/0.034002130274329337836748795229551D0/
  DATA wgk(11)/0.037116271483415543560330625367620D0/
  DATA wgk(12)/0.040083825504032382074839284467076D0/
  DATA wgk(13)/0.042872845020170049476895792439495D0/
  DATA wgk(14)/0.045502913049921788909870584752660D0/
  DATA wgk(15)/0.047982537138836713906392255756915D0/
  DATA wgk(16)/0.050277679080715671963325259433440D0/
  DATA wgk(17)/0.052362885806407475864366712137873D0/
  DATA wgk(18)/0.054251129888545490144543370459876D0/
  DATA wgk(19)/0.055950811220412317308240686382747D0/
  DATA wgk(20)/0.057437116361567832853582693939506D0/
  DATA wgk(21)/0.058689680022394207961974175856788D0/
  DATA wgk(22)/0.059720340324174059979099291932562D0/
  DATA wgk(23)/0.060539455376045862945360267517565D0/
  DATA wgk(24)/0.061128509717053048305859030416293D0/
  DATA wgk(25)/0.061471189871425316661544131965264D0/
  DATA wgk(26)/0.061580818067832935078759824240055D0/
  !
  !
  !           LIST OF MAJOR VARIABLES
  !           -----------------------
  !
  !           CENTR  - MID POINT OF THE INTERVAL
  !           HLGTH  - HALF-LENGTH OF THE INTERVAL
  !           ABSC   - ABSCISSA
  !           FVAL*  - FUNCTION VALUE
  !           RESG   - RESULT OF THE 25-POINT GAUSS FORMULA
  !           RESK   - RESULT OF THE 51-POINT KRONROD FORMULA
  !           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
  !                    I.E. TO I/(B-A)
  !
  !           MACHINE DEPENDENT CONSTANTS
  !           ---------------------------
  !
  !           EPMACH IS THE LARGEST RELATIVE SPACING.
  !           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
  !
  !***FIRST EXECUTABLE STATEMENT  DQK51
  epmach = D1MACH(4)
  uflow = D1MACH(1)
  !
  centr = 0.5D+00*(A+B)
  hlgth = 0.5D+00*(B-A)
  dhlgth = ABS(hlgth)
  !
  !           COMPUTE THE 51-POINT KRONROD APPROXIMATION TO
  !           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
  !
  fc = F(centr)
  resg = wg(13)*fc
  resk = wgk(26)*fc
  Resabs = ABS(resk)
  DO j = 1, 12
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
  DO j = 1, 13
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
  Resasc = wgk(26)*ABS(fc-reskh)
  DO j = 1, 25
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
END SUBROUTINE DQK51
