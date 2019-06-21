!** DQK61
SUBROUTINE DQK61(F,A,B,Result,Abserr,Resabs,Resasc)
  !> To compute I = Integral of F over (A,B) with error
  !                           estimate
  !                       J = Integral of ABS(F) over (A,B)
  !***
  ! **Library:**   SLATEC (QUADPACK)
  !***
  ! **Category:**  H2A1A2
  !***
  ! **Type:**      DOUBLE PRECISION (QK61-S, DQK61-D)
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
  !        Double precision version
  !
  !
  !        PARAMETERS
  !         ON ENTRY
  !           F      - Double precision
  !                    Function subprogram defining the integrand
  !                    function F(X). The actual name for F needs to be
  !                    declared E X T E R N A L in the calling program.
  !
  !           A      - Double precision
  !                    Lower limit of integration
  !
  !           B      - Double precision
  !                    Upper limit of integration
  !
  !         ON RETURN
  !           RESULT - Double precision
  !                    Approximation to the integral I
  !                    RESULT is computed by applying the 61-point
  !                    Kronrod rule (RESK) obtained by optimal addition of
  !                    abscissae to the 30-point Gauss rule (RESG).
  !
  !           ABSERR - Double precision
  !                    Estimate of the modulus of the absolute error,
  !                    which should equal or exceed ABS(I-RESULT)
  !
  !           RESABS - Double precision
  !                    Approximation to the integral J
  !
  !           RESASC - Double precision
  !                    Approximation to the integral of ABS(F-I/(B-A))
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
    REAL(DP) FUNCTION F(X)
      IMPORT DP
      REAL(DP) :: X
    END FUNCTION F
  END INTERFACE
  REAL(DP) :: A, Abserr, B, Resabs, Resasc, Result
  INTEGER :: j, jtw, jtwm1
  REAL(DP) :: dabsc, centr, dhlgth, epmach, fc, fsum, fval1, fval2, fv1(30), &
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
  !
  ! GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS
  ! AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
  ! BELL LABS, NOV. 1981.
  !
  REAL(DP), PARAMETER :: wg(15) = [ 0.007968192496166605615465883474674_DP, &
    0.018466468311090959142302131912047_DP, 0.028784707883323369349719179611292_DP, &
    0.038799192569627049596801936446348_DP, 0.048402672830594052902938140422808_DP, &
    0.057493156217619066481721689402056_DP, 0.065974229882180495128128515115962_DP, &
    0.073755974737705206268243850022191_DP, 0.080755895229420215354694938460530_DP, &
    0.086899787201082979802387530715126_DP, 0.092122522237786128717632707087619_DP, &
    0.096368737174644259639468626351810_DP, 0.099593420586795267062780282103569_DP, &
    0.101762389748405504596428952168554_DP, 0.102852652893558840341285636705415_DP ]
  REAL(DP), PARAMETER :: xgk(31) = [ 0.999484410050490637571325895705811_DP, &
    0.996893484074649540271630050918695_DP, 0.991630996870404594858628366109486_DP, &
    0.983668123279747209970032581605663_DP, 0.973116322501126268374693868423707_DP, &
    0.960021864968307512216871025581798_DP, 0.944374444748559979415831324037439_DP, &
    0.926200047429274325879324277080474_DP, 0.905573307699907798546522558925958_DP, &
    0.882560535792052681543116462530226_DP, 0.857205233546061098958658510658944_DP, &
    0.829565762382768397442898119732502_DP, 0.799727835821839083013668942322683_DP, &
    0.767777432104826194917977340974503_DP, 0.733790062453226804726171131369528_DP, &
    0.697850494793315796932292388026640_DP, 0.660061064126626961370053668149271_DP, &
    0.620526182989242861140477556431189_DP, 0.579345235826361691756024932172540_DP, &
    0.536624148142019899264169793311073_DP, 0.492480467861778574993693061207709_DP, &
    0.447033769538089176780609900322854_DP, 0.400401254830394392535476211542661_DP, &
    0.352704725530878113471037207089374_DP, 0.304073202273625077372677107199257_DP, &
    0.254636926167889846439805129817805_DP, 0.204525116682309891438957671002025_DP, &
    0.153869913608583546963794672743256_DP, 0.102806937966737030147096751318001_DP, &
    0.051471842555317695833025213166723_DP, 0.000000000000000000000000000000000_DP ]
  REAL(DP), PARAMETER :: wgk(31) = [ 0.001389013698677007624551591226760_DP, &
    0.003890461127099884051267201844516_DP, 0.006630703915931292173319826369750_DP, &
    0.009273279659517763428441146892024_DP, 0.011823015253496341742232898853251_DP, &
    0.014369729507045804812451432443580_DP, 0.016920889189053272627572289420322_DP, &
    0.019414141193942381173408951050128_DP, 0.021828035821609192297167485738339_DP, &
    0.024191162078080601365686370725232_DP, 0.026509954882333101610601709335075_DP, &
    0.028754048765041292843978785354334_DP, 0.030907257562387762472884252943092_DP, &
    0.032981447057483726031814191016854_DP, 0.034979338028060024137499670731468_DP, &
    0.036882364651821229223911065617136_DP, 0.038678945624727592950348651532281_DP, &
    0.040374538951535959111995279752468_DP, 0.041969810215164246147147541285970_DP, &
    0.043452539701356069316831728117073_DP, 0.044814800133162663192355551616723_DP, &
    0.046059238271006988116271735559374_DP, 0.047185546569299153945261478181099_DP, &
    0.048185861757087129140779492298305_DP, 0.049055434555029778887528165367238_DP, &
    0.049795683427074206357811569379942_DP, 0.050405921402782346840893085653585_DP, &
    0.050881795898749606492297473049805_DP, 0.051221547849258772170656282604944_DP, &
    0.051426128537459025933862879215781_DP, 0.051494729429451567558340433647099_DP ]
  !
  !           LIST OF MAJOR VARIABLES
  !           -----------------------
  !
  !           CENTR  - MID POINT OF THE INTERVAL
  !           HLGTH  - HALF-LENGTH OF THE INTERVAL
  !           DABSC  - ABSCISSA
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
  !* FIRST EXECUTABLE STATEMENT  DQK61
  epmach = D1MACH(4)
  uflow = D1MACH(1)
  !
  centr = 0.5_DP*(B+A)
  hlgth = 0.5_DP*(B-A)
  dhlgth = ABS(hlgth)
  !
  !           COMPUTE THE 61-POINT KRONROD APPROXIMATION TO THE
  !           INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
  !
  resg = 0._DP
  fc = F(centr)
  resk = wgk(31)*fc
  Resabs = ABS(resk)
  DO j = 1, 15
    jtw = j*2
    dabsc = hlgth*xgk(jtw)
    fval1 = F(centr-dabsc)
    fval2 = F(centr+dabsc)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1 + fval2
    resg = resg + wg(j)*fsum
    resk = resk + wgk(jtw)*fsum
    Resabs = Resabs + wgk(jtw)*(ABS(fval1)+ABS(fval2))
  END DO
  DO j = 1, 15
    jtwm1 = j*2 - 1
    dabsc = hlgth*xgk(jtwm1)
    fval1 = F(centr-dabsc)
    fval2 = F(centr+dabsc)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1 + fval2
    resk = resk + wgk(jtwm1)*fsum
    Resabs = Resabs + wgk(jtwm1)*(ABS(fval1)+ABS(fval2))
  END DO
  reskh = resk*0.5_DP
  Resasc = wgk(31)*ABS(fc-reskh)
  DO j = 1, 30
    Resasc = Resasc + wgk(j)*(ABS(fv1(j)-reskh)+ABS(fv2(j)-reskh))
  END DO
  Result = resk*hlgth
  Resabs = Resabs*dhlgth
  Resasc = Resasc*dhlgth
  Abserr = ABS((resk-resg)*hlgth)
  IF( Resasc/=0._DP .AND. Abserr/=0._DP )&
    Abserr = Resasc*MIN(1._DP,(0.2E+03_DP*Abserr/Resasc)**1.5_DP)
  IF( Resabs>uflow/(0.5E+02_DP*epmach) ) Abserr = MAX((epmach*0.5E+02_DP)*Resabs,&
    Abserr)
END SUBROUTINE DQK61
