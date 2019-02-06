!*==DQK61.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DQK61
      SUBROUTINE DQK61(F,A,B,Result,Abserr,Resabs,Resasc)
      IMPLICIT NONE
!*--DQK615
!***BEGIN PROLOGUE  DQK61
!***PURPOSE  To compute I = Integral of F over (A,B) with error
!                           estimate
!                       J = Integral of ABS(F) over (A,B)
!***LIBRARY   SLATEC (QUADPACK)
!***CATEGORY  H2A1A2
!***TYPE      DOUBLE PRECISION (QK61-S, DQK61-D)
!***KEYWORDS  61-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
!***AUTHOR  Piessens, Robert
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!           de Doncker, Elise
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!***DESCRIPTION
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
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH
!***REVISION HISTORY  (YYMMDD)
!   800101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  DQK61
!
      DOUBLE PRECISION A , dabsc , Abserr , B , centr , dhlgth , D1MACH , 
     &                 epmach , F , fc , fsum , fval1 , fval2 , fv1 , fv2 , 
     &                 hlgth , Resabs , Resasc , resg , resk , reskh , Result , 
     &                 uflow , wg , wgk , xgk
      INTEGER j , jtw , jtwm1
      EXTERNAL F
!
      DIMENSION fv1(30) , fv2(30) , xgk(31) , wgk(31) , wg(15)
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
      SAVE wg , xgk , wgk
      DATA wg(1)/0.007968192496166605615465883474674D0/
      DATA wg(2)/0.018466468311090959142302131912047D0/
      DATA wg(3)/0.028784707883323369349719179611292D0/
      DATA wg(4)/0.038799192569627049596801936446348D0/
      DATA wg(5)/0.048402672830594052902938140422808D0/
      DATA wg(6)/0.057493156217619066481721689402056D0/
      DATA wg(7)/0.065974229882180495128128515115962D0/
      DATA wg(8)/0.073755974737705206268243850022191D0/
      DATA wg(9)/0.080755895229420215354694938460530D0/
      DATA wg(10)/0.086899787201082979802387530715126D0/
      DATA wg(11)/0.092122522237786128717632707087619D0/
      DATA wg(12)/0.096368737174644259639468626351810D0/
      DATA wg(13)/0.099593420586795267062780282103569D0/
      DATA wg(14)/0.101762389748405504596428952168554D0/
      DATA wg(15)/0.102852652893558840341285636705415D0/
!
      DATA xgk(1)/0.999484410050490637571325895705811D0/
      DATA xgk(2)/0.996893484074649540271630050918695D0/
      DATA xgk(3)/0.991630996870404594858628366109486D0/
      DATA xgk(4)/0.983668123279747209970032581605663D0/
      DATA xgk(5)/0.973116322501126268374693868423707D0/
      DATA xgk(6)/0.960021864968307512216871025581798D0/
      DATA xgk(7)/0.944374444748559979415831324037439D0/
      DATA xgk(8)/0.926200047429274325879324277080474D0/
      DATA xgk(9)/0.905573307699907798546522558925958D0/
      DATA xgk(10)/0.882560535792052681543116462530226D0/
      DATA xgk(11)/0.857205233546061098958658510658944D0/
      DATA xgk(12)/0.829565762382768397442898119732502D0/
      DATA xgk(13)/0.799727835821839083013668942322683D0/
      DATA xgk(14)/0.767777432104826194917977340974503D0/
      DATA xgk(15)/0.733790062453226804726171131369528D0/
      DATA xgk(16)/0.697850494793315796932292388026640D0/
      DATA xgk(17)/0.660061064126626961370053668149271D0/
      DATA xgk(18)/0.620526182989242861140477556431189D0/
      DATA xgk(19)/0.579345235826361691756024932172540D0/
      DATA xgk(20)/0.536624148142019899264169793311073D0/
      DATA xgk(21)/0.492480467861778574993693061207709D0/
      DATA xgk(22)/0.447033769538089176780609900322854D0/
      DATA xgk(23)/0.400401254830394392535476211542661D0/
      DATA xgk(24)/0.352704725530878113471037207089374D0/
      DATA xgk(25)/0.304073202273625077372677107199257D0/
      DATA xgk(26)/0.254636926167889846439805129817805D0/
      DATA xgk(27)/0.204525116682309891438957671002025D0/
      DATA xgk(28)/0.153869913608583546963794672743256D0/
      DATA xgk(29)/0.102806937966737030147096751318001D0/
      DATA xgk(30)/0.051471842555317695833025213166723D0/
      DATA xgk(31)/0.000000000000000000000000000000000D0/
!
      DATA wgk(1)/0.001389013698677007624551591226760D0/
      DATA wgk(2)/0.003890461127099884051267201844516D0/
      DATA wgk(3)/0.006630703915931292173319826369750D0/
      DATA wgk(4)/0.009273279659517763428441146892024D0/
      DATA wgk(5)/0.011823015253496341742232898853251D0/
      DATA wgk(6)/0.014369729507045804812451432443580D0/
      DATA wgk(7)/0.016920889189053272627572289420322D0/
      DATA wgk(8)/0.019414141193942381173408951050128D0/
      DATA wgk(9)/0.021828035821609192297167485738339D0/
      DATA wgk(10)/0.024191162078080601365686370725232D0/
      DATA wgk(11)/0.026509954882333101610601709335075D0/
      DATA wgk(12)/0.028754048765041292843978785354334D0/
      DATA wgk(13)/0.030907257562387762472884252943092D0/
      DATA wgk(14)/0.032981447057483726031814191016854D0/
      DATA wgk(15)/0.034979338028060024137499670731468D0/
      DATA wgk(16)/0.036882364651821229223911065617136D0/
      DATA wgk(17)/0.038678945624727592950348651532281D0/
      DATA wgk(18)/0.040374538951535959111995279752468D0/
      DATA wgk(19)/0.041969810215164246147147541285970D0/
      DATA wgk(20)/0.043452539701356069316831728117073D0/
      DATA wgk(21)/0.044814800133162663192355551616723D0/
      DATA wgk(22)/0.046059238271006988116271735559374D0/
      DATA wgk(23)/0.047185546569299153945261478181099D0/
      DATA wgk(24)/0.048185861757087129140779492298305D0/
      DATA wgk(25)/0.049055434555029778887528165367238D0/
      DATA wgk(26)/0.049795683427074206357811569379942D0/
      DATA wgk(27)/0.050405921402782346840893085653585D0/
      DATA wgk(28)/0.050881795898749606492297473049805D0/
      DATA wgk(29)/0.051221547849258772170656282604944D0/
      DATA wgk(30)/0.051426128537459025933862879215781D0/
      DATA wgk(31)/0.051494729429451567558340433647099D0/
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
!***FIRST EXECUTABLE STATEMENT  DQK61
      epmach = D1MACH(4)
      uflow = D1MACH(1)
!
      centr = 0.5D+00*(B+A)
      hlgth = 0.5D+00*(B-A)
      dhlgth = ABS(hlgth)
!
!           COMPUTE THE 61-POINT KRONROD APPROXIMATION TO THE
!           INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
!
      resg = 0.0D+00
      fc = F(centr)
      resk = wgk(31)*fc
      Resabs = ABS(resk)
      DO j = 1 , 15
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
      ENDDO
      DO j = 1 , 15
        jtwm1 = j*2 - 1
        dabsc = hlgth*xgk(jtwm1)
        fval1 = F(centr-dabsc)
        fval2 = F(centr+dabsc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1 + fval2
        resk = resk + wgk(jtwm1)*fsum
        Resabs = Resabs + wgk(jtwm1)*(ABS(fval1)+ABS(fval2))
      ENDDO
      reskh = resk*0.5D+00
      Resasc = wgk(31)*ABS(fc-reskh)
      DO j = 1 , 30
        Resasc = Resasc + wgk(j)*(ABS(fv1(j)-reskh)+ABS(fv2(j)-reskh))
      ENDDO
      Result = resk*hlgth
      Resabs = Resabs*dhlgth
      Resasc = Resasc*dhlgth
      Abserr = ABS((resk-resg)*hlgth)
      IF ( Resasc/=0.0D+00.AND.Abserr/=0.0D+00 )
     &     Abserr = Resasc*MIN(0.1D+01,(0.2D+03*Abserr/Resasc)**1.5D+00)
      IF ( Resabs>uflow/(0.5D+02*epmach) ) Abserr = MAX((epmach*0.5D+02)*Resabs,
     &     Abserr)
      END SUBROUTINE DQK61
