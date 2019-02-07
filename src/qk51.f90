!*==QK51.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK QK51
SUBROUTINE QK51(F,A,B,Result,Abserr,Resabs,Resasc)
  IMPLICIT NONE
  !*--QK515
  !***BEGIN PROLOGUE  QK51
  !***PURPOSE  To compute I = Integral of F over (A,B) with error
  !                           estimate
  !                       J = Integral of ABS(F) over (A,B)
  !***LIBRARY   SLATEC (QUADPACK)
  !***CATEGORY  H2A1A2
  !***TYPE      SINGLE PRECISION (QK51-S, DQK51-D)
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
  !           Real version
  !
  !           PARAMETERS
  !            ON ENTRY
  !              F      - Real
  !                       Function subroutine defining the integrand
  !                       function F(X). The actual name for F needs to be
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
  !                       RESULT is computed by applying the 51-point
  !                       Kronrod rule (RESK) obtained by optimal addition
  !                       of abscissae to the 25-point Gauss rule (RESG).
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
  !***END PROLOGUE  QK51
  !
  REAL A , absc , Abserr , B , centr , dhlgth , epmach , F , fc , fsum , &
    fval1 , fval2 , fv1 , fv2 , hlgth , Resabs , Resasc , resg , resk , &
    reskh , Result , R1MACH , uflow , wg , wgk , xgk
  INTEGER j , jtw , jtwm1
  EXTERNAL F
  !
  DIMENSION fv1(25) , fv2(25) , xgk(26) , wgk(26) , wg(13)
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
  SAVE xgk , wgk , wg
  DATA xgk(1) , xgk(2) , xgk(3) , xgk(4) , xgk(5) , xgk(6) , xgk(7) , xgk(8)&
    , xgk(9) , xgk(10) , xgk(11) , xgk(12) , xgk(13) , xgk(14)&
    /0.9992621049926098E+00 , 0.9955569697904981E+00 , &
    0.9880357945340772E+00 , 0.9766639214595175E+00 , &
    0.9616149864258425E+00 , 0.9429745712289743E+00 , &
    0.9207471152817016E+00 , 0.8949919978782754E+00 , &
    0.8658470652932756E+00 , 0.8334426287608340E+00 , &
    0.7978737979985001E+00 , 0.7592592630373576E+00 , &
    0.7177664068130844E+00 , 0.6735663684734684E+00/
  DATA xgk(15) , xgk(16) , xgk(17) , xgk(18) , xgk(19) , xgk(20) , xgk(21) , &
    xgk(22) , xgk(23) , xgk(24) , xgk(25) , xgk(26)&
    /0.6268100990103174E+00 , 0.5776629302412230E+00 , &
    0.5263252843347192E+00 , 0.4730027314457150E+00 , &
    0.4178853821930377E+00 , 0.3611723058093878E+00 , &
    0.3030895389311078E+00 , 0.2438668837209884E+00 , &
    0.1837189394210489E+00 , 0.1228646926107104E+00 , &
    0.6154448300568508E-01 , 0.0E+00/
  DATA wgk(1) , wgk(2) , wgk(3) , wgk(4) , wgk(5) , wgk(6) , wgk(7) , wgk(8)&
    , wgk(9) , wgk(10) , wgk(11) , wgk(12) , wgk(13) , wgk(14)&
    /0.1987383892330316E-02 , 0.5561932135356714E-02 , &
    0.9473973386174152E-02 , 0.1323622919557167E-01 , &
    0.1684781770912830E-01 , 0.2043537114588284E-01 , &
    0.2400994560695322E-01 , 0.2747531758785174E-01 , &
    0.3079230016738749E-01 , 0.3400213027432934E-01 , &
    0.3711627148341554E-01 , 0.4008382550403238E-01 , &
    0.4287284502017005E-01 , 0.4550291304992179E-01/
  DATA wgk(15) , wgk(16) , wgk(17) , wgk(18) , wgk(19) , wgk(20) , wgk(21) , &
    wgk(22) , wgk(23) , wgk(24) , wgk(25) , wgk(26)&
    /0.4798253713883671E-01 , 0.5027767908071567E-01 , &
    0.5236288580640748E-01 , 0.5425112988854549E-01 , &
    0.5595081122041232E-01 , 0.5743711636156783E-01 , &
    0.5868968002239421E-01 , 0.5972034032417406E-01 , &
    0.6053945537604586E-01 , 0.6112850971705305E-01 , &
    0.6147118987142532E-01 , 0.6158081806783294E-01/
  DATA wg(1) , wg(2) , wg(3) , wg(4) , wg(5) , wg(6) , wg(7) , wg(8) , &
    wg(9) , wg(10) , wg(11) , wg(12) , wg(13)/0.1139379850102629E-01 , &
    0.2635498661503214E-01 , 0.4093915670130631E-01 , &
    0.5490469597583519E-01 , 0.6803833381235692E-01 , &
    0.8014070033500102E-01 , 0.9102826198296365E-01 , &
    0.1005359490670506E+00 , 0.1085196244742637E+00 , &
    0.1148582591457116E+00 , 0.1194557635357848E+00 , &
    0.1222424429903100E+00 , 0.1231760537267155E+00/
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
  !***FIRST EXECUTABLE STATEMENT  QK51
  epmach = R1MACH(4)
  uflow = R1MACH(1)
  !
  centr = 0.5E+00*(A+B)
  hlgth = 0.5E+00*(B-A)
  dhlgth = ABS(hlgth)
  !
  !           COMPUTE THE 51-POINT KRONROD APPROXIMATION TO
  !           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
  !
  fc = F(centr)
  resg = wg(13)*fc
  resk = wgk(26)*fc
  Resabs = ABS(resk)
  DO j = 1 , 12
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
  DO j = 1 , 13
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
  Resasc = wgk(26)*ABS(fc-reskh)
  DO j = 1 , 25
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
END SUBROUTINE QK51
