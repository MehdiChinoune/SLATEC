!** QK51
PURE SUBROUTINE QK51(F,A,B,Result,Abserr,Resabs,Resasc)
  !> To compute I = Integral of F over (A,B) with error estimate
  !  J = Integral of ABS(F) over (A,B)
  !***
  ! **Library:**   SLATEC (QUADPACK)
  !***
  ! **Category:**  H2A1A2
  !***
  ! **Type:**      SINGLE PRECISION (QK51-S, DQK51-D)
  !***
  ! **Keywords:**  51-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
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
    REAL(SP) PURE FUNCTION F(X)
      IMPORT SP
      REAL(SP), INTENT(IN) :: X
    END FUNCTION F
  END INTERFACE
  REAL(SP), INTENT(IN) :: A, B
  REAL(SP), INTENT(OUT) :: Abserr, Resabs, Resasc, Result
  !
  INTEGER :: j, jtw, jtwm1
  REAL(SP) :: absc, centr, dhlgth, epmach, fc, fsum, fval1, fval2, fv1(25), &
    fv2(25), hlgth, resg, resk, reskh, uflow
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
  REAL(SP), PARAMETER :: xgk(26) = [ 0.9992621049926098_SP, 0.9955569697904981_SP, &
    0.9880357945340772_SP, 0.9766639214595175_SP, 0.9616149864258425_SP, &
    0.9429745712289743_SP, 0.9207471152817016_SP, 0.8949919978782754_SP, &
    0.8658470652932756_SP, 0.8334426287608340_SP, 0.7978737979985001_SP, &
    0.7592592630373576_SP, 0.7177664068130844_SP, 0.6735663684734684_SP, &
    0.6268100990103174_SP, 0.5776629302412230_SP, 0.5263252843347192_SP, &
    0.4730027314457150_SP, 0.4178853821930377_SP, 0.3611723058093878_SP, &
    0.3030895389311078_SP, 0.2438668837209884_SP, 0.1837189394210489_SP, &
    0.1228646926107104_SP, 0.6154448300568508E-01_SP, 0._SP ]
  REAL(SP), PARAMETER :: wgk(26) = [ 0.1987383892330316E-02_SP, 0.5561932135356714E-02_SP, &
    0.9473973386174152E-02_SP, 0.1323622919557167E-01_SP, 0.1684781770912830E-01_SP, &
    0.2043537114588284E-01_SP, 0.2400994560695322E-01_SP, 0.2747531758785174E-01_SP, &
    0.3079230016738749E-01_SP, 0.3400213027432934E-01_SP, 0.3711627148341554E-01_SP, &
    0.4008382550403238E-01_SP, 0.4287284502017005E-01_SP, 0.4550291304992179E-01_SP, &
    0.4798253713883671E-01_SP, 0.5027767908071567E-01_SP, 0.5236288580640748E-01_SP, &
    0.5425112988854549E-01_SP, 0.5595081122041232E-01_SP, 0.5743711636156783E-01_SP, &
    0.5868968002239421E-01_SP, 0.5972034032417406E-01_SP, 0.6053945537604586E-01_SP, &
    0.6112850971705305E-01_SP, 0.6147118987142532E-01_SP, 0.6158081806783294E-01_SP ]
  REAL(SP), PARAMETER :: wg(13) = [ 0.1139379850102629E-01_SP, 0.2635498661503214E-01_SP, &
    0.4093915670130631E-01_SP, 0.5490469597583519E-01_SP, 0.6803833381235692E-01_SP, &
    0.8014070033500102E-01_SP, 0.9102826198296365E-01_SP, 0.1005359490670506E+00_SP, &
    0.1085196244742637E+00_SP, 0.1148582591457116E+00_SP, 0.1194557635357848E+00_SP, &
    0.1222424429903100E+00_SP, 0.1231760537267155E+00_SP ]
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
  !* FIRST EXECUTABLE STATEMENT  QK51
  epmach = R1MACH(4)
  uflow = R1MACH(1)
  !
  centr = 0.5_SP*(A+B)
  hlgth = 0.5_SP*(B-A)
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
  END DO
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
  END DO
  reskh = resk*0.5_SP
  Resasc = wgk(26)*ABS(fc-reskh)
  DO j = 1, 25
    Resasc = Resasc + wgk(j)*(ABS(fv1(j)-reskh)+ABS(fv2(j)-reskh))
  END DO
  Result = resk*hlgth
  Resabs = Resabs*dhlgth
  Resasc = Resasc*dhlgth
  Abserr = ABS((resk-resg)*hlgth)
  IF( Resasc/=0._SP .AND. Abserr/=0._SP )&
    Abserr = Resasc*MIN(1._SP,(200._SP*Abserr/Resasc)**1.5_SP)
  IF( Resabs>uflow/(50._SP*epmach) ) Abserr = MAX((epmach*50._SP)*Resabs,Abserr)
  !
END SUBROUTINE QK51