!** DQK31
SUBROUTINE DQK31(F,A,B,Result,Abserr,Resabs,Resasc)
  !>
  !  To compute I = Integral of F over (A,B) with error
  !                           estimate
  !                       J = Integral of ABS(F) over (A,B)
  !***
  ! **Library:**   SLATEC (QUADPACK)
  !***
  ! **Category:**  H2A1A2
  !***
  ! **Type:**      DOUBLE PRECISION (QK31-S, DQK31-D)
  !***
  ! **Keywords:**  31-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
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
  !                       RESULT is computed by applying the 31-POINT
  !                       GAUSS-KRONROD RULE (RESK), obtained by optimal
  !                       addition of abscissae to the 15-POINT GAUSS
  !                       RULE (RESG).
  !
  !              ABSERR - Double precision
  !                       Estimate of the modulus of the modulus,
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
  USE service, ONLY : D1MACH
  INTERFACE
    REAL(DP) FUNCTION F(X)
      IMPORT DP
      REAL(DP) :: X
    END FUNCTION F
  END INTERFACE
  REAL(DP) :: A, Abserr, B, Resabs, Resasc, Result
  INTEGER :: j, jtw, jtwm1
  REAL(DP) :: absc, centr, dhlgth, epmach, fc, fsum, fval1, fval2, fv1(15), &
    fv2(15), hlgth, resg, resk, reskh, uflow
  !
  !           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
  !           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
  !           CORRESPONDING WEIGHTS ARE GIVEN.
  !
  !           XGK    - ABSCISSAE OF THE 31-POINT KRONROD RULE
  !                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 15-POINT
  !                    GAUSS RULE
  !                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
  !                    ADDED TO THE 15-POINT GAUSS RULE
  !
  !           WGK    - WEIGHTS OF THE 31-POINT KRONROD RULE
  !
  !           WG     - WEIGHTS OF THE 15-POINT GAUSS RULE
  !
  !
  ! GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS
  ! AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
  ! BELL LABS, NOV. 1981.
  !
  REAL(DP), PARAMETER :: wg(8) = [ 0.030753241996117268354628393577204D0, &
    0.070366047488108124709267416450667D0, 0.107159220467171935011869546685869D0, &
    0.139570677926154314447804794511028D0, 0.166269205816993933553200860481209D0, &
    0.186161000015562211026800561866423D0, 0.198431485327111576456118326443839D0, &
    0.202578241925561272880620199967519D0 ]
  REAL(DP), PARAMETER :: xgk(16) = [ 0.998002298693397060285172840152271D0, &
    0.987992518020485428489565718586613D0, 0.967739075679139134257347978784337D0, &
    0.937273392400705904307758947710209D0, 0.897264532344081900882509656454496D0, &
    0.848206583410427216200648320774217D0, 0.790418501442465932967649294817947D0, &
    0.724417731360170047416186054613938D0, 0.650996741297416970533735895313275D0, &
    0.570972172608538847537226737253911D0, 0.485081863640239680693655740232351D0, &
    0.394151347077563369897207370981045D0, 0.299180007153168812166780024266389D0, &
    0.201194093997434522300628303394596D0, 0.101142066918717499027074231447392D0, &
    0.000000000000000000000000000000000D0 ]
  REAL(DP), PARAMETER :: wgk(16) = [ 0.005377479872923348987792051430128D0, &
    0.015007947329316122538374763075807D0, 0.025460847326715320186874001019653D0, &
    0.035346360791375846222037948478360D0, 0.044589751324764876608227299373280D0, &
    0.053481524690928087265343147239430D0, 0.062009567800670640285139230960803D0, &
    0.069854121318728258709520077099147D0, 0.076849680757720378894432777482659D0, &
    0.083080502823133021038289247286104D0, 0.088564443056211770647275443693774D0, &
    0.093126598170825321225486872747346D0, 0.096642726983623678505179907627589D0, &
    0.099173598721791959332393173484603D0, 0.100769845523875595044946662617570D0, &
    0.101330007014791549017374792767493D0 ]
  !
  !
  !           LIST OF MAJOR VARIABLES
  !           -----------------------
  !           CENTR  - MID POINT OF THE INTERVAL
  !           HLGTH  - HALF-LENGTH OF THE INTERVAL
  !           ABSC   - ABSCISSA
  !           FVAL*  - FUNCTION VALUE
  !           RESG   - RESULT OF THE 15-POINT GAUSS FORMULA
  !           RESK   - RESULT OF THE 31-POINT KRONROD FORMULA
  !           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
  !                    I.E. TO I/(B-A)
  !
  !           MACHINE DEPENDENT CONSTANTS
  !           ---------------------------
  !           EPMACH IS THE LARGEST RELATIVE SPACING.
  !           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
  !* FIRST EXECUTABLE STATEMENT  DQK31
  epmach = D1MACH(4)
  uflow = D1MACH(1)
  !
  centr = 0.5D+00*(A+B)
  hlgth = 0.5D+00*(B-A)
  dhlgth = ABS(hlgth)
  !
  !           COMPUTE THE 31-POINT KRONROD APPROXIMATION TO
  !           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
  !
  fc = F(centr)
  resg = wg(8)*fc
  resk = wgk(16)*fc
  Resabs = ABS(resk)
  DO j = 1, 7
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
  DO j = 1, 8
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
  reskh = resk*0.5D+00
  Resasc = wgk(16)*ABS(fc-reskh)
  DO j = 1, 15
    Resasc = Resasc + wgk(j)*(ABS(fv1(j)-reskh)+ABS(fv2(j)-reskh))
  END DO
  Result = resk*hlgth
  Resabs = Resabs*dhlgth
  Resasc = Resasc*dhlgth
  Abserr = ABS((resk-resg)*hlgth)
  IF ( Resasc/=0.0D+00.AND.Abserr/=0.0D+00 )&
    Abserr = Resasc*MIN(0.1D+01,(0.2D+03*Abserr/Resasc)**1.5D+00)
  IF ( Resabs>uflow/(0.5D+02*epmach) ) Abserr = MAX((epmach*0.5D+02)*Resabs,&
    Abserr)
END SUBROUTINE DQK31
