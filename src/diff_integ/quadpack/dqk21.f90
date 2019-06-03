!** DQK21
SUBROUTINE DQK21(F,A,B,Result,Abserr,Resabs,Resasc)
  !>
  !  To compute I = Integral of F over (A,B), with error
  !                           estimate
  !                       J = Integral of ABS(F) over (A,B)
  !***
  ! **Library:**   SLATEC (QUADPACK)
  !***
  ! **Category:**  H2A1A2
  !***
  ! **Type:**      DOUBLE PRECISION (QK21-S, DQK21-D)
  !***
  ! **Keywords:**  21-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
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
  !                       Declared E X T E R N A L in the driver program.
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
  !                       RESULT is computed by applying the 21-POINT
  !                       KRONROD RULE (RESK) obtained by optimal addition
  !                       of abscissae to the 10-POINT GAUSS RULE (RESG).
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
  REAL(DP) :: absc, centr, dhlgth, epmach, fc, fsum, fval1, fval2, fv1(10), &
    fv2(10), hlgth, resg, resk, reskh, uflow
  !
  !           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
  !           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
  !           CORRESPONDING WEIGHTS ARE GIVEN.
  !
  !           XGK    - ABSCISSAE OF THE 21-POINT KRONROD RULE
  !                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 10-POINT
  !                    GAUSS RULE
  !                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
  !                    ADDED TO THE 10-POINT GAUSS RULE
  !
  !           WGK    - WEIGHTS OF THE 21-POINT KRONROD RULE
  !
  !           WG     - WEIGHTS OF THE 10-POINT GAUSS RULE
  !
  !
  ! GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS
  ! AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
  ! BELL LABS, NOV. 1981.
  !
  REAL(DP), PARAMETER :: wg(5) = [ 0.066671344308688137593568809893332D0, &
    0.149451349150580593145776339657697D0, 0.219086362515982043995534934228163D0, &
    0.269266719309996355091226921569469D0, 0.295524224714752870173892994651338D0 ]
  REAL(DP), PARAMETER :: xgk(11) = [ 0.995657163025808080735527280689003D0, &
    0.973906528517171720077964012084452D0, 0.930157491355708226001207180059508D0, &
    0.865063366688984510732096688423493D0, 0.780817726586416897063717578345042D0, &
    0.679409568299024406234327365114874D0, 0.562757134668604683339000099272694D0, &
    0.433395394129247190799265943165784D0, 0.294392862701460198131126603103866D0, &
    0.148874338981631210884826001129720D0, 0.000000000000000000000000000000000D0 ]
  REAL(DP), PARAMETER :: wgk(11) = [ 0.011694638867371874278064396062192D0, &
    0.032558162307964727478818972459390D0, 0.054755896574351996031381300244580D0, &
    0.075039674810919952767043140916190D0, 0.093125454583697605535065465083366D0, &
    0.109387158802297641899210590325805D0, 0.123491976262065851077958109831074D0, &
    0.134709217311473325928054001771707D0, 0.142775938577060080797094273138717D0, &
    0.147739104901338491374841515972068D0, 0.149445554002916905664936468389821D0 ]
  !
  !
  !           LIST OF MAJOR VARIABLES
  !           -----------------------
  !
  !           CENTR  - MID POINT OF THE INTERVAL
  !           HLGTH  - HALF-LENGTH OF THE INTERVAL
  !           ABSC   - ABSCISSA
  !           FVAL*  - FUNCTION VALUE
  !           RESG   - RESULT OF THE 10-POINT GAUSS FORMULA
  !           RESK   - RESULT OF THE 21-POINT KRONROD FORMULA
  !           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
  !                    I.E. TO I/(B-A)
  !
  !
  !           MACHINE DEPENDENT CONSTANTS
  !           ---------------------------
  !
  !           EPMACH IS THE LARGEST RELATIVE SPACING.
  !           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
  !
  !* FIRST EXECUTABLE STATEMENT  DQK21
  epmach = D1MACH(4)
  uflow = D1MACH(1)
  !
  centr = 0.5D+00*(A+B)
  hlgth = 0.5D+00*(B-A)
  dhlgth = ABS(hlgth)
  !
  !           COMPUTE THE 21-POINT KRONROD APPROXIMATION TO
  !           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
  !
  resg = 0.0D+00
  fc = F(centr)
  resk = wgk(11)*fc
  Resabs = ABS(resk)
  DO j = 1, 5
    jtw = 2*j
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
  DO j = 1, 5
    jtwm1 = 2*j - 1
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
  Resasc = wgk(11)*ABS(fc-reskh)
  DO j = 1, 10
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
END SUBROUTINE DQK21
