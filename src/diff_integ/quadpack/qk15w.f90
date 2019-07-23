!** QK15W
PURE SUBROUTINE QK15W(F,W,P1,P2,P3,P4,Kp,A,B,Result,Abserr,Resabs,Resasc)
  !> To compute I = Integral of F*W over (A,B), with error estimate
  !  J = Integral of ABS(F*W) over (A,B)
  !***
  ! **Library:**   SLATEC (QUADPACK)
  !***
  ! **Category:**  H2A2A2
  !***
  ! **Type:**      SINGLE PRECISION (QK15W-S, DQK15W-D)
  !***
  ! **Keywords:**  15-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
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
  !             ON ENTRY
  !              F      - Real
  !                       Function subprogram defining the integrand
  !                       function F(X). The actual name for F needs to be
  !                       declared E X T E R N A L in the driver program.
  !
  !              W      - Real
  !                       Function subprogram defining the integrand
  !                       WEIGHT function W(X). The actual name for W
  !                       needs to be declared E X T E R N A L in the
  !                       calling program.
  !
  !              P1, P2, P3, P4 - Real
  !                       Parameters in the WEIGHT function
  !
  !              KP     - Integer
  !                       Key for indicating the type of WEIGHT function
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
  !                       RESULT is computed by applying the 15-point
  !                       Kronrod rule (RESK) obtained by optimal addition
  !                       of abscissae to the 7-point Gauss rule (RESG).
  !
  !              ABSERR - Real
  !                       Estimate of the modulus of the absolute error,
  !                       which should equal or exceed ABS(I-RESULT)
  !
  !              RESABS - Real
  !                       Approximation to the integral of ABS(F)
  !
  !              RESASC - Real
  !                       Approximation to the integral of ABS(F-I/(B-A))
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  R1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   810101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : R1MACH
  INTERFACE
    REAL(SP) PURE FUNCTION F(X)
      IMPORT SP
      REAL(SP), INTENT(IN) :: X
    END FUNCTION F
    REAL(SP) PURE FUNCTION W(X,P1,P2,P3,P4,Kp)
      IMPORT SP
      REAL(SP), INTENT(IN) :: X, P1, P2, P3, P4
      INTEGER, INTENT(IN) :: Kp
    END FUNCTION
  END INTERFACE
  INTEGER, INTENT(IN) :: Kp
  REAL(SP), INTENT(IN) :: A, B, P1, P2, P3, P4
  REAL(SP), INTENT(OUT) :: Abserr, Resabs, Resasc, Result
  !
  INTEGER :: j, jtw, jtwm1
  REAL(SP) :: absc, absc1, absc2, centr, dhlgth, epmach, fc, fsum, fval1, &
    fval2, fv1(7), fv2(7), hlgth, resg, resk, reskh, uflow
  !
  !           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
  !           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
  !           CORRESPONDING WEIGHTS ARE GIVEN.
  !
  !           XGK    - ABSCISSAE OF THE 15-POINT GAUSS-KRONROD RULE
  !                    XGK(2), XGK(4), ... ABSCISSAE OF THE 7-POINT
  !                    GAUSS RULE
  !                    XGK(1), XGK(3), ... ABSCISSAE WHICH ARE OPTIMALLY
  !                    ADDED TO THE 7-POINT GAUSS RULE
  !
  !           WGK    - WEIGHTS OF THE 15-POINT GAUSS-KRONROD RULE
  !
  !           WG     - WEIGHTS OF THE 7-POINT GAUSS RULE
  !
  REAL(SP), PARAMETER :: xgk(8) = [ 0.9914553711208126_SP, 0.9491079123427585_SP, &
    0.8648644233597691_SP, 0.7415311855993944_SP, 0.5860872354676911_SP, &
    0.4058451513773972_SP, 0.2077849550078985_SP, 0.0000000000000000_SP ]
  REAL(SP), PARAMETER :: wgk(8) = [ 0.2293532201052922E-01_SP, 0.6309209262997855E-01_SP, &
    0.1047900103222502_SP, 0.1406532597155259_SP, 0.1690047266392679_SP, &
    0.1903505780647854_SP, 0.2044329400752989_SP, 0.2094821410847278_SP ]
  REAL(SP), PARAMETER :: wg(4) = [ 0.1294849661688697_SP, 0.2797053914892767_SP, &
    0.3818300505051889_SP, 0.4179591836734694_SP ]
  !
  !           LIST OF MAJOR VARIABLES
  !           -----------------------
  !
  !           CENTR  - MID POINT OF THE INTERVAL
  !           HLGTH  - HALF-LENGTH OF THE INTERVAL
  !           ABSC*  - ABSCISSA
  !           FVAL*  - FUNCTION VALUE
  !           RESG   - RESULT OF THE 7-POINT GAUSS FORMULA
  !           RESK   - RESULT OF THE 15-POINT KRONROD FORMULA
  !           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F*W OVER (A,B),
  !                    I.E. TO I/(B-A)
  !
  !           MACHINE DEPENDENT CONSTANTS
  !           ---------------------------
  !
  !           EPMACH IS THE LARGEST RELATIVE SPACING.
  !           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
  !
  !* FIRST EXECUTABLE STATEMENT  QK15W
  epmach = R1MACH(4)
  uflow = R1MACH(1)
  !
  centr = 0.5_SP*(A+B)
  hlgth = 0.5_SP*(B-A)
  dhlgth = ABS(hlgth)
  !
  !           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO THE
  !           INTEGRAL, AND ESTIMATE THE ERROR.
  !
  fc = F(centr)*W(centr,P1,P2,P3,P4,Kp)
  resg = wg(4)*fc
  resk = wgk(8)*fc
  Resabs = ABS(resk)
  DO j = 1, 3
    jtw = j*2
    absc = hlgth*xgk(jtw)
    absc1 = centr - absc
    absc2 = centr + absc
    fval1 = F(absc1)*W(absc1,P1,P2,P3,P4,Kp)
    fval2 = F(absc2)*W(absc2,P1,P2,P3,P4,Kp)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1 + fval2
    resg = resg + wg(j)*fsum
    resk = resk + wgk(jtw)*fsum
    Resabs = Resabs + wgk(jtw)*(ABS(fval1)+ABS(fval2))
  END DO
  DO j = 1, 4
    jtwm1 = j*2 - 1
    absc = hlgth*xgk(jtwm1)
    absc1 = centr - absc
    absc2 = centr + absc
    fval1 = F(absc1)*W(absc1,P1,P2,P3,P4,Kp)
    fval2 = F(absc2)*W(absc2,P1,P2,P3,P4,Kp)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1 + fval2
    resk = resk + wgk(jtwm1)*fsum
    Resabs = Resabs + wgk(jtwm1)*(ABS(fval1)+ABS(fval2))
  END DO
  reskh = resk*0.5_SP
  Resasc = wgk(8)*ABS(fc-reskh)
  DO j = 1, 7
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
END SUBROUTINE QK15W