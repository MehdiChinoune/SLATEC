!** QK21
PURE SUBROUTINE QK21(F,A,B,Result,Abserr,Resabs,Resasc)
  !> To compute I = Integral of F over (A,B) with error estimate
  !  J = Integral of ABS(F) over (A,B)
  !***
  ! **Library:**   SLATEC (QUADPACK)
  !***
  ! **Category:**  H2A1A2
  !***
  ! **Type:**      SINGLE PRECISION (QK21-S, DQK21-D)
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
  !           Real version
  !
  !           PARAMETERS
  !            ON ENTRY
  !              F      - Real
  !                       Function subprogram defining the integrand
  !                       FUNCTION F(X). The actual name for F needs to be
  !                       Declared E X T E R N A L in the driver program.
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
  !                       RESULT is computed by applying the 21-POINT
  !                       KRONROD RULE (RESK) obtained by optimal addition
  !                       of abscissae to the 10-POINT GAUSS RULE (RESG).
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
  REAL(SP) :: absc, centr, dhlgth, epmach, fc, fsum, fval1, fval2, fv1(10), &
    fv2(10), hlgth, resg, resk, reskh, uflow
  INTEGER :: j, jtw, jtwm1
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
  REAL(SP), PARAMETER :: xgk(11) = [ 0.9956571630258081_SP, 0.9739065285171717_SP, &
    0.9301574913557082_SP, 0.8650633666889845_SP, 0.7808177265864169_SP, &
    0.6794095682990244_SP, 0.5627571346686047_SP, 0.4333953941292472_SP, &
    0.2943928627014602_SP, 0.1488743389816312_SP, 0.0000000000000000_SP ]
  REAL(SP), PARAMETER :: wgk(11) = [ 0.1169463886737187E-01_SP, 0.3255816230796473E-01_SP, &
    0.5475589657435200E-01_SP, 0.7503967481091995E-01_SP, 0.9312545458369761E-01_SP, &
    0.1093871588022976_SP, 0.1234919762620659_SP, 0.1347092173114733_SP, &
    0.1427759385770601_SP, 0.1477391049013385_SP, 0.1494455540029169_SP ]
  REAL(SP), PARAMETER :: wg(5) = [ 0.6667134430868814E-01_SP, 0.1494513491505806_SP, &
    0.2190863625159820_SP, 0.2692667193099964_SP, 0.2955242247147529_SP ]
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
  !* FIRST EXECUTABLE STATEMENT  QK21
  epmach = R1MACH(4)
  uflow = R1MACH(1)
  !
  centr = 0.5_SP*(A+B)
  hlgth = 0.5_SP*(B-A)
  dhlgth = ABS(hlgth)
  !
  !           COMPUTE THE 21-POINT KRONROD APPROXIMATION TO
  !           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
  !
  resg = 0._SP
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
  reskh = resk*0.5_SP
  Resasc = wgk(11)*ABS(fc-reskh)
  DO j = 1, 10
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
END SUBROUTINE QK21