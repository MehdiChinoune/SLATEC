!** DQC25C
SUBROUTINE DQC25C(F,A,B,C,Result,Abserr,Krul,Neval)
  !>
  !  To compute I = Integral of F*W over (A,B) with
  !            error estimate, where W(X) = 1/(X-C)
  !***
  ! **Library:**   SLATEC (QUADPACK)
  !***
  ! **Category:**  H2A2A2, J4
  !***
  ! **Type:**      DOUBLE PRECISION (QC25C-S, DQC25C-D)
  !***
  ! **Keywords:**  25-POINT CLENSHAW-CURTIS INTEGRATION, QUADPACK, QUADRATURE
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
  !        Integration rules for the computation of CAUCHY
  !        PRINCIPAL VALUE integrals
  !        Standard fortran subroutine
  !        Double precision version
  !
  !        PARAMETERS
  !           F      - Double precision
  !                    Function subprogram defining the integrand function
  !                    F(X). The actual name for F needs to be declared
  !                    E X T E R N A L  in the driver program.
  !
  !           A      - Double precision
  !                    Left end point of the integration interval
  !
  !           B      - Double precision
  !                    Right end point of the integration interval, B.GT.A
  !
  !           C      - Double precision
  !                    Parameter in the WEIGHT function
  !
  !           RESULT - Double precision
  !                    Approximation to the integral
  !                    result is computed by using a generalized
  !                    Clenshaw-Curtis method if C lies within ten percent
  !                    of the integration interval. In the other case the
  !                    15-point Kronrod rule obtained by optimal addition
  !                    of abscissae to the 7-point Gauss rule, is applied.
  !
  !           ABSERR - Double precision
  !                    Estimate of the modulus of the absolute error,
  !                    which should equal or exceed ABS(I-RESULT)
  !
  !           KRUL   - Integer
  !                    Key which is decreased by 1 if the 15-point
  !                    Gauss-Kronrod scheme has been used
  !
  !           NEVAL  - Integer
  !                    Number of integrand evaluations
  !
  ! ......................................................................
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  DQCHEB, DQK15W, DQWGTC

  !* REVISION HISTORY  (YYMMDD)
  !   810101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)

  !
  INTERFACE
    REAL(8) FUNCTION F(X)
      REAL(8) :: X
    END FUNCTION F
  END INTERFACE
  INTEGER :: Krul, Neval
  REAL(8) :: A, Abserr, B, C, Result
  INTEGER :: i, isym, k, kp
  REAL(8) :: ak22, amom0, amom1, amom2, cc, centr, cheb12(13), cheb24(25), &
    fval(25), hlgth, p2, p3, p4, resabs, resasc, res12, res24, u
  !
  !           THE VECTOR X CONTAINS THE VALUES COS(K*PI/24),
  !           K = 1, ..., 11, TO BE USED FOR THE CHEBYSHEV SERIES
  !           EXPANSION OF F
  !
  REAL(8), PARAMETER :: x(11) = [ 0.9914448613738104D+00, 0.9659258262890683D+00, &
    0.9238795325112868D+00, 0.8660254037844386D+00, 0.7933533402912352D+00, &
    0.7071067811865475D+00, 0.6087614290087206D+00, 0.5000000000000000D+00, &
    0.3826834323650898D+00, 0.2588190451025208D+00, 0.1305261922200516D+00 ]
  !
  !           LIST OF MAJOR VARIABLES
  !           ----------------------
  !           FVAL   - VALUE OF THE FUNCTION F AT THE POINTS
  !                    COS(K*PI/24),  K = 0, ..., 24
  !           CHEB12 - CHEBYSHEV SERIES EXPANSION COEFFICIENTS,
  !                    FOR THE FUNCTION F, OF DEGREE 12
  !           CHEB24 - CHEBYSHEV SERIES EXPANSION COEFFICIENTS,
  !                    FOR THE FUNCTION F, OF DEGREE 24
  !           RES12  - APPROXIMATION TO THE INTEGRAL CORRESPONDING
  !                    TO THE USE OF CHEB12
  !           RES24  - APPROXIMATION TO THE INTEGRAL CORRESPONDING
  !                    TO THE USE OF CHEB24
  !           DQWGTC - EXTERNAL FUNCTION SUBPROGRAM DEFINING
  !                    THE WEIGHT FUNCTION
  !           HLGTH  - HALF-LENGTH OF THE INTERVAL
  !           CENTR  - MID POINT OF THE INTERVAL
  !
  !
  !           CHECK THE POSITION OF C.
  !
  !* FIRST EXECUTABLE STATEMENT  DQC25C
  cc = (0.2D+01*C-B-A)/(B-A)
  IF ( ABS(cc)<0.11D+01 ) THEN
    !
    !           USE THE GENERALIZED CLENSHAW-CURTIS METHOD.
    !
    hlgth = 0.5D+00*(B-A)
    centr = 0.5D+00*(B+A)
    Neval = 25
    fval(1) = 0.5D+00*F(hlgth+centr)
    fval(13) = F(centr)
    fval(25) = 0.5D+00*F(centr-hlgth)
    DO i = 2, 12
      u = hlgth*x(i-1)
      isym = 26 - i
      fval(i) = F(u+centr)
      fval(isym) = F(centr-u)
    END DO
    !
    !           COMPUTE THE CHEBYSHEV SERIES EXPANSION.
    !
    CALL DQCHEB(x,fval,cheb12,cheb24)
    !
    !           THE MODIFIED CHEBYSHEV MOMENTS ARE COMPUTED BY FORWARD
    !           RECURSION, USING AMOM0 AND AMOM1 AS STARTING VALUES.
    !
    amom0 = LOG(ABS((0.1D+01-cc)/(0.1D+01+cc)))
    amom1 = 0.2D+01 + cc*amom0
    res12 = cheb12(1)*amom0 + cheb12(2)*amom1
    res24 = cheb24(1)*amom0 + cheb24(2)*amom1
    DO k = 3, 13
      amom2 = 0.2D+01*cc*amom1 - amom0
      ak22 = (k-2)*(k-2)
      IF ( (k/2)*2==k ) amom2 = amom2 - 0.4D+01/(ak22-0.1D+01)
      res12 = res12 + cheb12(k)*amom2
      res24 = res24 + cheb24(k)*amom2
      amom0 = amom1
      amom1 = amom2
    END DO
    DO k = 14, 25
      amom2 = 0.2D+01*cc*amom1 - amom0
      ak22 = (k-2)*(k-2)
      IF ( (k/2)*2==k ) amom2 = amom2 - 0.4D+01/(ak22-0.1D+01)
      res24 = res24 + cheb24(k)*amom2
      amom0 = amom1
      amom1 = amom2
    END DO
    Result = res24
    Abserr = ABS(res24-res12)
  ELSE
    !
    !           APPLY THE 15-POINT GAUSS-KRONROD SCHEME.
    !
    Krul = Krul - 1
    CALL DQK15W(F,DQWGTC,C,p2,p3,p4,kp,A,B,Result,Abserr,resabs,resasc)
    Neval = 15
    IF ( resasc==Abserr ) Krul = Krul + 1
  END IF
END SUBROUTINE DQC25C
