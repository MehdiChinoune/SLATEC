!** QC25F
PURE SUBROUTINE QC25F(F,A,B,Omega,Integr,Nrmom,Maxp1,Ksave,Result,Abserr,Neval,&
    Resabs,Resasc,Momcom,Chebmo)
  !> To compute the integral I=Integral of F(X) over (A,B) Where W(X) = COS(OMEGA*X)
  !  or W(X)=SIN(OMEGA*X) and to compute J = Integral of ABS(F) over (A,B).
  !  For small value of OMEGA or small intervals (A,B) the 15-point GAUSS-KRONRO
  !  Rule is used. Otherwise a generalized CLENSHAW-CURTIS method is used.
  !***
  ! **Library:**   SLATEC (QUADPACK)
  !***
  ! **Category:**  H2A2A2
  !***
  ! **Type:**      SINGLE PRECISION (QC25F-S, DQC25F-D)
  !***
  ! **Keywords:**  CLENSHAW-CURTIS METHOD, GAUSS-KRONROD RULES,
  !             INTEGRATION RULES FOR FUNCTIONS WITH COS OR SIN FACTOR,
  !             QUADPACK, QUADRATURE
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
  !        Integration rules for functions with COS or SIN factor
  !        Standard fortran subroutine
  !        Real version
  !
  !        PARAMETERS
  !         ON ENTRY
  !           F      - Real
  !                    Function subprogram defining the integrand
  !                    function F(X). The actual name for F needs to
  !                    be declared E X T E R N A L in the calling program.
  !
  !           A      - Real
  !                    Lower limit of integration
  !
  !           B      - Real
  !                    Upper limit of integration
  !
  !           OMEGA  - Real
  !                    Parameter in the WEIGHT function
  !
  !           INTEGR - Integer
  !                    Indicates which WEIGHT function is to be used
  !                       INTEGR = 1   W(X) = COS(OMEGA*X)
  !                       INTEGR = 2   W(X) = SIN(OMEGA*X)
  !
  !           NRMOM  - Integer
  !                    The length of interval (A,B) is equal to the length
  !                    of the original integration interval divided by
  !                    2**NRMOM (we suppose that the routine is used in an
  !                    adaptive integration process, otherwise set
  !                    NRMOM = 0). NRMOM must be zero at the first call.
  !
  !           MAXP1  - Integer
  !                    Gives an upper bound on the number of Chebyshev
  !                    moments which can be stored, i.e. for the
  !                    intervals of lengths ABS(BB-AA)*2**(-L),
  !                    L = 0,1,2, ..., MAXP1-2.
  !
  !           KSAVE  - Integer
  !                    Key which is one when the moments for the
  !                    current interval have been computed
  !
  !         ON RETURN
  !           RESULT - Real
  !                    Approximation to the integral I
  !
  !           ABSERR - Real
  !                    Estimate of the modulus of the absolute
  !                    error, which should equal or exceed ABS(I-RESULT)
  !
  !           NEVAL  - Integer
  !                    Number of integrand evaluations
  !
  !           RESABS - Real
  !                    Approximation to the integral J
  !
  !           RESASC - Real
  !                    Approximation to the integral of ABS(F-I/(B-A))
  !
  !         ON ENTRY AND RETURN
  !           MOMCOM - Integer
  !                    For each interval length we need to compute the
  !                    Chebyshev moments. MOMCOM counts the number of
  !                    intervals for which these moments have already been
  !                    computed. If NRMOM<MOMCOM or KSAVE = 1, the
  !                    Chebyshev moments for the interval (A,B) have
  !                    already been computed and stored, otherwise we
  !                    compute them and we increase MOMCOM.
  !
  !           CHEBMO - Real
  !                    Array of dimension at least (MAXP1,25) containing
  !                    the modified Chebyshev moments for the first MOMCOM
  !                    MOMCOM interval lengths
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  QCHEB, QK15W, QWGTF, R1MACH, SGTSL

  !* REVISION HISTORY  (YYMMDD)
  !   810101  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : R1MACH
  USE lapack, ONLY : SGTSV
  !
  INTERFACE
    REAL(SP) PURE FUNCTION F(X)
      IMPORT SP
      REAL(SP), INTENT(IN) :: X
    END FUNCTION F
  END INTERFACE
  INTEGER, INTENT(IN) :: Integr, Ksave, Maxp1, Nrmom
  INTEGER, INTENT(INOUT) :: Momcom
  INTEGER, INTENT(OUT) :: Neval
  REAL(SP), INTENT(IN) :: A, B, Omega
  REAL(SP), INTENT(INOUT) :: Chebmo(Maxp1,25)
  REAL(SP), INTENT(OUT) :: Abserr, Resabs, Resasc, Result
  !
  INTEGER :: i, iers, isym, j, k, m, noequ, noeq1
  REAL(SP) :: ac, an, an2, as, asap, ass, centr, cheb12(13), cheb24(25), conc, &
    cons, cospar, d(25), d1(25), d2(25), estc, ests, fval(25), hlgth, oflow, &
    parint, par2, par22, p2, p3, p4, resc12, resc24, ress12, ress24, sinpar
  REAL(SP), TARGET :: v(28)
  REAL(SP), POINTER :: v2(:,:)
  !
  !           THE VECTOR X CONTAINS THE VALUES COS(K*PI/24)
  !           K = 1, ...,11, TO BE USED FOR THE CHEBYSHEV EXPANSION OF F
  !
  REAL(SP), PARAMETER :: x(11) = [ 0.9914448613738104_SP, 0.9659258262890683_SP, &
    0.9238795325112868_SP, 0.8660254037844386_SP, 0.7933533402912352_SP, &
    0.7071067811865475_SP, 0.6087614290087206_SP, 0.5000000000000000_SP, &
    0.3826834323650898_SP, 0.2588190451025208_SP, 0.1305261922200516_SP ]
  !
  !           LIST OF MAJOR VARIABLES
  !           -----------------------
  !
  !           CENTR  - MID POINT OF THE INTEGRATION INTERVAL
  !           HLGTH  - HALF-LENGTH OF THE INTEGRATION INTERVAL
  !           FVAL   - VALUE OF THE FUNCTION F AT THE POINTS
  !                    (B-A)*0.5*COS(K*PI/12) + (B+A)*0.5,
  !                    K = 0, ..., 24
  !           CHEB12 - COEFFICIENTS OF THE CHEBYSHEV SERIES EXPANSION
  !                    OF DEGREE 12, FOR THE FUNCTION F, IN THE
  !                    INTERVAL (A,B)
  !           CHEB24 - COEFFICIENTS OF THE CHEBYSHEV SERIES EXPANSION
  !                    OF DEGREE 24, FOR THE FUNCTION F, IN THE
  !                    INTERVAL (A,B)
  !           RESC12 - APPROXIMATION TO THE INTEGRAL OF
  !                    COS(0.5*(B-A)*OMEGA*X)*F(0.5*(B-A)*X+0.5*(B+A))
  !                    OVER (-1,+1), USING THE CHEBYSHEV SERIES
  !                    EXPANSION OF DEGREE 12
  !           RESC24 - APPROXIMATION TO THE SAME INTEGRAL, USING THE
  !                    CHEBYSHEV SERIES EXPANSION OF DEGREE 24
  !           RESS12 - THE ANALOGUE OF RESC12 FOR THE SINE
  !           RESS24 - THE ANALOGUE OF RESC24 FOR THE SINE
  !
  !
  !           MACHINE DEPENDENT CONSTANT
  !           --------------------------
  !
  !           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
  !
  !* FIRST EXECUTABLE STATEMENT  QC25F
  oflow = R1MACH(2)
  !
  centr = 0.5_SP*(B+A)
  hlgth = 0.5_SP*(B-A)
  parint = Omega*hlgth
  !
  !           COMPUTE THE INTEGRAL USING THE 15-POINT GAUSS-KRONROD
  !           FORMULA IF THE VALUE OF THE PARAMETER IN THE INTEGRAND
  !           IS SMALL.
  !
  IF( ABS(parint)>2._SP ) THEN
    !
    !           COMPUTE THE INTEGRAL USING THE GENERALIZED CLENSHAW-
    !           CURTIS METHOD.
    !
    conc = hlgth*COS(centr*Omega)
    cons = hlgth*SIN(centr*Omega)
    Resasc = oflow
    Neval = 25
    !
    !           CHECK WHETHER THE CHEBYSHEV MOMENTS FOR THIS INTERVAL
    !           HAVE ALREADY BEEN COMPUTED.
    !
    IF( Nrmom>=Momcom .AND. Ksave/=1 ) THEN
      !
      !           COMPUTE A NEW SET OF CHEBYSHEV MOMENTS.
      !
      m = Momcom + 1
      par2 = parint*parint
      par22 = par2 + 2._SP
      sinpar = SIN(parint)
      cospar = COS(parint)
      !
      !           COMPUTE THE CHEBYSHEV MOMENTS WITH RESPECT TO COSINE.
      !
      v(1) = 2._SP*sinpar/parint
      v(2) = (8._SP*cospar+(par2+par2-8._SP)*sinpar/parint)/par2
      v(3) = (0.32E+02_SP*(par2-0.12E+02_SP)*cospar+(2._SP*((par2-0.80E+02_SP)*par2&
        +0.192E+03_SP)*sinpar)/parint)/(par2*par2)
      ac = 8._SP*cospar
      as = 0.24E+02_SP*parint*sinpar
      IF( ABS(parint)>0.24E+02 ) THEN
        !
        !           COMPUTE THE CHEBYSHEV MOMENTS BY MEANS OF FORWARD
        !           RECURSION.
        !
        an = 4._SP
        DO i = 4, 13
          an2 = an*an
          v(i) = ((an2-4._SP)*(2._SP*(par22-an2-an2)*v(i-1)-ac)&
            +as-par2*(an+1._SP)*(an+2._SP)*v(i-2))&
            /(par2*(an-1._SP)*(an-2._SP))
          an = an + 2._SP
        END DO
      ELSE
        !
        !           COMPUTE THE CHEBYSHEV MOMENTS AS THE
        !           SOLUTIONS OF A BOUNDARY VALUE PROBLEM WITH 1
        !           INITIAL VALUE (V(3)) AND 1 END VALUE (COMPUTED
        !           USING AN ASYMPTOTIC FORMULA).
        !
        noequ = 25
        noeq1 = noequ - 1
        an = 6._SP
        d1(1) = 0._SP
        DO k = 1, noeq1
          an2 = an*an
          d(k) = -2._SP*(an2-4._SP)*(par22-an2-an2)
          d2(k) = (an-1._SP)*(an-2._SP)*par2
          d1(k+1) = (an+3._SP)*(an+4._SP)*par2
          v(k+3) = as - (an2-4._SP)*ac
          an = an + 2._SP
        END DO
        d2(noequ) = 0._SP
        an2 = an*an
        d(noequ) = -2._SP*(an2-4._SP)*(par22-an2-an2)
        v(noequ+3) = as - (an2-4._SP)*ac
        v(4) = v(4) - 0.56E+02_SP*par2*v(3)
        ass = parint*sinpar
        asap = (((((0.210E+03_SP*par2-1._SP)*cospar-(0.105E+03_SP*par2-0.63E+02_SP)&
          *ass)/an2-(1._SP-0.15E+02_SP*par2)*cospar+0.15E+02_SP*ass)&
          /an2-cospar+3._SP*ass)/an2-cospar)/an2
        v(noequ+3) = v(noequ+3) - 2._SP*asap*par2*(an-1._SP)*(an-2._SP)
        !
        !           SOLVE THE TRIDIAGONAL SYSTEM BY MEANS OF GAUSSIAN
        !           ELIMINATION WITH PARTIAL PIVOTING.
        !
        v2(1:25,1:1) => v(4:28)
        CALL SGTSV(noequ,1,d1,d,d2,v2,noequ,iers)
      END IF
      DO j = 1, 13
        Chebmo(m,2*j-1) = v(j)
      END DO
      !
      !           COMPUTE THE CHEBYSHEV MOMENTS WITH RESPECT TO SINE.
      !
      v(1) = 2._SP*(sinpar-parint*cospar)/par2
      v(2) = (18._SP-48._SP/par2)*sinpar/par2 + (-2._SP+48._SP/par2)&
        *cospar/parint
      ac = -24._SP*parint*cospar
      as = -8._SP*sinpar
      IF( ABS(parint)>24._SP ) THEN
        !
        !           COMPUTE THE CHEBYSHEV MOMENTS BY MEANS OF
        !           FORWARD RECURSION.
        !
        an = 3._SP
        DO i = 3, 12
          an2 = an*an
          v(i) = ((an2-4._SP)*(2._SP*(par22-an2-an2)*v(i-1)+as)&
            +ac-par2*(an+1._SP)*(an+2._SP)*v(i-2))&
            /(par2*(an-1._SP)*(an-2._SP))
          an = an + 2._SP
        END DO
      ELSE
        !
        !           COMPUTE THE CHEBYSHEV MOMENTS AS THE
        !           SOLUTIONS OF A BOUNDARY VALUE PROBLEM WITH 1
        !           INITIAL VALUE (V(2)) AND 1 END VALUE (COMPUTED
        !           USING AN ASYMPTOTIC FORMULA).
        !
        an = 5._SP
        DO k = 1, noeq1
          an2 = an*an
          d(k) = -2._SP*(an2-4._SP)*(par22-an2-an2)
          d2(k) = (an-1._SP)*(an-2._SP)*par2
          d1(k+1) = (an+3._SP)*(an+4._SP)*par2
          v(k+2) = ac + (an2-4._SP)*as
          an = an + 2._SP
        END DO
        an2 = an*an
        d(noequ) = -2._SP*(an2-4._SP)*(par22-an2-an2)
        v(noequ+2) = ac + (an2-4._SP)*as
        v(3) = v(3) - 0.42E+02_SP*par2*v(2)
        ass = parint*cospar
        asap = (((((0.105E+03_SP*par2-0.63E+02_SP)*ass+(0.210E+03_SP*par2-1._SP)*&
          sinpar)/an2+(0.15E+02_SP*par2-1._SP)*sinpar-0.15E+02_SP*ass)&
          /an2-3._SP*ass-sinpar)/an2-sinpar)/an2
        v(noequ+2) = v(noequ+2) - 2._SP*asap*par2*(an-1._SP)*(an-2._SP)
        !
        !           SOLVE THE TRIDIAGONAL SYSTEM BY MEANS OF GAUSSIAN
        !           ELIMINATION WITH PARTIAL PIVOTING.
        !
        v2(1:25,1:1) => v(3:27)
        CALL SGTSV(noequ,1,d1,d,d2,v2,noequ,iers)
      END IF
      DO j = 1, 12
        Chebmo(m,2*j) = v(j)
      END DO
    END IF
    IF( Nrmom<Momcom ) m = Nrmom + 1
    IF( Momcom<Maxp1-1 .AND. Nrmom>=Momcom ) Momcom = Momcom + 1
    !
    !           COMPUTE THE COEFFICIENTS OF THE CHEBYSHEV EXPANSIONS
    !           OF DEGREES 12 AND 24 OF THE FUNCTION F.
    !
    fval(1) = 0.5_SP*F(centr+hlgth)
    fval(13) = F(centr)
    fval(25) = 0.5_SP*F(centr-hlgth)
    DO i = 2, 12
      isym = 26 - i
      fval(i) = F(hlgth*x(i-1)+centr)
      fval(isym) = F(centr-hlgth*x(i-1))
    END DO
    CALL QCHEB(x,fval,cheb12,cheb24)
    !
    !           COMPUTE THE INTEGRAL AND ERROR ESTIMATES.
    !
    resc12 = cheb12(13)*Chebmo(m,13)
    ress12 = 0._SP
    k = 11
    DO j = 1, 6
      resc12 = resc12 + cheb12(k)*Chebmo(m,k)
      ress12 = ress12 + cheb12(k+1)*Chebmo(m,k+1)
      k = k - 2
    END DO
    resc24 = cheb24(25)*Chebmo(m,25)
    ress24 = 0._SP
    Resabs = ABS(cheb24(25))
    k = 23
    DO j = 1, 12
      resc24 = resc24 + cheb24(k)*Chebmo(m,k)
      ress24 = ress24 + cheb24(k+1)*Chebmo(m,k+1)
      Resabs = ABS(cheb24(k)) + ABS(cheb24(k+1))
      k = k - 2
    END DO
    estc = ABS(resc24-resc12)
    ests = ABS(ress24-ress12)
    Resabs = Resabs*ABS(hlgth)
    IF( Integr==2 ) THEN
      Result = conc*ress24 + cons*resc24
      Abserr = ABS(conc*ests) + ABS(cons*estc)
    ELSE
      Result = conc*resc24 - cons*ress24
      Abserr = ABS(conc*estc) + ABS(cons*ests)
    END IF
  ELSE
    CALL QK15W(F,QWGTF,Omega,p2,p3,p4,Integr,A,B,Result,Abserr,Resabs,Resasc)
    Neval = 15
  END IF
  !
END SUBROUTINE QC25F