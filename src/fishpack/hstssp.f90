!DECK HSTSSP
SUBROUTINE HSTSSP(A,B,M,Mbdcnd,Bda,Bdb,C,D,N,Nbdcnd,Bdc,Bdd,Elmbda,F,&
    Idimf,Pertrb,Ierror,W)
  IMPLICIT NONE
  REAL A, a1, a2, a3, B, Bda, Bdb, Bdc, Bdd, C, D, deltar, &
    deltht, dlrsq, dlthsq, dum, Elmbda, F, Pertrb, pi
  REAL PIMACH, W
  INTEGER i, Idimf, ierr1, Ierror, isw, iwb, iwc, iwr, iws, j, &
    jsw, k, lp, M, mb, Mbdcnd, mm1, N, Nbdcnd, np
  !***BEGIN PROLOGUE  HSTSSP
  !***PURPOSE  Solve the standard five-point finite difference
  !            approximation on a staggered grid to the Helmholtz
  !            equation in spherical coordinates and on the surface of
  !            the unit sphere (radius of 1).
  !***LIBRARY   SLATEC (FISHPACK)
  !***CATEGORY  I2B1A1A
  !***TYPE      SINGLE PRECISION (HSTSSP-S)
  !***KEYWORDS  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE, SPHERICAL
  !***AUTHOR  Adams, J., (NCAR)
  !           Swarztrauber, P. N., (NCAR)
  !           Sweet, R., (NCAR)
  !***DESCRIPTION
  !
  !     HSTSSP solves the standard five-point finite difference
  !     approximation on a staggered grid to the Helmholtz equation in
  !     spherical coordinates and on the surface of the unit sphere
  !     (radius of 1)
  !
  !             (1/SIN(THETA))(d/dTHETA)(SIN(THETA)(dU/dTHETA)) +
  !
  !       (1/SIN(THETA)**2)(d/dPHI)(dU/dPHI) + LAMBDA*U = F(THETA,PHI)
  !
  !     where THETA is colatitude and PHI is longitude.
  !
  !    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  !
  !    * * * * * * * *    Parameter Description     * * * * * * * * * *
  !
  !            * * * * * *   On Input    * * * * * *
  !
  !   A,B
  !     The range of THETA (colatitude), i.e. A .LE. THETA .LE. B.  A
  !     must be less than B and A must be non-negative.  A and B are in
  !     radians.  A = 0 corresponds to the north pole and B = PI
  !     corresponds to the south pole.
  !
  !
  !                  * * *  IMPORTANT  * * *
  !
  !     If B is equal to PI, then B must be computed using the statement
  !
  !     B = PIMACH(DUM)
  !
  !     This insures that B in the user's program is equal to PI in this
  !     program which permits several tests of the input parameters that
  !     otherwise would not be possible.
  !
  !                  * * * * * * * * * * * *
  !
  !
  !
  !   M
  !     The number of grid points in the interval (A,B).  The grid points
  !     in the THETA-direction are given by THETA(I) = A + (I-0.5)DTHETA
  !     for I=1,2,...,M where DTHETA =(B-A)/M.  M must be greater than 2.
  !
  !   MBDCND
  !     Indicates the type of boundary conditions at THETA = A and
  !     THETA = B.
  !
  !     = 1  If the solution is specified at THETA = A and THETA = B.
  !          (see note 3 below)
  !
  !     = 2  If the solution is specified at THETA = A and the derivative
  !          of the solution with respect to THETA is specified at
  !          THETA = B (see notes 2 and 3 below).
  !
  !     = 3  If the derivative of the solution with respect to THETA is
  !          specified at THETA = A (see notes 1, 2 below) and THETA = B.
  !
  !     = 4  If the derivative of the solution with respect to THETA is
  !          specified at THETA = A (see notes 1 and 2 below) and the
  !          solution is specified at THETA = B.
  !
  !     = 5  If the solution is unspecified at THETA = A = 0 and the
  !          solution is specified at THETA = B.  (see note 3 below)
  !
  !     = 6  If the solution is unspecified at THETA = A = 0 and the
  !          derivative of the solution with respect to THETA is
  !          specified at THETA = B (see note 2 below).
  !
  !     = 7  If the solution is specified at THETA = A and the
  !          solution is unspecified at THETA = B = PI. (see note 3 below)
  !
  !     = 8  If the derivative of the solution with respect to
  !          THETA is specified at THETA = A (see note 1 below)
  !          and the solution is unspecified at THETA = B = PI.
  !
  !     = 9  If the solution is unspecified at THETA = A = 0 and
  !          THETA = B = PI.
  !
  !     NOTES:  1.  If A = 0, do not use MBDCND = 3, 4, or 8,
  !                 but instead use MBDCND = 5, 6, or 9.
  !
  !             2.  If B = PI, do not use MBDCND = 2, 3, or 6,
  !                 but instead use MBDCND = 7, 8, or 9.
  !
  !             3.  When the solution is specified at THETA = 0 and/or
  !                 THETA = PI and the other boundary conditions are
  !                 combinations of unspecified, normal derivative, or
  !                 periodicity a singular system results.  The unique
  !                 solution is determined by extrapolation to the
  !                 specification of the solution at either THETA = 0 or
  !                 THETA = PI.  But in these cases the right side of the
  !                 system will be perturbed by the constant PERTRB.
  !
  !   BDA
  !     A one-dimensional array of length N that specifies the boundary
  !     values (if any) of the solution at THETA = A.  When
  !     MBDCND = 1, 2, or 7,
  !
  !              BDA(J) = U(A,PHI(J)),              J=1,2,...,N.
  !
  !     When MBDCND = 3, 4, or 8,
  !
  !              BDA(J) = (d/dTHETA)U(A,PHI(J)),    J=1,2,...,N.
  !
  !     When MBDCND has any other value, BDA is a dummy variable.
  !
  !   BDB
  !     A one-dimensional array of length N that specifies the boundary
  !     values of the solution at THETA = B.  When MBDCND = 1,4, or 5,
  !
  !              BDB(J) = U(B,PHI(J)),              J=1,2,...,N.
  !
  !     When MBDCND = 2,3, or 6,
  !
  !              BDB(J) = (d/dTHETA)U(B,PHI(J)),    J=1,2,...,N.
  !
  !     When MBDCND has any other value, BDB is a dummy variable.
  !
  !   C,D
  !     The range of PHI (longitude), i.e. C .LE. PHI .LE. D.
  !     C must be less than D.  If D-C = 2*PI, periodic boundary
  !     conditions are usually prescribed.
  !
  !   N
  !     The number of unknowns in the interval (C,D).  The unknowns in
  !     the PHI-direction are given by PHI(J) = C + (J-0.5)DPHI,
  !     J=1,2,...,N, where DPHI = (D-C)/N.  N must be greater than 2.
  !
  !   NBDCND
  !     Indicates the type of boundary conditions at PHI = C
  !     and PHI = D.
  !
  !     = 0  If the solution is periodic in PHI, i.e.
  !          U(I,J) = U(I,N+J).
  !
  !     = 1  If the solution is specified at PHI = C and PHI = D
  !          (see note below).
  !
  !     = 2  If the solution is specified at PHI = C and the derivative
  !          of the solution with respect to PHI is specified at
  !          PHI = D (see note below).
  !
  !     = 3  If the derivative of the solution with respect to PHI is
  !          specified at PHI = C and PHI = D.
  !
  !     = 4  If the derivative of the solution with respect to PHI is
  !          specified at PHI = C and the solution is specified at
  !          PHI = D (see note below).
  !
  !     NOTE:  When NBDCND = 1, 2, or 4, do not use MBDCND = 5, 6, 7, 8,
  !     or 9 (the former indicates that the solution is specified at
  !     a pole; the latter indicates the solution is unspecified).  Use
  !     instead MBDCND = 1 or 2.
  !
  !   BDC
  !     A one dimensional array of length M that specifies the boundary
  !     values of the solution at PHI = C.   When NBDCND = 1 or 2,
  !
  !              BDC(I) = U(THETA(I),C),              I=1,2,...,M.
  !
  !     When NBDCND = 3 or 4,
  !
  !              BDC(I) = (d/dPHI)U(THETA(I),C),       I=1,2,...,M.
  !
  !     When NBDCND = 0, BDC is a dummy variable.
  !
  !   BDD
  !     A one-dimensional array of length M that specifies the boundary
  !     values of the solution at PHI = D.  When NBDCND = 1 or 4,
  !
  !              BDD(I) = U(THETA(I),D),              I=1,2,...,M.
  !
  !     When NBDCND = 2 or 3,
  !
  !              BDD(I) = (d/dPHI)U(THETA(I),D),      I=1,2,...,M.
  !
  !     When NBDCND = 0, BDD is a dummy variable.
  !
  !   ELMBDA
  !     The constant LAMBDA in the Helmholtz equation.  If LAMBDA is
  !     greater than 0, a solution may not exist.  However, HSTSSP will
  !     attempt to find a solution.
  !
  !   F
  !     A two-dimensional array that specifies the values of the right
  !     side of the Helmholtz equation.  For I=1,2,...,M and J=1,2,...,N
  !
  !              F(I,J) = F(THETA(I),PHI(J)) .
  !
  !     F must be dimensioned at least M X N.
  !
  !   IDIMF
  !     The row (or first) dimension of the array F as it appears in the
  !     program calling HSTSSP.  This parameter is used to specify the
  !     variable dimension of F.  IDIMF must be at least M.
  !
  !   W
  !     A one-dimensional array that must be provided by the user for
  !     work space.  W may require up to 13M + 4N + M*INT(log2(N))
  !     locations.  The actual number of locations used is computed by
  !     HSTSSP and is returned in the location W(1).
  !
  !
  !            * * * * * *   On Output   * * * * * *
  !
  !   F
  !     Contains the solution U(I,J) of the finite difference
  !     approximation for the grid point (THETA(I),PHI(J)) for
  !     I=1,2,...,M, J=1,2,...,N.
  !
  !   PERTRB
  !     If a combination of periodic, derivative, or unspecified
  !     boundary conditions is specified for a Poisson equation
  !     (LAMBDA = 0), a solution may not exist.  PERTRB is a con-
  !     stant, calculated and subtracted from F, which ensures
  !     that a solution exists.  HSTSSP then computes this
  !     solution, which is a least squares solution to the
  !     original approximation.  This solution plus any constant is also
  !     a solution; hence, the solution is not unique.  The value of
  !     PERTRB should be small compared to the right side F.
  !     Otherwise, a solution is obtained to an essentially different
  !     problem.  This comparison should always be made to insure that
  !     a meaningful solution has been obtained.
  !
  !   IERROR
  !     An error flag that indicates invalid input parameters.
  !      Except for numbers 0 and 14, a solution is not attempted.
  !
  !     =  0  No error
  !
  !     =  1  A .LT. 0 or B .GT. PI
  !
  !     =  2  A .GE. B
  !
  !     =  3  MBDCND .LT. 1 or MBDCND .GT. 9
  !
  !     =  4  C .GE. D
  !
  !     =  5  N .LE. 2
  !
  !     =  6  NBDCND .LT. 0 or NBDCND .GT. 4
  !
  !     =  7  A .GT. 0 and MBDCND = 5, 6, or 9
  !
  !     =  8  A = 0 and MBDCND = 3, 4, or 8
  !
  !     =  9  B .LT. PI and MBDCND .GE. 7
  !
  !     = 10  B = PI and MBDCND = 2,3, or 6
  !
  !     = 11  MBDCND .GE. 5 and NDBCND = 1, 2, or 4
  !
  !     = 12  IDIMF .LT. M
  !
  !     = 13  M .LE. 2
  !
  !     = 14  LAMBDA .GT. 0
  !
  !     Since this is the only means of indicating a possibly
  !     incorrect call to HSTSSP, the user should test IERROR after
  !     the call.
  !
  !   W
  !     W(1) contains the required length of W.
  !
  ! *Long Description:
  !
  !    * * * * * * *   Program Specifications    * * * * * * * * * * * *
  !
  !    Dimension of   BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N),
  !    Arguments      W(see argument list)
  !
  !    Latest         June 1, 1977
  !    Revision
  !
  !    Subprograms    HSTSSP,POISTG,POSTG2,GENBUN,POISD2,POISN2,POISP2,
  !    Required       COSGEN,MERGE,TRIX,TRI3,PIMACH
  !
  !    Special        NONE
  !    Conditions
  !
  !    Common         NONE
  !    Blocks
  !
  !    I/O            NONE
  !
  !    Precision      Single
  !
  !    Specialist     Roland Sweet
  !
  !    Language       FORTRAN
  !
  !    History        Written by Roland Sweet at NCAR in April, 1977
  !
  !    Algorithm      This subroutine defines the finite-difference
  !                   equations, incorporates boundary data, adjusts the
  !                   right side when the system is singular and calls
  !                   either POISTG or GENBUN which solves the linear
  !                   system of equations.
  !
  !    Space          8427(decimal) = 20353(octal) locations on the
  !    Required       NCAR Control Data 7600
  !
  !     Timing and        The execution time T on the NCAR Control Data
  !     Accuracy       7600 for subroutine HSTSSP is roughly proportional
  !                    to M*N*log2(N).  Some typical values are listed in
  !                    the table below.
  !                       The solution process employed results in a loss
  !                    of no more than four significant digits for N and M
  !                    as large as 64.  More detailed information about
  !                    accuracy can be found in the documentation for
  !                    subroutine POISTG which is the routine that
  !                    actually solves the finite difference equations.
  !
  !
  !                       M(=N)    MBDCND    NBDCND    T(MSECS)
  !                       -----    ------    ------    --------
  !
  !                        32       1-9       1-4         56
  !                        64       1-9       1-4        230
  !
  !    Portability     American National Standards Institute FORTRAN.
  !                    The machine dependent constant PI is defined in
  !                    function PIMACH.
  !
  !    Required       COS
  !    Resident
  !    Routines
  !
  !    Reference      Schumann, U. and R. Sweet,'A Direct Method For
  !                   The Solution Of Poisson's Equation With Neumann
  !                   Boundary Conditions On A Staggered Grid Of
  !                   Arbitrary Size,' J. Comp. Phys. 20(1976),
  !                   pp. 171-182.
  !
  !    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  !
  !***REFERENCES  U. Schumann and R. Sweet, A direct method for the
  !                 solution of Poisson's equation with Neumann boundary
  !                 conditions on a staggered grid of arbitrary size,
  !                 Journal of Computational Physics 20, (1976),
  !                 pp. 171-182.
  !***ROUTINES CALLED  GENBUN, PIMACH, POISTG
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  HSTSSP
  !
  !
  DIMENSION F(Idimf,*), Bda(*), Bdb(*), Bdc(*), Bdd(*), W(*)
  !***FIRST EXECUTABLE STATEMENT  HSTSSP
  Ierror = 0
  pi = PIMACH(dum)
  IF ( A<0..OR.B>pi ) Ierror = 1
  IF ( A>=B ) Ierror = 2
  IF ( Mbdcnd<=0.OR.Mbdcnd>9 ) Ierror = 3
  IF ( C>=D ) Ierror = 4
  IF ( N<=2 ) Ierror = 5
  IF ( Nbdcnd<0.OR.Nbdcnd>=5 ) Ierror = 6
  IF ( A>0..AND.(Mbdcnd==5.OR.Mbdcnd==6.OR.Mbdcnd==9) ) Ierror = 7
  IF ( A==0..AND.(Mbdcnd==3.OR.Mbdcnd==4.OR.Mbdcnd==8) ) Ierror = 8
  IF ( B<pi.AND.Mbdcnd>=7 ) Ierror = 9
  IF ( B==pi.AND.(Mbdcnd==2.OR.Mbdcnd==3.OR.Mbdcnd==6) ) Ierror = 10
  IF ( Mbdcnd>=5.AND.(Nbdcnd==1.OR.Nbdcnd==2.OR.Nbdcnd==4) ) Ierror = 11
  IF ( Idimf<M ) Ierror = 12
  IF ( M<=2 ) Ierror = 13
  IF ( Ierror/=0 ) RETURN
  deltar = (B-A)/M
  dlrsq = deltar**2
  deltht = (D-C)/N
  dlthsq = deltht**2
  np = Nbdcnd + 1
  isw = 1
  jsw = 1
  mb = Mbdcnd
  IF ( Elmbda==0. ) THEN
    SELECT CASE (Mbdcnd)
      CASE (2)
        IF ( A==0. ) THEN
          mb = 6
          jsw = 2
        ENDIF
      CASE (3,6,8,9)
      CASE (4)
        IF ( B==pi ) THEN
          mb = 8
          jsw = 2
        ENDIF
      CASE DEFAULT
        IF ( A==0..AND.B==pi ) THEN
          mb = 9
          jsw = 2
        ENDIF
    END SELECT
  ENDIF
  !
  !     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
  !
  iwb = M
  iwc = iwb + M
  iwr = iwc + M
  iws = iwr + M
  DO i = 1, M
    j = iwr + i
    W(j) = SIN(A+(i-0.5)*deltar)
    W(i) = SIN((A+(i-1)*deltar))/dlrsq
  ENDDO
  mm1 = M - 1
  DO i = 1, mm1
    k = iwc + i
    W(k) = W(i+1)
    j = iwr + i
    k = iwb + i
    W(k) = Elmbda*W(j) - (W(i)+W(i+1))
  ENDDO
  W(iwr) = SIN(B)/dlrsq
  W(iwc) = Elmbda*W(iws) - (W(M)+W(iwr))
  DO i = 1, M
    j = iwr + i
    a1 = W(j)
    DO j = 1, N
      F(i,j) = a1*F(i,j)
    ENDDO
  ENDDO
  !
  !     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES.
  !
  SELECT CASE (mb)
    CASE (3,4,8)
      a1 = deltar*W(1)
      W(iwb+1) = W(iwb+1) + W(1)
      DO j = 1, N
        F(1,j) = F(1,j) + a1*Bda(j)
      ENDDO
    CASE (5,6,9)
    CASE DEFAULT
      a1 = 2.*W(1)
      W(iwb+1) = W(iwb+1) - W(1)
      DO j = 1, N
        F(1,j) = F(1,j) - a1*Bda(j)
      ENDDO
  END SELECT
  SELECT CASE (mb)
    CASE (2,3,6)
      a1 = deltar*W(iwr)
      W(iwc) = W(iwc) + W(iwr)
      DO j = 1, N
        F(M,j) = F(M,j) - a1*Bdb(j)
      ENDDO
    CASE (7,8,9)
    CASE DEFAULT
      a1 = 2.*W(iwr)
      W(iwc) = W(iwc) - W(iwr)
      DO j = 1, N
        F(M,j) = F(M,j) - a1*Bdb(j)
      ENDDO
  END SELECT
  !
  !     ENTER BOUNDARY DATA FOR PHI-BOUNDARIES.
  !
  a1 = 2./dlthsq
  SELECT CASE (np)
    CASE (1)
      GOTO 100
    CASE (4,5)
      a1 = 1./deltht
      DO i = 1, M
        j = iwr + i
        F(i,1) = F(i,1) + a1*Bdc(i)/W(j)
      ENDDO
    CASE DEFAULT
      DO i = 1, M
        j = iwr + i
        F(i,1) = F(i,1) - a1*Bdc(i)/W(j)
      ENDDO
  END SELECT
  a1 = 2./dlthsq
  SELECT CASE (np)
    CASE (1)
    CASE (3,4)
      a1 = 1./deltht
      DO i = 1, M
        j = iwr + i
        F(i,N) = F(i,N) - a1*Bdd(i)/W(j)
      ENDDO
    CASE DEFAULT
      DO i = 1, M
        j = iwr + i
        F(i,N) = F(i,N) - a1*Bdd(i)/W(j)
      ENDDO
  END SELECT
  !
  !     ADJUST RIGHT SIDE OF SINGULAR PROBLEMS TO INSURE EXISTENCE OF A
  !     SOLUTION.
  !
  100  Pertrb = 0.
  IF ( Elmbda<0 ) THEN
  ELSEIF ( Elmbda==0 ) THEN
    SELECT CASE (mb)
      CASE (1,2,4,5,7)
      CASE DEFAULT
        SELECT CASE (np)
          CASE (2,3,5)
          CASE DEFAULT
            isw = 2
            DO j = 1, N
              DO i = 1, M
                Pertrb = Pertrb + F(i,j)
              ENDDO
            ENDDO
            a1 = N*(COS(A)-COS(B))/(2.*SIN(0.5*deltar))
            Pertrb = Pertrb/a1
            DO i = 1, M
              j = iwr + i
              a1 = Pertrb*W(j)
              DO j = 1, N
                F(i,j) = F(i,j) - a1
              ENDDO
            ENDDO
            a2 = 0.
            a3 = 0.
            DO j = 1, N
              a2 = a2 + F(1,j)
              a3 = a3 + F(M,j)
            ENDDO
            a2 = a2/W(iwr+1)
            a3 = a3/W(iws)
        END SELECT
    END SELECT
  ELSE
    Ierror = 14
  ENDIF
  !
  !     MULTIPLY I-TH EQUATION THROUGH BY  R(I)*DELTHT**2
  !
  DO i = 1, M
    j = iwr + i
    a1 = dlthsq*W(j)
    W(i) = a1*W(i)
    j = iwc + i
    W(j) = a1*W(j)
    j = iwb + i
    W(j) = a1*W(j)
    DO j = 1, N
      F(i,j) = a1*F(i,j)
    ENDDO
  ENDDO
  lp = Nbdcnd
  W(1) = 0.
  W(iwr) = 0.
  !
  !     CALL POISTG OR GENBUN TO SOLVE THE SYSTEM OF EQUATIONS.
  !
  IF ( Nbdcnd==0 ) THEN
    CALL GENBUN(lp,N,1,M,W,W(iwb+1),W(iwc+1),Idimf,F,ierr1,W(iwr+1))
  ELSE
    CALL POISTG(lp,N,1,M,W,W(iwb+1),W(iwc+1),Idimf,F,ierr1,W(iwr+1))
  ENDIF
  W(1) = W(iwr+1) + 3*M
  IF ( isw==2.AND.jsw==2 ) THEN
    IF ( mb/=8 ) THEN
      a1 = 0.
      DO j = 1, N
        a1 = a1 + F(1,j)
      ENDDO
      a1 = (a1-dlrsq*a2/16.)/N
      IF ( Nbdcnd==3 ) a1 = a1 + (Bdd(1)-Bdc(1))/(D-C)
      a1 = Bda(1) - a1
    ELSE
      a1 = 0.
      DO j = 1, N
        a1 = a1 + F(M,j)
      ENDDO
      a1 = (a1-dlrsq*a3/16.)/N
      IF ( Nbdcnd==3 ) a1 = a1 + (Bdd(M)-Bdc(M))/(D-C)
      a1 = Bdb(1) - a1
    ENDIF
    DO i = 1, M
      DO j = 1, N
        F(i,j) = F(i,j) + a1
      ENDDO
    ENDDO
  ENDIF
END SUBROUTINE HSTSSP
