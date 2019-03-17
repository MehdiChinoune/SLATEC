!DECK HSTCYL
SUBROUTINE HSTCYL(A,B,M,Mbdcnd,Bda,Bdb,C,D,N,Nbdcnd,Bdc,Bdd,Elmbda,F,&
    Idimf,Pertrb,Ierror,W)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  HSTCYL
  !***PURPOSE  Solve the standard five-point finite difference
  !            approximation on a staggered grid to the modified
  !            Helmholtz equation in cylindrical coordinates.
  !***LIBRARY   SLATEC (FISHPACK)
  !***CATEGORY  I2B1A1A
  !***TYPE      SINGLE PRECISION (HSTCYL-S)
  !***KEYWORDS  CYLINDRICAL, ELLIPTIC, FISHPACK, HELMHOLTZ, PDE
  !***AUTHOR  Adams, J., (NCAR)
  !           Swarztrauber, P. N., (NCAR)
  !           Sweet, R., (NCAR)
  !***DESCRIPTION
  !
  !      HSTCYL solves the standard five-point finite difference
  !      approximation on a staggered grid to the modified Helmholtz
  !      equation in cylindrical coordinates
  !
  !          (1/R)(d/dR)(R(dU/dR)) + (d/dZ)(dU/dZ)C
  !                      + LAMBDA*(1/R**2)*U = F(R,Z)
  !
  !      This two-dimensional modified Helmholtz equation results
  !      from the Fourier transform of a three-dimensional Poisson
  !      equation.
  !
  !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  !
  !     * * * * * * * *    Parameter Description     * * * * * * * * * *
  !
  !             * * * * * *   On Input    * * * * * *
  !
  !    A,B
  !      The range of R, i.e. A .LE. R .LE. B.  A must be less than B and
  !      A must be non-negative.
  !
  !    M
  !      The number of grid points in the interval (A,B).  The grid points
  !      in the R-direction are given by R(I) = A + (I-0.5)DR for
  !      I=1,2,...,M where DR =(B-A)/M.  M must be greater than 2.
  !
  !    MBDCND
  !      Indicates the type of boundary conditions at R = A and R = B.
  !
  !      = 1  If the solution is specified at R = A (see note below) and
  !           R = B.
  !
  !      = 2  If the solution is specified at R = A (see note below) and
  !           the derivative of the solution with respect to R is
  !           specified at R = B.
  !
  !      = 3  If the derivative of the solution with respect to R is
  !           specified at R = A (see note below) and R = B.
  !
  !      = 4  If the derivative of the solution with respect to R is
  !           specified at R = A (see note below) and the solution is
  !           specified at R = B.
  !
  !      = 5  If the solution is unspecified at R = A = 0 and the solution
  !           is specified at R = B.
  !
  !      = 6  If the solution is unspecified at R = A = 0 and the
  !           derivative of the solution with respect to R is specified at
  !           R = B.
  !
  !      NOTE:  If A = 0, do not use MBDCND = 1,2,3, or 4, but instead
  !             use MBDCND = 5 or 6.  The resulting approximation gives
  !             the only meaningful boundary condition, i.e. dU/dR = 0.
  !             (see D. Greenspan, 'Introductory Numerical Analysis Of
  !             Elliptic Boundary Value Problems,' Harper and Row, 1965,
  !             Chapter 5.)
  !
  !    BDA
  !      A one-dimensional array of length N that specifies the boundary
  !      values (if any) of the solution at R = A.  When MBDCND = 1 or 2,
  !
  !               BDA(J) = U(A,Z(J)),          J=1,2,...,N.
  !
  !      When MBDCND = 3 or 4,
  !
  !               BDA(J) = (d/dR)U(A,Z(J)),    J=1,2,...,N.
  !
  !      When MBDCND = 5 or 6, BDA is a dummy variable.
  !
  !    BDB
  !      A one-dimensional array of length N that specifies the boundary
  !      values of the solution at R = B.  When MBDCND = 1,4, or 5,
  !
  !               BDB(J) = U(B,Z(J)),          J=1,2,...,N.
  !
  !      When MBDCND = 2,3, or 6,
  !
  !               BDB(J) = (d/dR)U(B,Z(J)),    J=1,2,...,N.
  !
  !    C,D
  !      The range of Z, i.e. C .LE. Z .LE. D.  C must be less
  !      than D.
  !
  !    N
  !      The number of unknowns in the interval (C,D).  The unknowns in
  !      the Z-direction are given by Z(J) = C + (J-0.5)DZ,
  !      J=1,2,...,N, where DZ = (D-C)/N.  N must be greater than 2.
  !
  !    NBDCND
  !      Indicates the type of boundary conditions at Z = C
  !      and Z = D.
  !
  !      = 0  If the solution is periodic in Z, i.e.
  !           U(I,J) = U(I,N+J).
  !
  !      = 1  If the solution is specified at Z = C and Z = D.
  !
  !      = 2  If the solution is specified at Z = C and the derivative
  !           of the solution with respect to Z is specified at
  !           Z = D.
  !
  !      = 3  If the derivative of the solution with respect to Z is
  !           specified at Z = C and Z = D.
  !
  !      = 4  If the derivative of the solution with respect to Z is
  !           specified at Z = C and the solution is specified at
  !           Z = D.
  !
  !    BDC
  !      A one dimensional array of length M that specifies the boundary
  !      values of the solution at Z = C.   When NBDCND = 1 or 2,
  !
  !               BDC(I) = U(R(I),C),              I=1,2,...,M.
  !
  !      When NBDCND = 3 or 4,
  !
  !               BDC(I) = (d/dZ)U(R(I),C),         I=1,2,...,M.
  !
  !      When NBDCND = 0, BDC is a dummy variable.
  !
  !    BDD
  !      A one-dimensional array of length M that specifies the boundary
  !      values of the solution at Z = D.  when NBDCND = 1 or 4,
  !
  !               BDD(I) = U(R(I),D),              I=1,2,...,M.
  !
  !      When NBDCND = 2 or 3,
  !
  !               BDD(I) = (d/dZ)U(R(I),D),        I=1,2,...,M.
  !
  !      When NBDCND = 0, BDD is a dummy variable.
  !
  !    ELMBDA
  !      The constant LAMBDA in the modified Helmholtz equation.  If
  !      LAMBDA is greater than 0, a solution may not exist.  However,
  !      HSTCYL will attempt to find a solution.  LAMBDA must be zero
  !      when MBDCND = 5 or 6.
  !
  !    F
  !      A two-dimensional array that specifies the values of the right
  !      side of the modified Helmholtz equation.  For I=1,2,...,M
  !      and J=1,2,...,N
  !
  !               F(I,J) = F(R(I),Z(J)) .
  !
  !      F must be dimensioned at least M X N.
  !
  !    IDIMF
  !      The row (or first) dimension of the array F as it appears in the
  !      program calling HSTCYL.  This parameter is used to specify the
  !      variable dimension of F.  IDIMF must be at least M.
  !
  !    W
  !      A one-dimensional array that must be provided by the user for
  !      work space.  W may require up to 13M + 4N + M*INT(log2(N))
  !      locations.  The actual number of locations used is computed by
  !      HSTCYL and is returned in the location W(1).
  !
  !
  !             * * * * * *   On Output   * * * * * *
  !
  !    F
  !      Contains the solution U(I,J) of the finite difference
  !      approximation for the grid point (R(I),Z(J)) for
  !      I=1,2,...,M, J=1,2,...,N.
  !
  !    PERTRB
  !      If a combination of periodic, derivative, or unspecified
  !      boundary conditions is specified for a Poisson equation
  !      (LAMBDA = 0), a solution may not exist.  PERTRB is a con-
  !      stant, calculated and subtracted from F, which ensures
  !      that a solution exists.  HSTCYL then computes this
  !      solution, which is a least squares solution to the
  !      original approximation.  This solution plus any constant is also
  !      a solution; hence, the solution is not unique.  The value of
  !      PERTRB should be small compared to the right side F.
  !      Otherwise, a solution is obtained to an essentially different
  !      problem.  This comparison should always be made to insure that
  !      a meaningful solution has been obtained.
  !
  !    IERROR
  !      An error flag that indicates invalid input parameters.
  !      Except for numbers 0 and 11, a solution is not attempted.
  !
  !      =  0  No error
  !
  !      =  1  A .LT. 0
  !
  !      =  2  A .GE. B
  !
  !      =  3  MBDCND .LT. 1 or MBDCND .GT. 6
  !
  !      =  4  C .GE. D
  !
  !      =  5  N .LE. 2
  !
  !      =  6  NBDCND .LT. 0 or NBDCND .GT. 4
  !
  !      =  7  A = 0 and MBDCND = 1,2,3, or 4
  !
  !      =  8  A .GT. 0 and MBDCND .GE. 5
  !
  !      =  9  M .LE. 2
  !
  !      = 10  IDIMF .LT. M
  !
  !      = 11  LAMBDA .GT. 0
  !
  !      = 12  A=0, MBDCND .GE. 5, ELMBDA .NE. 0
  !
  !      Since this is the only means of indicating a possibly
  !      incorrect call to HSTCYL, the user should test IERROR after
  !      the call.
  !
  !    W
  !      W(1) contains the required length of W.
  !
  ! *Long Description:
  !
  !     * * * * * * *   Program Specifications    * * * * * * * * * * * *
  !
  !     Dimension OF   BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N),
  !     Arguments      W(see argument list)
  !
  !     Latest         June 1, 1977
  !     Revision
  !
  !     Subprograms    HSTCYL,POISTG,POSTG2,GENBUN,POISD2,POISN2,POISP2,
  !     Required       COSGEN,MERGE,TRIX,TRI3,PIMACH
  !
  !     Special        NONE
  !     Conditions
  !
  !     Common         NONE
  !     Blocks
  !
  !     I/O            NONE
  !
  !     Precision      Single
  !
  !     Specialist     Roland Sweet
  !
  !     Language       FORTRAN
  !
  !     History        Written by Roland Sweet at NCAR in March, 1977
  !
  !     Algorithm      This subroutine defines the finite-difference
  !                    equations, incorporates boundary data, adjusts the
  !                    right side when the system is singular and calls
  !                    either POISTG or GENBUN which solves the linear
  !                    system of equations.
  !
  !     Space          8228(decimal) = 20044(octal) locations on the
  !     Required       NCAR Control Data 7600
  !
  !     Timing and        The execution time T on the NCAR Control Data
  !     Accuracy       7600 for subroutine HSTCYL is roughly proportional
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
  !                        32       1-6       1-4         56
  !                        64       1-6       1-4        230
  !
  !     Portability    American National Standards Institute Fortran.
  !                    The machine dependent constant PI is defined in
  !                    function PIMACH.
  !
  !     Required       COS
  !     Resident
  !     Routines
  !
  !     Reference      Schumann, U. and R. Sweet,'A Direct Method For
  !                    The Solution of Poisson's Equation With Neumann
  !                    Boundary Conditions On A Staggered Grid Of
  !                    Arbitrary Size,' J. Comp. Phys. 20(1976),
  !                    pp. 171-182.
  !
  !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  !
  !***REFERENCES  U. Schumann and R. Sweet, A direct method for the
  !                 solution of Poisson's equation with Neumann boundary
  !                 conditions on a staggered grid of arbitrary size,
  !                 Journal of Computational Physics 20, (1976),
  !                 pp. 171-182.
  !***ROUTINES CALLED  GENBUN, POISTG
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  HSTCYL
  REAL A, a1, B, Bda, Bdb, Bdc, Bdd, C, D, deltar, deltht, &
    dlrsq, dlthsq, Elmbda, F, Pertrb, W
  INTEGER i, Idimf, ierr1, Ierror, iwb, iwc, iwr, j, k, lp, M, &
    Mbdcnd, N, Nbdcnd, np
  DIMENSION F(Idimf,*), Bda(*), Bdb(*), Bdc(*), Bdd(*), W(*)
  !***FIRST EXECUTABLE STATEMENT  HSTCYL
  Ierror = 0
  IF ( A<0. ) Ierror = 1
  IF ( A>=B ) Ierror = 2
  IF ( Mbdcnd<=0.OR.Mbdcnd>=7 ) Ierror = 3
  IF ( C>=D ) Ierror = 4
  IF ( N<=2 ) Ierror = 5
  IF ( Nbdcnd<0.OR.Nbdcnd>=5 ) Ierror = 6
  IF ( A==0..AND.Mbdcnd/=5.AND.Mbdcnd/=6 ) Ierror = 7
  IF ( A>0..AND.Mbdcnd>=5 ) Ierror = 8
  IF ( Idimf<M ) Ierror = 10
  IF ( M<=2 ) Ierror = 9
  IF ( A==0..AND.Mbdcnd>=5.AND.Elmbda/=0. ) Ierror = 12
  IF ( Ierror/=0 ) RETURN
  deltar = (B-A)/M
  dlrsq = deltar**2
  deltht = (D-C)/N
  dlthsq = deltht**2
  np = Nbdcnd + 1
  !
  !     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
  !
  iwb = M
  iwc = iwb + M
  iwr = iwc + M
  DO i = 1, M
    j = iwr + i
    W(j) = A + (i-0.5)*deltar
    W(i) = (A+(i-1)*deltar)/(dlrsq*W(j))
    k = iwc + i
    W(k) = (A+i*deltar)/(dlrsq*W(j))
    k = iwb + i
    W(k) = Elmbda/W(j)**2 - 2./dlrsq
  ENDDO
  !
  !     ENTER BOUNDARY DATA FOR R-BOUNDARIES.
  !
  SELECT CASE (Mbdcnd)
    CASE (3,4)
      a1 = deltar*W(1)
      W(iwb+1) = W(iwb+1) + W(1)
      DO j = 1, N
        F(1,j) = F(1,j) + a1*Bda(j)
      ENDDO
    CASE (5,6)
    CASE DEFAULT
      a1 = 2.*W(1)
      W(iwb+1) = W(iwb+1) - W(1)
      DO j = 1, N
        F(1,j) = F(1,j) - a1*Bda(j)
      ENDDO
  END SELECT
  SELECT CASE (Mbdcnd)
    CASE (2,3,6)
      W(iwc) = W(iwc) + W(iwr)
      a1 = deltar*W(iwr)
      DO j = 1, N
        F(M,j) = F(M,j) - a1*Bdb(j)
      ENDDO
    CASE DEFAULT
      W(iwc) = W(iwc) - W(iwr)
      a1 = 2.*W(iwr)
      DO j = 1, N
        F(M,j) = F(M,j) - a1*Bdb(j)
      ENDDO
  END SELECT
  !
  !     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES.
  !
  a1 = 2./dlthsq
  SELECT CASE (np)
    CASE (1)
      GOTO 100
    CASE (4,5)
      a1 = 1./deltht
      DO i = 1, M
        F(i,1) = F(i,1) + a1*Bdc(i)
      ENDDO
    CASE DEFAULT
      DO i = 1, M
        F(i,1) = F(i,1) - a1*Bdc(i)
      ENDDO
  END SELECT
  a1 = 2./dlthsq
  SELECT CASE (np)
    CASE (1)
    CASE (3,4)
      a1 = 1./deltht
      DO i = 1, M
        F(i,N) = F(i,N) - a1*Bdd(i)
      ENDDO
    CASE DEFAULT
      DO i = 1, M
        F(i,N) = F(i,N) - a1*Bdd(i)
      ENDDO
  END SELECT
  !
  !     ADJUST RIGHT SIDE OF SINGULAR PROBLEMS TO INSURE EXISTENCE OF A
  !     SOLUTION.
  !
  100  Pertrb = 0.
  IF ( Elmbda<0 ) THEN
  ELSEIF ( Elmbda==0 ) THEN
    SELECT CASE (Mbdcnd)
      CASE (1,2,4,5)
      CASE DEFAULT
        SELECT CASE (np)
          CASE (2,3,5)
          CASE DEFAULT
            DO i = 1, M
              a1 = 0.
              DO j = 1, N
                a1 = a1 + F(i,j)
              ENDDO
              j = iwr + i
              Pertrb = Pertrb + a1*W(j)
            ENDDO
            Pertrb = Pertrb/(M*N*0.5*(A+B))
            DO i = 1, M
              DO j = 1, N
                F(i,j) = F(i,j) - Pertrb
              ENDDO
            ENDDO
        END SELECT
    END SELECT
  ELSE
    Ierror = 11
  ENDIF
  !
  !     MULTIPLY I-TH EQUATION THROUGH BY  DELTHT**2
  !
  DO i = 1, M
    W(i) = W(i)*dlthsq
    j = iwc + i
    W(j) = W(j)*dlthsq
    j = iwb + i
    W(j) = W(j)*dlthsq
    DO j = 1, N
      F(i,j) = F(i,j)*dlthsq
    ENDDO
  ENDDO
  lp = Nbdcnd
  W(1) = 0.
  W(iwr) = 0.
  !
  !     CALL GENBUN TO SOLVE THE SYSTEM OF EQUATIONS.
  !
  IF ( Nbdcnd==0 ) THEN
    CALL GENBUN(lp,N,1,M,W,W(iwb+1),W(iwc+1),Idimf,F,ierr1,W(iwr+1))
  ELSE
    CALL POISTG(lp,N,1,M,W,W(iwb+1),W(iwc+1),Idimf,F,ierr1,W(iwr+1))
  ENDIF
  W(1) = W(iwr+1) + 3*M
END SUBROUTINE HSTCYL
