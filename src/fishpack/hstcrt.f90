!DECK HSTCRT
SUBROUTINE HSTCRT(A,B,M,Mbdcnd,Bda,Bdb,C,D,N,Nbdcnd,Bdc,Bdd,Elmbda,F,&
    Idimf,Pertrb,Ierror,W)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  HSTCRT
  !***PURPOSE  Solve the standard five-point finite difference
  !            approximation on a staggered grid to the Helmholtz equation
  !            in Cartesian coordinates.
  !***LIBRARY   SLATEC (FISHPACK)
  !***CATEGORY  I2B1A1A
  !***TYPE      SINGLE PRECISION (HSTCRT-S)
  !***KEYWORDS  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE
  !***AUTHOR  Adams, J., (NCAR)
  !           Swarztrauber, P. N., (NCAR)
  !           Sweet, R., (NCAR)
  !***DESCRIPTION
  !
  !      HSTCRT solves the standard five-point finite difference
  !      approximation on a staggered grid to the Helmholtz equation in
  !      Cartesian coordinates
  !
  !      (d/dX)(dU/dX) + (d/dY)(dU/dY) + LAMBDA*U = F(X,Y)
  !
  !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  !
  !     * * * * * * * *    Parameter Description     * * * * * * * * * *
  !
  !             * * * * * *   On Input    * * * * * *
  !
  !    A,B
  !      The range of X, i.e. A .LE. X .LE. B.  A must be less than B.
  !
  !    M
  !      The number of grid points in the interval (A,B).  The grid points
  !      in the X-direction are given by X(I) = A + (I-0.5)dX for
  !      I=1,2,...,M where dX =(B-A)/M.  M must be greater than 2.
  !
  !    MBDCND
  !      Indicates the type of boundary conditions at X = A and X = B.
  !
  !      = 0  If the solution is periodic in X,
  !           U(M+I,J) = U(I,J).
  !
  !      = 1  If the solution is specified at X = A and X = B.
  !
  !      = 2  If the solution is specified at X = A and the derivative
  !           of the solution with respect to X is specified at X = B.
  !
  !      = 3  If the derivative of the solution with respect to X is
  !           specified at X = A  and X = B.
  !
  !      = 4  If the derivative of the solution with respect to X is
  !           specified at X = A  and the solution is specified at X = B.
  !
  !    BDA
  !      A one-dimensional array of length N that specifies the boundary
  !      values (if any) of the solution at X = A.  When MBDCND = 1 or 2,
  !
  !               BDA(J) = U(A,Y(J)),          J=1,2,...,N.
  !
  !      When MBDCND = 3 or 4,
  !
  !               BDA(J) = (d/dX)U(A,Y(J)),    J=1,2,...,N.
  !
  !    BDB
  !      A one-dimensional array of length N that specifies the boundary
  !      values of the solution at X = B.  When MBDCND = 1 or 4
  !
  !               BDB(J) = U(B,Y(J)),          J=1,2,...,N.
  !
  !      When MBDCND = 2 or 3
  !
  !               BDB(J) = (d/dX)U(B,Y(J)),    J=1,2,...,N.
  !
  !    C,D
  !      The range of Y, i.e. C .LE. Y .LE. D.  C must be less
  !      than D.
  !
  !    N
  !      The number of unknowns in the interval (C,D).  The unknowns in
  !      the Y-direction are given by Y(J) = C + (J-0.5)DY,
  !      J=1,2,...,N, where DY = (D-C)/N.  N must be greater than 2.
  !
  !    NBDCND
  !      Indicates the type of boundary conditions at Y = C
  !      and Y = D.
  !
  !      = 0  If the solution is periodic in Y, i.e.
  !           U(I,J) = U(I,N+J).
  !
  !      = 1  If the solution is specified at Y = C and Y = D.
  !
  !      = 2  If the solution is specified at Y = C and the derivative
  !           of the solution with respect to Y is specified at Y = D.
  !
  !      = 3  If the derivative of the solution with respect to Y is
  !           specified at Y = C and Y = D.
  !
  !      = 4  If the derivative of the solution with respect to Y is
  !           specified at Y = C and the solution is specified at Y = D.
  !
  !    BDC
  !      A one dimensional array of length M that specifies the boundary
  !      values of the solution at Y = C.   When NBDCND = 1 or 2,
  !
  !               BDC(I) = U(X(I),C),              I=1,2,...,M.
  !
  !      When NBDCND = 3 or 4,
  !
  !               BDC(I) = (d/dY)U(X(I),C),     I=1,2,...,M.
  !
  !      When NBDCND = 0, BDC is a dummy variable.
  !
  !    BDD
  !      A one-dimensional array of length M that specifies the boundary
  !      values of the solution at Y = D.  When NBDCND = 1 or 4,
  !
  !               BDD(I) = U(X(I),D),              I=1,2,...,M.
  !
  !      When NBDCND = 2 or 3,
  !
  !               BDD(I) = (d/dY)U(X(I),D),    I=1,2,...,M.
  !
  !      When NBDCND = 0, BDD is a dummy variable.
  !
  !    ELMBDA
  !      The constant LAMBDA in the Helmholtz equation.  If LAMBDA is
  !      greater than 0, a solution may not exist.  However, HSTCRT will
  !      attempt to find a solution.
  !
  !    F
  !      A two-dimensional array that specifies the values of the right
  !      side of the Helmholtz equation.  For I=1,2,...,M and J=1,2,...,N
  !
  !               F(I,J) = F(X(I),Y(J)) .
  !
  !      F must be dimensioned at least M X N.
  !
  !    IDIMF
  !      The row (or first) dimension of the array F as it appears in the
  !      program calling HSTCRT.  This parameter is used to specify the
  !      variable dimension of F.  IDIMF must be at least M.
  !
  !    W
  !      A one-dimensional array that must be provided by the user for
  !      work space.  W may require up to 13M + 4N + M*INT(log2(N))
  !      locations.  The actual number of locations used is computed by
  !      HSTCRT and is returned in the location W(1).
  !
  !
  !             * * * * * *   On Output   * * * * * *
  !
  !    F
  !      Contains the solution U(I,J) of the finite difference
  !      approximation for the grid point (X(I),Y(J)) for
  !      I=1,2,...,M, J=1,2,...,N.
  !
  !    PERTRB
  !      If a combination of periodic or derivative boundary conditions is
  !      specified for a Poisson equation (LAMBDA = 0), a solution may not
  !      exist.  PERTRB is a constant, calculated and subtracted from F,
  !      which ensures that a solution exists.  HSTCRT then computes this
  !      solution, which is a least squares solution to the original
  !      approximation.  This solution plus any constant is also a
  !      solution; hence, the solution is not unique.  The value of PERTRB
  !      should be small compared to the right side F.  Otherwise, a
  !      solution is obtained to an essentially different problem.  This
  !      comparison should always be made to insure that a meaningful
  !      solution has been obtained.
  !
  !    IERROR
  !      An error flag that indicates invalid input parameters.
  !       Except for numbers 0 and  6, a solution is not attempted.
  !
  !      =  0  No error
  !
  !      =  1  A .GE. B
  !
  !      =  2  MBDCND .LT. 0 or MBDCND .GT. 4
  !
  !      =  3  C .GE. D
  !
  !      =  4  N .LE. 2
  !
  !      =  5  NBDCND .LT. 0 or NBDCND .GT. 4
  !
  !      =  6  LAMBDA .GT. 0
  !
  !      =  7  IDIMF .LT. M
  !
  !      =  8  M .LE. 2
  !
  !      Since this is the only means of indicating a possibly
  !      incorrect call to HSTCRT, the user should test IERROR after
  !      the call.
  !
  !    W
  !      W(1) contains the required length of W.
  !
  ! *Long Description:
  !
  !     * * * * * * *   Program Specifications    * * * * * * * * * * * *
  !
  !     Dimension of   BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N),
  !     Arguments      W(See argument list)
  !
  !     Latest         June 1, 1977
  !     Revision
  !
  !     Subprograms    HSTCRT,POISTG,POSTG2,GENBUN,POISD2,POISN2,POISP2,
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
  !     History        Written by Roland Sweet at NCAR in January, 1977
  !
  !     Algorithm      This subroutine defines the finite-difference
  !                    equations, incorporates boundary data, adjusts the
  !                    right side when the system is singular and calls
  !                    either POISTG or GENBUN which solves the linear
  !                    system of equations.
  !
  !     Space          8131(decimal) = 17703(octal) locations on the
  !     Required       NCAR Control Data 7600
  !
  !     Timing and        The execution time T on the NCAR Control Data
  !     Accuracy       7600 for subroutine HSTCRT is roughly proportional
  !                    to M*N*log2(N).  Some typical values are listed in
  !                    the table below.
  !                       The solution process employed results in a loss
  !                    of no more than FOUR significant digits for N and M
  !                    as large as 64.  More detailed information about
  !                    accuracy can be found in the documentation for
  !                    subroutine POISTG which is the routine that
  !                    actually solves the finite difference equations.
  !
  !
  !                       M(=N)    MBDCND    NBDCND    T(MSECS)
  !                       -----    ------    ------    --------
  !
  !                        32       1-4       1-4         56
  !                        64       1-4       1-4        230
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
  !                    The Solution Of Poisson's Equation With Neumann
  !                    Boundary Conditions On A Staggered Grid Of
  !                    Arbitrary Size,' J. COMP. PHYS. 20(1976),
  !                    PP. 171-182.
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
  !***END PROLOGUE  HSTCRT
  REAL A, B, Bda, Bdb, Bdc, Bdd, C, D, deltax, deltay, delxsq, &
    delysq, Elmbda, F, Pertrb, s, st2, twdelx, twdely, twdysq
  REAL W
  INTEGER i, id2, id3, id4, Idimf, ierr1, Ierror, j, M, Mbdcnd, &
    mp, mperod, N, Nbdcnd, np, nperod
  DIMENSION F(Idimf,*), Bda(*), Bdb(*), Bdc(*), Bdd(*), W(*)
  !***FIRST EXECUTABLE STATEMENT  HSTCRT
  Ierror = 0
  IF ( A>=B ) Ierror = 1
  IF ( Mbdcnd<0.OR.Mbdcnd>4 ) Ierror = 2
  IF ( C>=D ) Ierror = 3
  IF ( N<=2 ) Ierror = 4
  IF ( Nbdcnd<0.OR.Nbdcnd>4 ) Ierror = 5
  IF ( Idimf<M ) Ierror = 7
  IF ( M<=2 ) Ierror = 8
  IF ( Ierror/=0 ) RETURN
  nperod = Nbdcnd
  mperod = 0
  IF ( Mbdcnd>0 ) mperod = 1
  deltax = (B-A)/M
  twdelx = 1./deltax
  delxsq = 2./deltax**2
  deltay = (D-C)/N
  twdely = 1./deltay
  delysq = deltay**2
  twdysq = 2./delysq
  np = Nbdcnd + 1
  mp = Mbdcnd + 1
  !
  !     DEFINE THE A,B,C COEFFICIENTS IN W-ARRAY.
  !
  id2 = M
  id3 = id2 + M
  id4 = id3 + M
  s = (deltay/deltax)**2
  st2 = 2.*s
  DO i = 1, M
    W(i) = s
    j = id2 + i
    W(j) = -st2 + Elmbda*delysq
    j = id3 + i
    W(j) = s
  ENDDO
  !
  !     ENTER BOUNDARY DATA FOR X-BOUNDARIES.
  !
  SELECT CASE (mp)
    CASE (1)
      GOTO 100
    CASE (4,5)
      DO j = 1, N
        F(1,j) = F(1,j) + Bda(j)*twdelx
      ENDDO
      W(id2+1) = W(id2+1) + W(1)
    CASE DEFAULT
      DO j = 1, N
        F(1,j) = F(1,j) - Bda(j)*delxsq
      ENDDO
      W(id2+1) = W(id2+1) - W(1)
  END SELECT
  SELECT CASE (mp)
    CASE (1)
    CASE (3,4)
      DO j = 1, N
        F(M,j) = F(M,j) - Bdb(j)*twdelx
      ENDDO
      W(id3) = W(id3) + W(1)
    CASE DEFAULT
      DO j = 1, N
        F(M,j) = F(M,j) - Bdb(j)*delxsq
      ENDDO
      W(id3) = W(id3) - W(1)
  END SELECT
  !
  !     ENTER BOUNDARY DATA FOR Y-BOUNDARIES.
  !
  100 CONTINUE
  SELECT CASE (np)
    CASE (1)
      GOTO 200
    CASE (4,5)
      DO i = 1, M
        F(i,1) = F(i,1) + Bdc(i)*twdely
      ENDDO
    CASE DEFAULT
      DO i = 1, M
        F(i,1) = F(i,1) - Bdc(i)*twdysq
      ENDDO
  END SELECT
  SELECT CASE (np)
    CASE (1)
    CASE (3,4)
      DO i = 1, M
        F(i,N) = F(i,N) - Bdd(i)*twdely
      ENDDO
    CASE DEFAULT
      DO i = 1, M
        F(i,N) = F(i,N) - Bdd(i)*twdysq
      ENDDO
  END SELECT
  200 CONTINUE
  DO i = 1, M
    DO j = 1, N
      F(i,j) = F(i,j)*delysq
    ENDDO
  ENDDO
  IF ( mperod/=0 ) THEN
    W(1) = 0.
    W(id4) = 0.
  ENDIF
  Pertrb = 0.
  IF ( Elmbda<0 ) THEN
  ELSEIF ( Elmbda==0 ) THEN
    SELECT CASE (mp)
      CASE (2,3,5)
      CASE DEFAULT
        SELECT CASE (np)
          CASE (2,3,5)
          CASE DEFAULT
            !
            !     FOR SINGULAR PROBLEMS MUST ADJUST DATA TO INSURE THAT A SOLUTION
            !     WILL EXIST.
            !
            s = 0.
            DO j = 1, N
              DO i = 1, M
                s = s + F(i,j)
              ENDDO
            ENDDO
            Pertrb = s/(M*N)
            DO j = 1, N
              DO i = 1, M
                F(i,j) = F(i,j) - Pertrb
              ENDDO
            ENDDO
            Pertrb = Pertrb/delysq
        END SELECT
    END SELECT
  ELSE
    Ierror = 6
  ENDIF
  !
  !     SOLVE THE EQUATION.
  !
  IF ( nperod==0 ) THEN
    CALL GENBUN(nperod,N,mperod,M,W(1),W(id2+1),W(id3+1),Idimf,F,ierr1,&
      W(id4+1))
  ELSE
    CALL POISTG(nperod,N,mperod,M,W(1),W(id2+1),W(id3+1),Idimf,F,ierr1,&
      W(id4+1))
  ENDIF
  W(1) = W(id4+1) + 3*M
END SUBROUTINE HSTCRT
