!** HWSCRT
SUBROUTINE HWSCRT(A,B,M,Mbdcnd,Bda,Bdb,C,D,N,Nbdcnd,Bdc,Bdd,Elmbda,F,&
    Idimf,Pertrb,Ierror,W)
  !>
  !  Solves the standard five-point finite difference
  !            approximation to the Helmholtz equation in Cartesian
  !            coordinates.
  !***
  ! **Library:**   SLATEC (FISHPACK)
  !***
  ! **Category:**  I2B1A1A
  !***
  ! **Type:**      SINGLE PRECISION (HWSCRT-S)
  !***
  ! **Keywords:**  CARTESIAN, ELLIPTIC, FISHPACK, HELMHOLTZ, PDE
  !***
  ! **Author:**  Adams, J., (NCAR)
  !           Swarztrauber, P. N., (NCAR)
  !           Sweet, R., (NCAR)
  !***
  ! **Description:**
  !
  !     Subroutine HWSCRT solves the standard five-point finite
  !     difference approximation to the Helmholtz equation in Cartesian
  !     coordinates:
  !
  !          (d/dX)(dU/dX) + (d/dY)(dU/dY) + LAMBDA*U = F(X,Y).
  !
  !
  !
  !     * * * * * * * *    Parameter Description     * * * * * * * * * *
  !
  !             * * * * * *   On Input    * * * * * *
  !
  !     A,B
  !       The range of X, i.e., A .LE. X .LE. B.  A must be less than B.
  !
  !     M
  !       The number of panels into which the interval (A,B) is
  !       subdivided.  Hence, there will be M+1 grid points in the
  !       X-direction given by X(I) = A+(I-1)DX for I = 1,2,...,M+1,
  !       where DX = (B-A)/M is the panel width. M must be greater than 3.
  !
  !     MBDCND
  !       Indicates the type of boundary conditions at X = A and X = B.
  !
  !       = 0  If the solution is periodic in X, i.e., U(I,J) = U(M+I,J).
  !       = 1  If the solution is specified at X = A and X = B.
  !       = 2  If the solution is specified at X = A and the derivative of
  !            the solution with respect to X is specified at X = B.
  !       = 3  If the derivative of the solution with respect to X is
  !            specified at X = A and X = B.
  !       = 4  If the derivative of the solution with respect to X is
  !            specified at X = A and the solution is specified at X = B.
  !
  !     BDA
  !       A one-dimensional array of length N+1 that specifies the values
  !       of the derivative of the solution with respect to X at X = A.
  !       When MBDCND = 3 or 4,
  !
  !            BDA(J) = (d/dX)U(A,Y(J)), J = 1,2,...,N+1  .
  !
  !       When MBDCND has any other value, BDA is a dummy variable.
  !
  !     BDB
  !       A one-dimensional array of length N+1 that specifies the values
  !       of the derivative of the solution with respect to X at X = B.
  !       When MBDCND = 2 or 3,
  !
  !            BDB(J) = (d/dX)U(B,Y(J)), J = 1,2,...,N+1  .
  !
  !       When MBDCND has any other value BDB is a dummy variable.
  !
  !     C,D
  !       The range of Y, i.e., C .LE. Y .LE. D.  C must be less than D.
  !
  !     N
  !       The number of panels into which the interval (C,D) is
  !       subdivided.  Hence, there will be N+1 grid points in the
  !       Y-direction given by Y(J) = C+(J-1)DY for J = 1,2,...,N+1, where
  !       DY = (D-C)/N is the panel width.  N must be greater than 3.
  !
  !     NBDCND
  !       Indicates the type of boundary conditions at Y = C and Y = D.
  !
  !       = 0  If the solution is periodic in Y, i.e., U(I,J) = U(I,N+J).
  !       = 1  If the solution is specified at Y = C and Y = D.
  !       = 2  If the solution is specified at Y = C and the derivative of
  !            the solution with respect to Y is specified at Y = D.
  !       = 3  If the derivative of the solution with respect to Y is
  !            specified at Y = C and Y = D.
  !       = 4  If the derivative of the solution with respect to Y is
  !            specified at Y = C and the solution is specified at Y = D.
  !
  !     BDC
  !       A one-dimensional array of length M+1 that specifies the values
  !       of the derivative of the solution with respect to Y at Y = C.
  !       When NBDCND = 3 or 4,
  !
  !            BDC(I) = (d/dY)U(X(I),C), I = 1,2,...,M+1  .
  !
  !       When NBDCND has any other value, BDC is a dummy variable.
  !
  !     BDD
  !       A one-dimensional array of length M+1 that specifies the values
  !       of the derivative of the solution with respect to Y at Y = D.
  !       When NBDCND = 2 or 3,
  !
  !            BDD(I) = (d/dY)U(X(I),D), I = 1,2,...,M+1  .
  !
  !       When NBDCND has any other value, BDD is a dummy variable.
  !
  !     ELMBDA
  !       The constant LAMBDA in the Helmholtz equation.  If
  !       LAMBDA .GT. 0, a solution may not exist.  However, HWSCRT will
  !       attempt to find a solution.
  !
  !     F
  !       A two-dimensional array which specifies the values of the right
  !       side of the Helmholtz equation and boundary values (if any).
  !       For I = 2,3,...,M and J = 2,3,...,N
  !
  !            F(I,J) = F(X(I),Y(J)).
  !
  !       On the boundaries F is defined by
  !
  !            MBDCND     F(1,J)        F(M+1,J)
  !            ------     ---------     --------
  !
  !              0        F(A,Y(J))     F(A,Y(J))
  !              1        U(A,Y(J))     U(B,Y(J))
  !              2        U(A,Y(J))     F(B,Y(J))     J = 1,2,...,N+1
  !              3        F(A,Y(J))     F(B,Y(J))
  !              4        F(A,Y(J))     U(B,Y(J))
  !
  !
  !            NBDCND     F(I,1)        F(I,N+1)
  !            ------     ---------     --------
  !
  !              0        F(X(I),C)     F(X(I),C)
  !              1        U(X(I),C)     U(X(I),D)
  !              2        U(X(I),C)     F(X(I),D)     I = 1,2,...,M+1
  !              3        F(X(I),C)     F(X(I),D)
  !              4        F(X(I),C)     U(X(I),D)
  !
  !       F must be dimensioned at least (M+1)*(N+1).
  !
  !       NOTE:
  !
  !       If the table calls for both the solution U and the right side F
  !       at a corner then the solution must be specified.
  !
  !     IDIMF
  !       The row (or first) dimension of the array F as it appears in the
  !       program calling HWSCRT.  This parameter is used to specify the
  !       variable dimension of F.  IDIMF must be at least M+1  .
  !
  !     W
  !       A one-dimensional array that must be provided by the user for
  !       work space.  W may require up to 4*(N+1) +
  !       (13 + INT(log2(N+1)))*(M+1) locations.  The actual number of
  !       locations used is computed by HWSCRT and is returned in location
  !       W(1).
  !
  !
  !             * * * * * *   On Output     * * * * * *
  !
  !     F
  !       Contains the solution U(I,J) of the finite difference
  !       approximation for the grid point (X(I),Y(J)), I = 1,2,...,M+1,
  !       J = 1,2,...,N+1  .
  !
  !     PERTRB
  !       If a combination of periodic or derivative boundary conditions
  !       is specified for a Poisson equation (LAMBDA = 0), a solution may
  !       not exist.  PERTRB is a constant, calculated and subtracted from
  !       F, which ensures that a solution exists.  HWSCRT then computes
  !       this solution, which is a least squares solution to the original
  !       approximation.  This solution plus any constant is also a
  !       solution.  Hence, the solution is not unique.  The value of
  !       PERTRB should be small compared to the right side F.  Otherwise,
  !       a solution is obtained to an essentially different problem.
  !       This comparison should always be made to insure that a
  !       meaningful solution has been obtained.
  !
  !     IERROR
  !       An error flag that indicates invalid input parameters.  Except
  !       for numbers 0 and 6, a solution is not attempted.
  !
  !       = 0  No error.
  !       = 1  A .GE. B.
  !       = 2  MBDCND .LT. 0 or MBDCND .GT. 4  .
  !       = 3  C .GE. D.
  !       = 4  N .LE. 3
  !       = 5  NBDCND .LT. 0 or NBDCND .GT. 4  .
  !       = 6  LAMBDA .GT. 0  .
  !       = 7  IDIMF .LT. M+1  .
  !       = 8  M .LE. 3
  !
  !       Since this is the only means of indicating a possibly incorrect
  !       call to HWSCRT, the user should test IERROR after the call.
  !
  !     W
  !       W(1) contains the required length of W.
  !
  !- Long Description:
  !
  !     * * * * * * *   Program Specifications    * * * * * * * * * * * *
  !
  !
  !     Dimension of   BDA(N+1),BDB(N+1),BDC(M+1),BDD(M+1),F(IDIMF,N+1),
  !     Arguments      W(see argument list)
  !
  !     Latest         June 1, 1976
  !     Revision
  !
  !     Subprograms    HWSCRT,GENBUN,POISD2,POISN2,POISP2,COSGEN,MERGE,
  !     Required       TRIX,TRI3,PIMACH
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
  !     History        Standardized September 1, 1973
  !                    Revised April 1, 1976
  !
  !     Algorithm      The routine defines the finite difference
  !                    equations, incorporates boundary data, and adjusts
  !                    the right side of singular systems and then calls
  !                    GENBUN to solve the system.
  !
  !     Space          13110(octal) = 5704(decimal) locations on the NCAR
  !     Required       Control Data 7600
  !
  !     Timing and        The execution time T on the NCAR Control Data
  !     Accuracy       7600 for subroutine HWSCRT is roughly proportional
  !                    to M*N*log2(N), but also depends on the input
  !                    parameters NBDCND and MBDCND.  Some typical values
  !                    are listed in the table below.
  !                       The solution process employed results in a loss
  !                    of no more than three significant digits for N and
  !                    M as large as 64.  More detailed information about
  !                    accuracy can be found in the documentation for
  !                    subroutine GENBUN which is the routine that
  !                    solves the finite difference equations.
  !
  !
  !                       M(=N)    MBDCND    NBDCND    T(MSECS)
  !                       -----    ------    ------    --------
  !
  !                        32        0         0          31
  !                        32        1         1          23
  !                        32        3         3          36
  !                        64        0         0         128
  !                        64        1         1          96
  !                        64        3         3         142
  !
  !     Portability    American National Standards Institute FORTRAN.
  !                    The machine dependent constant PI is defined in
  !                    function PIMACH.
  !
  !     Reference      Swarztrauber, P. and R. Sweet, 'Efficient FORTRAN
  !                    Subprograms for The Solution Of Elliptic Equations'
  !                    NCAR TN/IA-109, July, 1975, 138 pp.
  !
  !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  !
  !***
  ! **References:**  P. N. Swarztrauber and R. Sweet, Efficient Fortran
  !                 subprograms for the solution of elliptic equations,
  !                 NCAR TN/IA-109, July 1975, 138 pp.
  !***
  ! **Routines called:**  GENBUN

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER i, id2, id3, id4, Idimf, ierr1, Ierror, j, M, Mbdcnd, &
    mp, mp1, mperod, mskip, msp1, mstart, mstm1, mstop, munk, N
  REAL A, a1, a2, B, Bda(*), Bdb(*), Bdc(*), Bdd(*), C, D, deltax, deltay, &
    delxsq, delysq, Elmbda, F(Idimf,*), Pertrb, s, s1, st2
  REAL twdelx, twdely, W(*)
  INTEGER Nbdcnd, np, np1, nperod, nskip, nsp1, nstart, nstm1, nstop, nunk
  !* FIRST EXECUTABLE STATEMENT  HWSCRT
  Ierror = 0
  IF ( A>=B ) Ierror = 1
  IF ( Mbdcnd<0.OR.Mbdcnd>4 ) Ierror = 2
  IF ( C>=D ) Ierror = 3
  IF ( N<=3 ) Ierror = 4
  IF ( Nbdcnd<0.OR.Nbdcnd>4 ) Ierror = 5
  IF ( Idimf<M+1 ) Ierror = 7
  IF ( M<=3 ) Ierror = 8
  IF ( Ierror/=0 ) RETURN
  nperod = Nbdcnd
  mperod = 0
  IF ( Mbdcnd>0 ) mperod = 1
  deltax = (B-A)/M
  twdelx = 2./deltax
  delxsq = 1./deltax**2
  deltay = (D-C)/N
  twdely = 2./deltay
  delysq = 1./deltay**2
  np = Nbdcnd + 1
  np1 = N + 1
  mp = Mbdcnd + 1
  mp1 = M + 1
  nstart = 1
  nstop = N
  nskip = 1
  SELECT CASE (np)
    CASE (1,5)
    CASE (3)
      nstart = 2
      GOTO 100
    CASE (4)
      GOTO 100
    CASE DEFAULT
      nstart = 2
  END SELECT
  GOTO 200
  100  nstop = np1
  nskip = 2
  200  nunk = nstop - nstart + 1
  !
  !     ENTER BOUNDARY DATA FOR X-BOUNDARIES.
  !
  mstart = 1
  mstop = M
  mskip = 1
  SELECT CASE (mp)
    CASE (1)
      GOTO 500
    CASE (3)
      mstart = 2
      mstop = mp1
      mskip = 2
    CASE (4)
      mstop = mp1
      mskip = 2
      GOTO 300
    CASE (5)
      GOTO 300
    CASE DEFAULT
      mstart = 2
  END SELECT
  DO j = nstart, nstop
    F(2,j) = F(2,j) - F(1,j)*delxsq
  END DO
  GOTO 400
  300 CONTINUE
  DO j = nstart, nstop
    F(1,j) = F(1,j) + Bda(j)*twdelx
  END DO
  400 CONTINUE
  IF ( mskip==2 ) THEN
    DO j = nstart, nstop
      F(mp1,j) = F(mp1,j) - Bdb(j)*twdelx
    END DO
  ELSE
    DO j = nstart, nstop
      F(M,j) = F(M,j) - F(mp1,j)*delxsq
    END DO
  END IF
  500  munk = mstop - mstart + 1
  !
  !     ENTER BOUNDARY DATA FOR Y-BOUNDARIES.
  !
  SELECT CASE (np)
    CASE (1)
      GOTO 600
    CASE (4,5)
      DO i = mstart, mstop
        F(i,1) = F(i,1) + Bdc(i)*twdely
      END DO
    CASE DEFAULT
      DO i = mstart, mstop
        F(i,2) = F(i,2) - F(i,1)*delysq
      END DO
  END SELECT
  IF ( nskip==2 ) THEN
    DO i = mstart, mstop
      F(i,np1) = F(i,np1) - Bdd(i)*twdely
    END DO
  ELSE
    DO i = mstart, mstop
      F(i,N) = F(i,N) - F(i,np1)*delysq
    END DO
  END IF
  !
  !    MULTIPLY RIGHT SIDE BY DELTAY**2.
  !
  600  delysq = deltay*deltay
  DO i = mstart, mstop
    DO j = nstart, nstop
      F(i,j) = F(i,j)*delysq
    END DO
  END DO
  !
  !     DEFINE THE A,B,C COEFFICIENTS IN W-ARRAY.
  !
  id2 = munk
  id3 = id2 + munk
  id4 = id3 + munk
  s = delysq*delxsq
  st2 = 2.*s
  DO i = 1, munk
    W(i) = s
    j = id2 + i
    W(j) = -st2 + Elmbda*delysq
    j = id3 + i
    W(j) = s
  END DO
  IF ( mp/=1 ) THEN
    W(1) = 0.
    W(id4) = 0.
  END IF
  SELECT CASE (mp)
    CASE (1,2)
    CASE (4)
      W(id2) = st2
      W(id3+1) = st2
    CASE (5)
      W(id3+1) = st2
    CASE DEFAULT
      W(id2) = st2
  END SELECT
  Pertrb = 0.
  IF ( Elmbda<0 ) THEN
  ELSEIF ( Elmbda==0 ) THEN
    IF ( (Nbdcnd==0.OR.Nbdcnd==3).AND.(Mbdcnd==0.OR.Mbdcnd==3) ) THEN
      !
      !     FOR SINGULAR PROBLEMS MUST ADJUST DATA TO INSURE THAT A SOLUTION
      !     WILL EXIST.
      !
      a1 = 1.
      a2 = 1.
      IF ( Nbdcnd==3 ) a2 = 2.
      IF ( Mbdcnd==3 ) a1 = 2.
      s1 = 0.
      msp1 = mstart + 1
      mstm1 = mstop - 1
      nsp1 = nstart + 1
      nstm1 = nstop - 1
      DO j = nsp1, nstm1
        s = 0.
        DO i = msp1, mstm1
          s = s + F(i,j)
        END DO
        s1 = s1 + s*a1 + F(mstart,j) + F(mstop,j)
      END DO
      s1 = a2*s1
      s = 0.
      DO i = msp1, mstm1
        s = s + F(i,nstart) + F(i,nstop)
      END DO
      s1 = s1 + s*a1 + F(mstart,nstart) + F(mstart,nstop) + F(mstop,nstart)&
        + F(mstop,nstop)
      s = (2.+(nunk-2)*a2)*(2.+(munk-2)*a1)
      Pertrb = s1/s
      DO j = nstart, nstop
        DO i = mstart, mstop
          F(i,j) = F(i,j) - Pertrb
        END DO
      END DO
      Pertrb = Pertrb/delysq
    END IF
  ELSE
    Ierror = 6
  END IF
  !
  !     SOLVE THE EQUATION.
  !
  CALL GENBUN(nperod,nunk,mperod,munk,W(1),W(id2+1),W(id3+1),Idimf,&
    F(mstart,nstart),ierr1,W(id4+1))
  W(1) = W(id4+1) + 3*munk
  !
  !     FILL IN IDENTICAL VALUES WHEN HAVE PERIODIC BOUNDARY CONDITIONS.
  !
  IF ( Nbdcnd==0 ) THEN
    DO i = mstart, mstop
      F(i,np1) = F(i,1)
    END DO
  END IF
  IF ( Mbdcnd==0 ) THEN
    DO j = nstart, nstop
      F(mp1,j) = F(1,j)
    END DO
    IF ( Nbdcnd==0 ) F(mp1,np1) = F(1,np1)
  END IF
END SUBROUTINE HWSCRT
