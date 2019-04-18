!** HWSCYL
SUBROUTINE HWSCYL(A,B,M,Mbdcnd,Bda,Bdb,C,D,N,Nbdcnd,Bdc,Bdd,Elmbda,F,&
    Idimf,Pertrb,Ierror,W)
  !>
  !  Solve a standard finite difference approximation
  !            to the Helmholtz equation in cylindrical coordinates.
  !***
  ! **Library:**   SLATEC (FISHPACK)
  !***
  ! **Category:**  I2B1A1A
  !***
  ! **Type:**      SINGLE PRECISION (HWSCYL-S)
  !***
  ! **Keywords:**  CYLINDRICAL, ELLIPTIC, FISHPACK, HELMHOLTZ, PDE
  !***
  ! **Author:**  Adams, J., (NCAR)
  !           Swarztrauber, P. N., (NCAR)
  !           Sweet, R., (NCAR)
  !***
  ! **Description:**
  !
  !     Subroutine HWSCYL solves a finite difference approximation to the
  !     Helmholtz equation in cylindrical coordinates:
  !
  !          (1/R)(d/dR)(R(dU/dR)) + (d/dZ)(dU/dZ)
  !
  !                                + (LAMBDA/R**2)U = F(R,Z)
  !
  !     This modified Helmholtz equation results from the Fourier
  !     transform of the three-dimensional Poisson equation.
  !
  !     * * * * * * * *    Parameter Description     * * * * * * * * * *
  !
  !             * * * * * *   On Input    * * * * * *
  !
  !     A,B
  !       The range of R, i.e., A .LE. R .LE. B.  A must be less than B
  !       and A must be non-negative.
  !
  !     M
  !       The number of panels into which the interval (A,B) is
  !       subdivided.  Hence, there will be M+1 grid points in the
  !       R-direction given by R(I) = A+(I-1)DR, for I = 1,2,...,M+1,
  !       where DR = (B-A)/M is the panel width. M must be greater than 3.
  !
  !     MBDCND
  !       Indicates the type of boundary conditions at R = A and R = B.
  !
  !       = 1  If the solution is specified at R = A and R = B.
  !       = 2  If the solution is specified at R = A and the derivative of
  !            the solution with respect to R is specified at R = B.
  !       = 3  If the derivative of the solution with respect to R is
  !            specified at R = A (see note below) and R = B.
  !       = 4  If the derivative of the solution with respect to R is
  !            specified at R = A (see note below) and the solution is
  !            specified at R = B.
  !       = 5  If the solution is unspecified at R = A = 0 and the
  !            solution is specified at R = B.
  !       = 6  If the solution is unspecified at R = A = 0 and the
  !            derivative of the solution with respect to R is specified
  !            at R = B.
  !
  !       NOTE:  If A = 0, do not use MBDCND = 3 or 4, but instead use
  !              MBDCND = 1,2,5, or 6  .
  !
  !     BDA
  !       A one-dimensional array of length N+1 that specifies the values
  !       of the derivative of the solution with respect to R at R = A.
  !       When MBDCND = 3 or 4,
  !
  !            BDA(J) = (d/dR)U(A,Z(J)), J = 1,2,...,N+1  .
  !
  !       When MBDCND has any other value, BDA is a dummy variable.
  !
  !     BDB
  !       A one-dimensional array of length N+1 that specifies the values
  !       of the derivative of the solution with respect to R at R = B.
  !       When MBDCND = 2,3, or 6,
  !
  !            BDB(J) = (d/dR)U(B,Z(J)), J = 1,2,...,N+1  .
  !
  !       When MBDCND has any other value, BDB is a dummy variable.
  !
  !     C,D
  !       The range of Z, i.e., C .LE. Z .LE. D.  C must be less than D.
  !
  !     N
  !       The number of panels into which the interval (C,D) is
  !       subdivided.  Hence, there will be N+1 grid points in the
  !       Z-direction given by Z(J) = C+(J-1)DZ, for J = 1,2,...,N+1,
  !       where DZ = (D-C)/N is the panel width. N must be greater than 3.
  !
  !     NBDCND
  !       Indicates the type of boundary conditions at Z = C and Z = D.
  !
  !       = 0  If the solution is periodic in Z, i.e., U(I,1) = U(I,N+1).
  !       = 1  If the solution is specified at Z = C and Z = D.
  !       = 2  If the solution is specified at Z = C and the derivative of
  !            the solution with respect to Z is specified at Z = D.
  !       = 3  If the derivative of the solution with respect to Z is
  !            specified at Z = C and Z = D.
  !       = 4  If the derivative of the solution with respect to Z is
  !            specified at Z = C and the solution is specified at Z = D.
  !
  !     BDC
  !       A one-dimensional array of length M+1 that specifies the values
  !       of the derivative of the solution with respect to Z at Z = C.
  !       When NBDCND = 3 or 4,
  !
  !            BDC(I) = (d/dZ)U(R(I),C), I = 1,2,...,M+1  .
  !
  !       When NBDCND has any other value, BDC is a dummy variable.
  !
  !     BDD
  !       A one-dimensional array of length M+1 that specifies the values
  !       of the derivative of the solution with respect to Z at Z = D.
  !       When NBDCND = 2 or 3,
  !
  !            BDD(I) = (d/dZ)U(R(I),D), I = 1,2,...,M+1  .
  !
  !       When NBDCND has any other value, BDD is a dummy variable.
  !
  !     ELMBDA
  !       The constant LAMBDA in the Helmholtz equation.  If
  !       LAMBDA .GT. 0, a solution may not exist.  However, HWSCYL will
  !       attempt to find a solution.  LAMBDA must be zero when
  !       MBDCND = 5 or 6  .
  !
  !     F
  !       A two-dimensional array that specifies the values of the right
  !       side of the Helmholtz equation and boundary data (if any).  For
  !       I = 2,3,...,M and J = 2,3,...,N
  !
  !            F(I,J) = F(R(I),Z(J)).
  !
  !       On the boundaries F is defined by
  !
  !            MBDCND   F(1,J)            F(M+1,J)
  !            ------   ---------         ---------
  !
  !              1      U(A,Z(J))         U(B,Z(J))
  !              2      U(A,Z(J))         F(B,Z(J))
  !              3      F(A,Z(J))         F(B,Z(J))   J = 1,2,...,N+1
  !              4      F(A,Z(J))         U(B,Z(J))
  !              5      F(0,Z(J))         U(B,Z(J))
  !              6      F(0,Z(J))         F(B,Z(J))
  !
  !            NBDCND   F(I,1)            F(I,N+1)
  !            ------   ---------         ---------
  !
  !              0      F(R(I),C)         F(R(I),C)
  !              1      U(R(I),C)         U(R(I),D)
  !              2      U(R(I),C)         F(R(I),D)   I = 1,2,...,M+1
  !              3      F(R(I),C)         F(R(I),D)
  !              4      F(R(I),C)         U(R(I),D)
  !
  !       F must be dimensioned at least (M+1)*(N+1).
  !
  !       NOTE
  !
  !       If the table calls for both the solution U and the right side F
  !       at a corner then the solution must be specified.
  !
  !     IDIMF
  !       The row (or first) dimension of the array F as it appears in the
  !       program calling HWSCYL.  This parameter is used to specify the
  !       variable dimension of F.  IDIMF must be at least M+1  .
  !
  !     W
  !       A one-dimensional array that must be provided by the user for
  !       work space.  W may require up to 4*(N+1) +
  !       (13 + INT(log2(N+1)))*(M+1) locations.  The actual number of
  !       locations used is computed by HWSCYL and is returned in location
  !       W(1).
  !
  !
  !             * * * * * *   On Output     * * * * * *
  !
  !     F
  !       Contains the solution U(I,J) of the finite difference
  !       approximation for the grid point (R(I),Z(J)), I = 1,2,...,M+1,
  !       J = 1,2,...,N+1  .
  !
  !     PERTRB
  !       If one specifies a combination of periodic, derivative, and
  !       unspecified boundary conditions for a Poisson equation
  !       (LAMBDA = 0), a solution may not exist.  PERTRB is a constant,
  !       calculated and subtracted from F, which ensures that a solution
  !       exists.  HWSCYL then computes this solution, which is a least
  !       squares solution to the original approximation.  This solution
  !       plus any constant is also a solution.  Hence, the solution is
  !       not unique.  The value of PERTRB should be small compared to the
  !       right side F.  Otherwise, a solution is obtained to an
  !       essentially different problem.  This comparison should always
  !       be made to insure that a meaningful solution has been obtained.
  !
  !     IERROR
  !       An error flag which indicates invalid input parameters.  Except
  !       for numbers 0 and 11, a solution is not attempted.
  !
  !       =  0  No error.
  !       =  1  A .LT. 0  .
  !       =  2  A .GE. B.
  !       =  3  MBDCND .LT. 1 or MBDCND .GT. 6  .
  !       =  4  C .GE. D.
  !       =  5  N .LE. 3
  !       =  6  NBDCND .LT. 0 or NBDCND .GT. 4  .
  !       =  7  A = 0, MBDCND = 3 or 4  .
  !       =  8  A .GT. 0, MBDCND .GE. 5  .
  !       =  9  A = 0, LAMBDA .NE. 0, MBDCND .GE. 5  .
  !       = 10  IDIMF .LT. M+1  .
  !       = 11  LAMBDA .GT. 0  .
  !       = 12  M .LE. 3
  !
  !       Since this is the only means of indicating a possibly incorrect
  !       call to HWSCYL, the user should test IERROR after the call.
  !
  !     W
  !       W(1) contains the required length of W.
  !
  !- Long Description:
  !
  !     * * * * * * *   Program Specifications    * * * * * * * * * * * *
  !
  !     Dimension of   BDA(N+1),BDB(N+1),BDC(M+1),BDD(M+1),F(IDIMF,N+1),
  !     Arguments      W(see argument list)
  !
  !     Latest         June 1, 1976
  !     Revision
  !
  !     Subprograms    HWSCYL,GENBUN,POISD2,POISN2,POISP2,COSGEN,MERGE,
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
  !     Space          5818(decimal) = 13272(octal) locations on the NCAR
  !     Required       Control Data 7600
  !
  !     Timing and        The execution time T on the NCAR Control Data
  !     Accuracy       7600 for subroutine HWSCYL is roughly proportional
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
  !                        32        1         0          31
  !                        32        1         1          23
  !                        32        3         3          36
  !                        64        1         0         128
  !                        64        1         1          96
  !                        64        3         3         142
  !
  !     Portability    American National Standards Institute FORTRAN.
  !                    The machine dependent constant PI is defined in
  !                    function PIMACH.
  !
  !     Required       COS
  !     Resident
  !     Routines
  !
  !     Reference      Swarztrauber, P. and R. Sweet, 'Efficient FORTRAN
  !                    Subprograms for the Solution of Elliptic Equations'
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

  INTEGER i, id2, id3, id4, id5, id6, Idimf, ierr1, Ierror, ij, &
    istart, j, k, l, M, Mbdcnd, mp1, mstart, mstop, munk
  REAL A, a1, a2, B, Bda(*), Bdb(*), Bdc(*), Bdd(*), C, D, deltar, deltht, &
    dlrby2, dlrsq, dlthsq, Elmbda, F(Idimf,*), Pertrb, r, s
  REAL s1, s2, W(*)
  INTEGER N, Nbdcnd, np, np1, nsp1, nstart, nstm1, nstop, nunk
  !* FIRST EXECUTABLE STATEMENT  HWSCYL
  Ierror = 0
  IF ( A<0. ) Ierror = 1
  IF ( A>=B ) Ierror = 2
  IF ( Mbdcnd<=0.OR.Mbdcnd>=7 ) Ierror = 3
  IF ( C>=D ) Ierror = 4
  IF ( N<=3 ) Ierror = 5
  IF ( Nbdcnd<=-1.OR.Nbdcnd>=5 ) Ierror = 6
  IF ( A==0..AND.(Mbdcnd==3.OR.Mbdcnd==4) ) Ierror = 7
  IF ( A>0..AND.Mbdcnd>=5 ) Ierror = 8
  IF ( A==0..AND.Elmbda/=0..AND.Mbdcnd>=5 ) Ierror = 9
  IF ( Idimf<M+1 ) Ierror = 10
  IF ( M<=3 ) Ierror = 12
  IF ( Ierror/=0 ) RETURN
  mp1 = M + 1
  deltar = (B-A)/M
  dlrby2 = deltar/2.
  dlrsq = deltar**2
  np1 = N + 1
  deltht = (D-C)/N
  dlthsq = deltht**2
  np = Nbdcnd + 1
  !
  !     DEFINE RANGE OF INDICES I AND J FOR UNKNOWNS U(I,J).
  !
  mstart = 2
  mstop = M
  SELECT CASE (Mbdcnd)
    CASE (1)
    CASE (2)
      mstop = mp1
    CASE (3,6)
      mstart = 1
      mstop = mp1
    CASE DEFAULT
      mstart = 1
  END SELECT
  munk = mstop - mstart + 1
  nstart = 1
  nstop = N
  SELECT CASE (np)
    CASE (1,5)
    CASE (3)
      nstart = 2
      nstop = np1
    CASE (4)
      nstop = np1
    CASE DEFAULT
      nstart = 2
  END SELECT
  nunk = nstop - nstart + 1
  !
  !     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
  !
  id2 = munk
  id3 = id2 + munk
  id4 = id3 + munk
  id5 = id4 + munk
  id6 = id5 + munk
  istart = 1
  a1 = 2./dlrsq
  ij = 0
  IF ( Mbdcnd==3.OR.Mbdcnd==4 ) ij = 1
  IF ( Mbdcnd>4 ) THEN
    W(1) = 0.
    W(id2+1) = -2.*a1
    W(id3+1) = 2.*a1
    istart = 2
    ij = 1
  END IF
  DO i = istart, munk
    r = A + (i-ij)*deltar
    j = id5 + i
    W(j) = r
    j = id6 + i
    W(j) = 1./r**2
    W(i) = (r-dlrby2)/(r*dlrsq)
    j = id3 + i
    W(j) = (r+dlrby2)/(r*dlrsq)
    k = id6 + i
    j = id2 + i
    W(j) = -a1 + Elmbda*W(k)
  END DO
  SELECT CASE (Mbdcnd)
    CASE (1,5)
    CASE (3,6)
      W(id2) = a1
      W(id3+1) = a1*istart
    CASE (4)
      W(id3+1) = a1*istart
    CASE DEFAULT
      W(id2) = a1
  END SELECT
  !
  !     ENTER BOUNDARY DATA FOR R-BOUNDARIES.
  !
  SELECT CASE (Mbdcnd)
    CASE (3,4)
      a1 = 2.*deltar*W(1)
      DO j = nstart, nstop
        F(1,j) = F(1,j) + a1*Bda(j)
      END DO
    CASE (5,6)
    CASE DEFAULT
      a1 = W(1)
      DO j = nstart, nstop
        F(2,j) = F(2,j) - a1*F(1,j)
      END DO
  END SELECT
  SELECT CASE (Mbdcnd)
    CASE (2,3,6)
      a1 = 2.*deltar*W(id4)
      DO j = nstart, nstop
        F(mp1,j) = F(mp1,j) - a1*Bdb(j)
      END DO
    CASE DEFAULT
      a1 = W(id4)
      DO j = nstart, nstop
        F(M,j) = F(M,j) - a1*F(mp1,j)
      END DO
  END SELECT
  !
  !     ENTER BOUNDARY DATA FOR Z-BOUNDARIES.
  !
  a1 = 1./dlthsq
  l = id5 - mstart + 1
  SELECT CASE (np)
    CASE (1)
      GOTO 100
    CASE (4,5)
      a1 = 2./deltht
      DO i = mstart, mstop
        F(i,1) = F(i,1) + a1*Bdc(i)
      END DO
    CASE DEFAULT
      DO i = mstart, mstop
        F(i,2) = F(i,2) - a1*F(i,1)
      END DO
  END SELECT
  a1 = 1./dlthsq
  SELECT CASE (np)
    CASE (1)
    CASE (3,4)
      a1 = 2./deltht
      DO i = mstart, mstop
        F(i,np1) = F(i,np1) - a1*Bdd(i)
      END DO
    CASE DEFAULT
      DO i = mstart, mstop
        F(i,N) = F(i,N) - a1*F(i,np1)
      END DO
  END SELECT
  !
  !     ADJUST RIGHT SIDE OF SINGULAR PROBLEMS TO INSURE EXISTENCE OF A
  !     SOLUTION.
  !
  100  Pertrb = 0.
  IF ( Elmbda<0 ) THEN
  ELSEIF ( Elmbda==0 ) THEN
    W(id5+1) = .5*(W(id5+2)-dlrby2)
    SELECT CASE (Mbdcnd)
      CASE (1,2,4,5)
        GOTO 200
      CASE (3)
      CASE DEFAULT
        W(id5+1) = .5*W(id5+1)
    END SELECT
    SELECT CASE (np)
      CASE (1)
        a2 = 1.
      CASE (2,3,5)
        GOTO 200
      CASE DEFAULT
        a2 = 2.
    END SELECT
    k = id5 + munk
    W(k) = .5*(W(k-1)+dlrby2)
    s = 0.
    DO i = mstart, mstop
      s1 = 0.
      nsp1 = nstart + 1
      nstm1 = nstop - 1
      DO j = nsp1, nstm1
        s1 = s1 + F(i,j)
      END DO
      k = i + l
      s = s + (a2*s1+F(i,nstart)+F(i,nstop))*W(k)
    END DO
    s2 = M*A + (.75+(M-1)*(M+1))*dlrby2
    IF ( Mbdcnd==3 ) s2 = s2 + .25*dlrby2
    s1 = (2.+a2*(nunk-2))*s2
    Pertrb = s/s1
    DO i = mstart, mstop
      DO j = nstart, nstop
        F(i,j) = F(i,j) - Pertrb
      END DO
    END DO
  ELSE
    Ierror = 11
  END IF
  !
  !     MULTIPLY I-TH EQUATION THROUGH BY DELTHT**2 TO PUT EQUATION INTO
  !     CORRECT FORM FOR SUBROUTINE GENBUN.
  !
  200 CONTINUE
  DO i = mstart, mstop
    k = i - mstart + 1
    W(k) = W(k)*dlthsq
    j = id2 + k
    W(j) = W(j)*dlthsq
    j = id3 + k
    W(j) = W(j)*dlthsq
    DO j = nstart, nstop
      F(i,j) = F(i,j)*dlthsq
    END DO
  END DO
  W(1) = 0.
  W(id4) = 0.
  !
  !     CALL GENBUN TO SOLVE THE SYSTEM OF EQUATIONS.
  !
  CALL GENBUN(Nbdcnd,nunk,1,munk,W(1),W(id2+1),W(id3+1),Idimf,&
    F(mstart,nstart),ierr1,W(id4+1))
  W(1) = W(id4+1) + 3*munk
  IF ( Nbdcnd==0 ) THEN
    DO i = mstart, mstop
      F(i,np1) = F(i,1)
    END DO
  END IF
END SUBROUTINE HWSCYL
