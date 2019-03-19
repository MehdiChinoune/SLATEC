!** HSTPLR
SUBROUTINE HSTPLR(A,B,M,Mbdcnd,Bda,Bdb,C,D,N,Nbdcnd,Bdc,Bdd,Elmbda,F,&
    Idimf,Pertrb,Ierror,W)
  IMPLICIT NONE
  !>
  !***
  !  Solve the standard five-point finite difference
  !            approximation on a staggered grid to the Helmholtz equation
  !            in polar coordinates.
  !***
  ! **Library:**   SLATEC (FISHPACK)
  !***
  ! **Category:**  I2B1A1A
  !***
  ! **Type:**      SINGLE PRECISION (HSTPLR-S)
  !***
  ! **Keywords:**  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE, POLAR
  !***
  ! **Author:**  Adams, J., (NCAR)
  !           Swarztrauber, P. N., (NCAR)
  !           Sweet, R., (NCAR)
  !***
  ! **Description:**
  !
  !      HSTPLR solves the standard five-point finite difference
  !      approximation on a staggered grid to the Helmholtz equation in
  !      polar coordinates
  !
  !      (1/R)(d/DR)(R(dU/DR)) + (1/R**2)(d/dTHETA)(dU/dTHETA)
  !
  !                      + LAMBDA*U = F(R,THETA)
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
  !      = 1  If the solution is specified at R = A and R = B.
  !
  !      = 2  If the solution is specified at R = A and the derivative
  !           of the solution with respect to R is specified at R = B.
  !           (see note 1 below)
  !
  !      = 3  If the derivative of the solution with respect to R is
  !           specified at R = A (see note 2 below) and R = B.
  !
  !      = 4  If the derivative of the solution with respect to R is
  !           specified at R = A (see note 2 below) and the solution is
  !           specified at R = B.
  !
  !      = 5  If the solution is unspecified at R = A = 0 and the solution
  !           is specified at R = B.
  !
  !      = 6  If the solution is unspecified at R = A = 0 and the
  !           derivative of the solution with respect to R is specified at
  !           R = B.
  !
  !      NOTE 1:  If A = 0, MBDCND = 2, and NBDCND = 0 or 3, the system of
  !               equations to be solved is singular.  The unique solution
  !               is determined by extrapolation to the specification of
  !               U(0,THETA(1)).  But in this case the right side of the
  !               system will be perturbed by the constant PERTRB.
  !
  !      NOTE 2:  If A = 0, do not use MBDCND = 3 or 4, but instead use
  !               MBDCND = 1,2,5, or 6.
  !
  !    BDA
  !      A one-dimensional array of length N that specifies the boundary
  !      values (if any) of the solution at R = A.  When MBDCND = 1 or 2,
  !
  !               BDA(J) = U(A,THETA(J)),          J=1,2,...,N.
  !
  !      When MBDCND = 3 or 4,
  !
  !               BDA(J) = (d/dR)U(A,THETA(J)),    J=1,2,...,N.
  !
  !      When MBDCND = 5 or 6, BDA is a dummy variable.
  !
  !    BDB
  !      A one-dimensional array of length N that specifies the boundary
  !      values of the solution at R = B.  When MBDCND = 1,4, or 5,
  !
  !               BDB(J) = U(B,THETA(J)),          J=1,2,...,N.
  !
  !      When MBDCND = 2,3, or 6,
  !
  !               BDB(J) = (d/dR)U(B,THETA(J)),    J=1,2,...,N.
  !
  !    C,D
  !      The range of THETA, i.e. C .LE. THETA .LE. D.  C must be less
  !      than D.
  !
  !    N
  !      The number of unknowns in the interval (C,D).  The unknowns in
  !      the THETA-direction are given by THETA(J) = C + (J-0.5)DT,
  !      J=1,2,...,N, where DT = (D-C)/N.  N must be greater than 2.
  !
  !    NBDCND
  !      Indicates the type of boundary conditions at THETA = C
  !      and THETA = D.
  !
  !      = 0  If the solution is periodic in THETA, i.e.
  !           U(I,J) = U(I,N+J).
  !
  !      = 1  If the solution is specified at THETA = C and THETA = D
  !           (see note below).
  !
  !      = 2  If the solution is specified at THETA = C and the derivative
  !           of the solution with respect to THETA is specified at
  !           THETA = D (see note below).
  !
  !      = 3  If the derivative of the solution with respect to THETA is
  !           specified at THETA = C and THETA = D.
  !
  !      = 4  If the derivative of the solution with respect to THETA is
  !           specified at THETA = C and the solution is specified at
  !           THETA = d (see note below).
  !
  !      NOTE:  When NBDCND = 1, 2, or 4, do not use MBDCND = 5 or 6 (the
  !      former indicates that the solution is specified at R =  0; the
  !      latter indicates the solution is unspecified at R = 0).  Use
  !      instead MBDCND = 1 or 2.
  !
  !    BDC
  !      A one dimensional array of length M that specifies the boundary
  !      values of the solution at THETA = C.   When NBDCND = 1 or 2,
  !
  !               BDC(I) = U(R(I),C),              I=1,2,...,M.
  !
  !      When NBDCND = 3 or 4,
  !
  !               BDC(I) = (d/dTHETA)U(R(I),C),     I=1,2,...,M.
  !
  !      When NBDCND = 0, BDC is a dummy variable.
  !
  !    BDD
  !      A one-dimensional array of length M that specifies the boundary
  !      values of the solution at THETA = D.  When NBDCND = 1 or 4,
  !
  !               BDD(I) = U(R(I),D),              I=1,2,...,M.
  !
  !      When NBDCND = 2 or 3,
  !
  !               BDD(I) = (d/dTHETA)U(R(I),D),    I=1,2,...,M.
  !
  !      When NBDCND = 0, BDD is a dummy variable.
  !
  !    ELMBDA
  !      The constant LAMBDA in the Helmholtz equation.  If LAMBDA is
  !      greater than 0, a solution may not exist.  However, HSTPLR will
  !      attempt to find a solution.
  !
  !    F
  !      A two-dimensional array that specifies the values of the right
  !      side of the Helmholtz equation.  For I=1,2,...,M and J=1,2,...,N
  !
  !               F(I,J) = F(R(I),THETA(J)) .
  !
  !      F must be dimensioned at least M X N.
  !
  !    IDIMF
  !      The row (or first) dimension of the array F as it appears in the
  !      program calling HSTPLR.  This parameter is used to specify the
  !      variable dimension of F.  IDIMF must be at least M.
  !
  !    W
  !      A one-dimensional array that must be provided by the user for
  !      work space.  W may require up to 13M + 4N + M*INT(log2(N))
  !      locations.  The actual number of locations used is computed by
  !      HSTPLR and is returned in the location W(1).
  !
  !
  !             * * * * * *   On Output   * * * * * *
  !
  !    F
  !      Contains the solution U(I,J) of the finite difference
  !      approximation for the grid point (R(I),THETA(J)) for
  !      I=1,2,...,M, J=1,2,...,N.
  !
  !    PERTRB
  !      If a combination of periodic, derivative, or unspecified
  !      boundary conditions is specified for a Poisson equation
  !      (LAMBDA = 0), a solution may not exist.  PERTRB is a con-
  !      stant, calculated and subtracted from F, which ensures
  !      that a solution exists.  HSTPLR then computes this
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
  !      =  7  A = 0 and MBDCND = 3 or 4
  !
  !      =  8  A .GT. 0 and MBDCND .GE. 5
  !
  !      =  9  MBDCND .GE. 5 and NBDCND .NE. 0 or 3
  !
  !      = 10  IDIMF .LT. M
  !
  !      = 11  LAMBDA .GT. 0
  !
  !      = 12  M .LE. 2
  !
  !      Since this is the only means of indicating a possibly
  !      incorrect call to HSTPLR, the user should test IERROR after
  !      the call.
  !
  !    W
  !      W(1) contains the required length of W.
  !
  !- Long Description:
  !
  !     * * * * * * *   Program Specifications    * * * * * * * * * * * *
  !
  !     Dimension of   BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N),
  !     Arguments      W(see ARGUMENT LIST)
  !
  !     Latest         June 1, 1977
  !     Revision
  !
  !     Subprograms    HSTPLR,POISTG,POSTG2,GENBUN,POISD2,POISN2,POISP2,
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
  !     History        Written by Roland Sweet at NCAR in February, 1977
  !
  !     Algorithm      This subroutine defines the finite-difference
  !                    equations, incorporates boundary data, adjusts the
  !                    right side when the system is singular and calls
  !                    either POISTG or GENBUN which solves the linear
  !                    system of equations.
  !
  !     Space          8265(decimal) = 20111(octal) LOCATIONS ON THE
  !     Required       NCAR Control Data 7600
  !
  !     Timing and        The execution time T on the NCAR Control Data
  !     Accuracy       7600 for subroutine HSTPLR is roughly proportional
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
  !                    The Solution Of Poisson's Equation With Neumann
  !                    Boundary Conditions On A Staggered Grid of
  !                    Arbitrary Size,' J. Comp. Phys. 20(1976),
  !                    pp. 171-182.
  !
  !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  !
  !***
  ! **References:**  U. Schumann and R. Sweet, A direct method for the
  !                 solution of Poisson's equation with Neumann boundary
  !                 conditions on a staggered grid of arbitrary size,
  !                 Journal of Computational Physics 20, (1976),
  !                 pp. 171-182.
  !***
  ! **Routines called:**  GENBUN, POISTG

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  
  REAL A, a1, a2, B, Bda, Bdb, Bdc, Bdd, C, D, deltar, deltht, &
    dlrsq, dlthsq, Elmbda, F, Pertrb, W
  INTEGER i, Idimf, ierr1, Ierror, isw, iwb, iwc, iwr, j, k, lp, &
    M, mb, Mbdcnd, N, Nbdcnd, np
  DIMENSION F(Idimf,*)
  DIMENSION Bda(*), Bdb(*), Bdc(*), Bdd(*), W(*)
  !* FIRST EXECUTABLE STATEMENT  HSTPLR
  Ierror = 0
  IF ( A<0. ) Ierror = 1
  IF ( A>=B ) Ierror = 2
  IF ( Mbdcnd<=0.OR.Mbdcnd>=7 ) Ierror = 3
  IF ( C>=D ) Ierror = 4
  IF ( N<=2 ) Ierror = 5
  IF ( Nbdcnd<0.OR.Nbdcnd>=5 ) Ierror = 6
  IF ( A==0..AND.(Mbdcnd==3.OR.Mbdcnd==4) ) Ierror = 7
  IF ( A>0..AND.Mbdcnd>=5 ) Ierror = 8
  IF ( Mbdcnd>=5.AND.Nbdcnd/=0.AND.Nbdcnd/=3 ) Ierror = 9
  IF ( Idimf<M ) Ierror = 10
  IF ( M<=2 ) Ierror = 12
  IF ( Ierror/=0 ) RETURN
  deltar = (B-A)/M
  dlrsq = deltar**2
  deltht = (D-C)/N
  dlthsq = deltht**2
  np = Nbdcnd + 1
  isw = 1
  mb = Mbdcnd
  IF ( A==0..AND.Mbdcnd==2 ) mb = 6
  !
  !     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
  !
  iwb = M
  iwc = iwb + M
  iwr = iwc + M
  DO i = 1, M
    j = iwr + i
    W(j) = A + (i-0.5)*deltar
    W(i) = (A+(i-1)*deltar)/dlrsq
    k = iwc + i
    W(k) = (A+i*deltar)/dlrsq
    k = iwb + i
    W(k) = (Elmbda-2./dlrsq)*W(j)
  ENDDO
  DO i = 1, M
    j = iwr + i
    a1 = W(j)
    DO j = 1, N
      F(i,j) = a1*F(i,j)
    ENDDO
  ENDDO
  !
  !     ENTER BOUNDARY DATA FOR R-BOUNDARIES.
  !
  SELECT CASE (mb)
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
  SELECT CASE (mb)
    CASE (2,3,6)
      a1 = deltar*W(iwr)
      W(iwc) = W(iwc) + W(iwr)
      DO j = 1, N
        F(M,j) = F(M,j) - a1*Bdb(j)
      ENDDO
    CASE DEFAULT
      a1 = 2.*W(iwr)
      W(iwc) = W(iwc) - W(iwr)
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
      CASE (1,2,4,5)
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
            Pertrb = Pertrb/(M*N*0.5*(A+B))
            DO i = 1, M
              j = iwr + i
              a1 = Pertrb*W(j)
              DO j = 1, N
                F(i,j) = F(i,j) - a1
              ENDDO
            ENDDO
            a2 = 0.
            DO j = 1, N
              a2 = a2 + F(1,j)
            ENDDO
            a2 = a2/W(iwr+1)
        END SELECT
    END SELECT
  ELSE
    Ierror = 11
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
  IF ( lp==0 ) THEN
    CALL GENBUN(lp,N,1,M,W,W(iwb+1),W(iwc+1),Idimf,F,ierr1,W(iwr+1))
  ELSE
    CALL POISTG(lp,N,1,M,W,W(iwb+1),W(iwc+1),Idimf,F,ierr1,W(iwr+1))
  ENDIF
  W(1) = W(iwr+1) + 3*M
  IF ( A==0..AND.Mbdcnd==2.AND.isw==2 ) THEN
    a1 = 0.
    DO j = 1, N
      a1 = a1 + F(1,j)
    ENDDO
    a1 = (a1-dlrsq*a2/16.)/N
    IF ( Nbdcnd==3 ) a1 = a1 + (Bdd(1)-Bdc(1))/(D-C)
    a1 = Bda(1) - a1
    DO i = 1, M
      DO j = 1, N
        F(i,j) = F(i,j) + a1
      ENDDO
    ENDDO
  ENDIF
END SUBROUTINE HSTPLR
