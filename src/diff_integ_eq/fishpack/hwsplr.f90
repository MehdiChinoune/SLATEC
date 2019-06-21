!** HWSPLR
SUBROUTINE HWSPLR(A,B,M,Mbdcnd,Bda,Bdb,C,D,N,Nbdcnd,Bdc,Bdd,Elmbda,F,&
    Idimf,Pertrb,Ierror,W)
  !> Solve a finite difference approximation to the Helmholtz
  !            equation in polar coordinates.
  !***
  ! **Library:**   SLATEC (FISHPACK)
  !***
  ! **Category:**  I2B1A1A
  !***
  ! **Type:**      SINGLE PRECISION (HWSPLR-S)
  !***
  ! **Keywords:**  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE, POLAR
  !***
  ! **Author:**  Adams, J., (NCAR)
  !           Swarztrauber, P. N., (NCAR)
  !           Sweet, R., (NCAR)
  !***
  ! **Description:**
  !
  !     Subroutine HWSPLR solves a finite difference approximation to the
  !     Helmholtz equation in polar coordinates:
  !
  !          (1/R)(d/dR)(R(dU/dR)) + (1/R**2)(d/dTHETA)(dU/dTHETA)
  !
  !                                + LAMBDA*U = F(R,THETA).
  !
  !
  !
  !
  !     * * * * * * * *    Parameter Description     * * * * * * * * * *
  !
  !             * * * * * *   On Input    * * * * * *
  !
  !     A,B
  !       The range of R, i.e., A <= R <= B.  A must be less than B
  !       and A must be non-negative.
  !
  !     M
  !       The number of panels into which the interval (A,B) is
  !       subdivided.  Hence, there will be M+1 grid points in the
  !       R-direction given by R(I) = A+(I-1)DR, for I = 1,2,...,M+1,
  !       where DR = (B-A)/M is the panel width. M must be greater than 3.
  !
  !     MBDCND
  !       Indicates the type of boundary condition at R = A and R = B.
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
  !            BDA(J) = (d/dR)U(A,THETA(J)), J = 1,2,...,N+1  .
  !
  !       When MBDCND has any other value, BDA is a dummy variable.
  !
  !     BDB
  !       A one-dimensional array of length N+1 that specifies the values
  !       of the derivative of the solution with respect to R at R = B.
  !       When MBDCND = 2,3, or 6,
  !
  !            BDB(J) = (d/dR)U(B,THETA(J)), J = 1,2,...,N+1  .
  !
  !       When MBDCND has any other value, BDB is a dummy variable.
  !
  !     C,D
  !       The range of THETA, i.e., C <= THETA <= D.  C must be less
  !       than D.
  !
  !     N
  !       The number of panels into which the interval (C,D) is
  !       subdivided.  Hence, there will be N+1 grid points in the
  !       THETA-direction given by THETA(J) = C+(J-1)DTHETA for
  !       J = 1,2,...,N+1, where DTHETA = (D-C)/N is the panel width.  N
  !       must be greater than 3.
  !
  !     NBDCND
  !       Indicates the type of boundary conditions at THETA = C and
  !       at THETA = D.
  !
  !       = 0  If the solution is periodic in THETA, i.e.,
  !            U(I,J) = U(I,N+J).
  !       = 1  If the solution is specified at THETA = C and THETA = D
  !            (see note below).
  !       = 2  If the solution is specified at THETA = C and the
  !            derivative of the solution with respect to THETA is
  !            specified at THETA = D (see note below).
  !       = 4  If the derivative of the solution with respect to THETA is
  !            specified at THETA = C and the solution is specified at
  !            THETA = D (see note below).
  !
  !       NOTE:  When NBDCND = 1,2, or 4, do not use MBDCND = 5 or 6
  !              (the former indicates that the solution is specified at
  !              R = 0, the latter indicates the solution is unspecified
  !              at R = 0).  Use instead MBDCND = 1 or 2  .
  !
  !     BDC
  !       A one-dimensional array of length M+1 that specifies the values
  !       of the derivative of the solution with respect to THETA at
  !       THETA = C.  When NBDCND = 3 or 4,
  !
  !            BDC(I) = (d/dTHETA)U(R(I),C), I = 1,2,...,M+1  .
  !
  !       When NBDCND has any other value, BDC is a dummy variable.
  !
  !     BDD
  !       A one-dimensional array of length M+1 that specifies the values
  !       of the derivative of the solution with respect to THETA at
  !       THETA = D.  When NBDCND = 2 or 3,
  !
  !            BDD(I) = (d/dTHETA)U(R(I),D), I = 1,2,...,M+1  .
  !
  !       When NBDCND has any other value, BDD is a dummy variable.
  !
  !     ELMBDA
  !       The constant LAMBDA in the Helmholtz equation.  If
  !       LAMBDA < 0, a solution may not exist.  However, HWSPLR will
  !       attempt to find a solution.
  !
  !     F
  !       A two-dimensional array that specifies the values of the right
  !       side of the Helmholtz equation and boundary values (if any).
  !       For I = 2,3,...,M and J = 2,3,...,N
  !
  !            F(I,J) = F(R(I),THETA(J)).
  !
  !       On the boundaries F is defined by
  !
  !            MBDCND   F(1,J)            F(M+1,J)
  !            ------   -------------     -------------
  !
  !              1      U(A,THETA(J))     U(B,THETA(J))
  !              2      U(A,THETA(J))     F(B,THETA(J))
  !              3      F(A,THETA(J))     F(B,THETA(J))
  !              4      F(A,THETA(J))     U(B,THETA(J))   J = 1,2,...,N+1
  !              5      F(0,0)            U(B,THETA(J))
  !              6      F(0,0)            F(B,THETA(J))
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
  !
  !     IDIMF
  !       The row (or first) dimension of the array F as it appears in the
  !       program calling HWSPLR.  This parameter is used to specify the
  !       variable dimension of F.  IDIMF must be at least M+1  .
  !
  !     W
  !       A one-dimensional array that must be provided by the user for
  !       work space.  W may require up to 4*(N+1) +
  !       (13 + INT(log2(N+1)))*(M+1) locations.  The actual number of
  !       locations used is computed by HWSPLR and is returned in location
  !       W(1).
  !
  !
  !             * * * * * *   On Output     * * * * * *
  !
  !     F
  !       Contains the solution U(I,J) of the finite difference
  !       approximation for the grid point (R(I),THETA(J)),
  !       I = 1,2,...,M+1, J = 1,2,...,N+1  .
  !
  !     PERTRB
  !       If a combination of periodic, derivative, or unspecified
  !       boundary conditions is specified for a Poisson equation
  !       (LAMBDA = 0), a solution may not exist.  PERTRB is a constant,
  !       calculated and subtracted from F, which ensures that a solution
  !       exists.  HWSPLR then computes this solution, which is a least
  !       squares solution to the original approximation.  This solution
  !       plus any constant is also a solution.  Hence, the solution is
  !       not unique.  PERTRB should be small compared to the right side.
  !       Otherwise, a solution is obtained to an essentially different
  !       problem.  This comparison should always be made to insure that a
  !       meaningful solution has been obtained.
  !
  !     IERROR
  !       An error flag that indicates invalid input parameters.  Except
  !       for numbers 0 and 11, a solution is not attempted.
  !
  !       =  0  No error.
  !       =  1  A < 0  .
  !       =  2  A >= B.
  !       =  3  MBDCND < 1 or MBDCND > 6  .
  !       =  4  C >= D.
  !       =  5  N <= 3
  !       =  6  NBDCND < 0 or > 4  .
  !       =  7  A = 0, MBDCND = 3 or 4  .
  !       =  8  A > 0, MBDCND >= 5  .
  !       =  9  MBDCND >= 5, NBDCND /= 0 and NBDCND /= 3  .
  !       = 10  IDIMF < M+1  .
  !       = 11  LAMBDA > 0  .
  !       = 12  M <= 3
  !
  !       Since this is the only means of indicating a possibly incorrect
  !       call to HWSPLR, the user should test IERROR after the call.
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
  !     Subprograms    HWSPLR,GENBUN,POISD2,POISN2,POISP2,COSGEN,MERGE,
  !     Required       TRIX,TRI3,PIMACH
  !
  !     Special        None
  !     Conditions
  !
  !     Common         NONE
  !     Blocks
  !
  !     I/O
  !
  !     Precision      Single
  !
  !     Specialist     Roland Sweet
  !
  !     Language       FORTRAN
  !
  !     History        Standardized April 1, 1973
  !                    Revised January 1, 1976
  !
  !     Algorithm      The routine defines the finite difference
  !                    equations, incorporates boundary data, and adjusts
  !                    the right side of singular systems and then calls
  !                    GENBUN to solve the system.
  !
  !     Space          13430(octal) = 5912(decimal)  locations on the NCAR
  !     Required       Control Data 7600
  !
  !     Timing and        The execution time T on the NCAR Control Data
  !     Accuracy       7600 for subroutine HWSPLR is roughly proportional
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
  !                    Subprograms For The Solution Of Elliptic Equations'
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

  INTEGER :: Idimf, Ierror, M, Mbdcnd, N
  REAL(SP) :: A, B, C, D, Elmbda, Pertrb
  REAL(SP) :: Bda(N+1), Bdb(N+1), Bdc(M+1), Bdd(M+1), F(Idimf,N+1), W(:)
  INTEGER :: i, id2, id3, id4, id5, id6, ierr1, ij, ip, iwstor, j, k, l, lp, &
    mp1, mstart, mstop, munk, Nbdcnd, np, np1, nstart, nstop, nunk
  REAL(SP) :: a1, a2, deltar, deltht, dlrby2, dlrsq, dlthsq, r, s, s1, s2, ypole
  !* FIRST EXECUTABLE STATEMENT  HWSPLR
  Ierror = 0
  IF( A<0. ) Ierror = 1
  IF( A>=B ) Ierror = 2
  IF( Mbdcnd<=0 .OR. Mbdcnd>=7 ) Ierror = 3
  IF( C>=D ) Ierror = 4
  IF( N<=3 ) Ierror = 5
  IF( Nbdcnd<=-1 .OR. Nbdcnd>=5 ) Ierror = 6
  IF( A==0. .AND. (Mbdcnd==3 .OR. Mbdcnd==4) ) Ierror = 7
  IF( A>0. .AND. Mbdcnd>=5 ) Ierror = 8
  IF( Mbdcnd>=5 .AND. Nbdcnd/=0 .AND. Nbdcnd/=3 ) Ierror = 9
  IF( Idimf<M+1 ) Ierror = 10
  IF( M<=3 ) Ierror = 12
  IF( Ierror/=0 ) RETURN
  mp1 = M + 1
  deltar = (B-A)/M
  dlrby2 = deltar/2._SP
  dlrsq = deltar**2
  np1 = N + 1
  deltht = (D-C)/N
  dlthsq = deltht**2
  np = Nbdcnd + 1
  !
  !     DEFINE RANGE OF INDICES I AND J FOR UNKNOWNS U(I,J).
  !
  mstart = 2
  mstop = mp1
  SELECT CASE (Mbdcnd)
    CASE (2,6)
    CASE (3)
      mstart = 1
    CASE (4)
      mstart = 1
      mstop = M
    CASE (5)
      mstop = M
    CASE DEFAULT
      mstop = M
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
  a1 = 2._SP/dlrsq
  ij = 0
  IF( Mbdcnd==3 .OR. Mbdcnd==4 ) ij = 1
  DO i = 1, munk
    r = A + (i-ij)*deltar
    j = id5 + i
    W(j) = r
    j = id6 + i
    W(j) = 1._SP/r**2
    W(i) = (r-dlrby2)/(r*dlrsq)
    j = id3 + i
    W(j) = (r+dlrby2)/(r*dlrsq)
    j = id2 + i
    W(j) = -a1 + Elmbda
  END DO
  SELECT CASE (Mbdcnd)
    CASE (1,5)
    CASE (3)
      W(id2) = a1
      W(id3+1) = a1
    CASE (4)
      W(id3+1) = a1
    CASE DEFAULT
      W(id2) = a1
  END SELECT
  !
  !     ENTER BOUNDARY DATA FOR R-BOUNDARIES.
  !
  SELECT CASE (Mbdcnd)
    CASE (3,4)
      a1 = 2._SP*deltar*W(1)
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
      a1 = 2._SP*deltar*W(id4)
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
  !     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES.
  !
  a1 = 1._SP/dlthsq
  l = id5 - mstart + 1
  lp = id6 - mstart + 1
  SELECT CASE (np)
    CASE (1)
      GOTO 100
    CASE (4,5)
      a1 = 2._SP/deltht
      DO i = mstart, mstop
        j = i + lp
        F(i,1) = F(i,1) + a1*W(j)*Bdc(i)
      END DO
    CASE DEFAULT
      DO i = mstart, mstop
        j = i + lp
        F(i,2) = F(i,2) - a1*W(j)*F(i,1)
      END DO
  END SELECT
  a1 = 1._SP/dlthsq
  SELECT CASE (np)
    CASE (1)
    CASE (3,4)
      a1 = 2._SP/deltht
      DO i = mstart, mstop
        j = i + lp
        F(i,np1) = F(i,np1) - a1*W(j)*Bdd(i)
      END DO
    CASE DEFAULT
      DO i = mstart, mstop
        j = i + lp
        F(i,N) = F(i,N) - a1*W(j)*F(i,np1)
      END DO
  END SELECT
  !
  !     ADJUST RIGHT SIDE OF EQUATION FOR UNKNOWN AT POLE WHEN HAVE
  !     DERIVATIVE SPECIFIED BOUNDARY CONDITIONS.
  !
  100 CONTINUE
  IF( Mbdcnd>=5 .AND. Nbdcnd==3 ) F(1,1) = F(1,1) - (Bdd(2)-Bdc(2))&
    *4._SP/(N*deltht*dlrsq)
  !
  !     ADJUST RIGHT SIDE OF SINGULAR PROBLEMS TO INSURE EXISTENCE OF A
  !     SOLUTION.
  !
  Pertrb = 0._SP
  IF( Elmbda<0 ) THEN
  ELSEIF( Elmbda==0 ) THEN
    IF( Nbdcnd==0 .OR. Nbdcnd==3 ) THEN
      s2 = 0._SP
      SELECT CASE (Mbdcnd)
        CASE (1,2,4,5)
          GOTO 200
        CASE (6)
        CASE DEFAULT
          W(id5+1) = 0.5_SP*(W(id5+2)-dlrby2)
          s2 = 0.25_SP*deltar
      END SELECT
      a2 = 2._SP
      IF( Nbdcnd==0 ) a2 = 1._SP
      j = id5 + munk
      W(j) = 0.5_SP*(W(j-1)+dlrby2)
      s = 0._SP
      DO i = mstart, mstop
        s1 = 0._SP
        ij = nstart + 1
        k = nstop - 1
        DO j = ij, k
          s1 = s1 + F(i,j)
        END DO
        j = i + l
        s = s + (a2*s1+F(i,nstart)+F(i,nstop))*W(j)
      END DO
      s2 = M*A + deltar*((M-1)*(M+1)*.5_SP+.25_SP) + s2
      s1 = (2._SP+a2*(nunk-2))*s2
      IF( Mbdcnd/=3 ) THEN
        s2 = N*a2*deltar/8._SP
        s = s + F(1,1)*s2
        s1 = s1 + s2
      END IF
      Pertrb = s/s1
      DO i = mstart, mstop
        DO j = nstart, nstop
          F(i,j) = F(i,j) - Pertrb
        END DO
      END DO
    END IF
  ELSE
    Ierror = 11
  END IF
  !
  !     MULTIPLY I-TH EQUATION THROUGH BY (R(I)*DELTHT)**2.
  !
  200 CONTINUE
  DO i = mstart, mstop
    k = i - mstart + 1
    j = i + lp
    a1 = dlthsq/W(j)
    W(k) = a1*W(k)
    j = id2 + k
    W(j) = a1*W(j)
    j = id3 + k
    W(j) = a1*W(j)
    DO j = nstart, nstop
      F(i,j) = a1*F(i,j)
    END DO
  END DO
  W(1) = 0._SP
  W(id4) = 0._SP
  !
  !     CALL GENBUN TO SOLVE THE SYSTEM OF EQUATIONS.
  !
  CALL GENBUN(Nbdcnd,nunk,1,munk,W(1:munk),W(id2+1:id3),W(id3+1:id4),Idimf,&
    F(mstart,nstart),ierr1,W(id4+1:))
  iwstor = INT( W(id4+1) ) + 3*munk
  SELECT CASE (Mbdcnd)
    CASE (1,2,3,4)
      GOTO 400
    CASE (5)
    CASE DEFAULT
      !
      !     ADJUST THE SOLUTION AS NECESSARY FOR THE PROBLEMS WHERE A = 0.
      !
      IF( Elmbda==0. ) THEN
        ypole = 0._SP
        GOTO 300
      END IF
  END SELECT
  j = id5 + munk
  W(j) = W(id2)/W(id3)
  DO ip = 3, munk
    i = munk - ip + 2
    j = id5 + i
    lp = id2 + i
    k = id3 + i
    W(j) = W(i)/(W(lp)-W(k)*W(j+1))
  END DO
  W(id5+1) = -.5_SP*dlthsq/(W(id2+1)-W(id3+1)*W(id5+2))
  DO i = 2, munk
    j = id5 + i
    W(j) = -W(j)*W(j-1)
  END DO
  s = 0._SP
  DO j = nstart, nstop
    s = s + F(2,j)
  END DO
  a2 = nunk
  IF( Nbdcnd/=0 ) THEN
    s = s - 0.5_SP*(F(2,nstart)+F(2,nstop))
    a2 = a2 - 1._SP
  END IF
  ypole = (.25_SP*dlrsq*F(1,1)-s/a2)/(W(id5+1)-1._SP+Elmbda*dlrsq*.25_SP)
  DO i = mstart, mstop
    k = l + i
    DO j = nstart, nstop
      F(i,j) = F(i,j) + ypole*W(k)
    END DO
  END DO
  300 CONTINUE
  DO j = 1, np1
    F(1,j) = ypole
  END DO
  400 CONTINUE
  IF( Nbdcnd==0 ) THEN
    DO i = mstart, mstop
      F(i,np1) = F(i,1)
    END DO
  END IF
  W(1) = iwstor
END SUBROUTINE HWSPLR
