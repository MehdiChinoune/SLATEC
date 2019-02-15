!DECK HW3CRT
SUBROUTINE HW3CRT(Xs,Xf,L,Lbdcnd,Bdxs,Bdxf,Ys,Yf,M,Mbdcnd,Bdys,Bdyf,Zs,Zf,&
    N,Nbdcnd,Bdzs,Bdzf,Elmbda,Ldimf,Mdimf,F,Pertrb,Ierror,W)
  IMPLICIT NONE
  REAL Bdxf, Bdxs, Bdyf, Bdys, Bdzf, Bdzs, c1, c2, c3, dx, dy, &
    dz, Elmbda, F, Pertrb, s, s1, s2, twbydx, twbydy
  REAL twbydz, W, Xf, xlp, Xs, Yf, ylp, Ys, Zf, zlp, Zs
  INTEGER i, Ierror, ir, iwb, iwc, iww, j, k, L, Lbdcnd, Ldimf, &
    lp, lp1, lstart, lstop, lstpm1, lunk, M, Mbdcnd, Mdimf
  INTEGER mp, mp1, mstart, mstop, mstpm1, munk, N, Nbdcnd, np, &
    np1, nperod, nstart, nstop, nstpm1, nunk
  !***BEGIN PROLOGUE  HW3CRT
  !***PURPOSE  Solve the standard seven-point finite difference
  !            approximation to the Helmholtz equation in Cartesian
  !            coordinates.
  !***LIBRARY   SLATEC (FISHPACK)
  !***CATEGORY  I2B1A1A
  !***TYPE      SINGLE PRECISION (HW3CRT-S)
  !***KEYWORDS  CARTESIAN, ELLIPTIC, FISHPACK, HELMHOLTZ, PDE
  !***AUTHOR  Adams, J., (NCAR)
  !           Swarztrauber, P. N., (NCAR)
  !           Sweet, R., (NCAR)
  !***DESCRIPTION
  !
  !     Subroutine HW3CRT solves the standard seven-point finite
  !     difference approximation to the Helmholtz equation in Cartesian
  !     coordinates:
  !
  !         (d/dX)(dU/dX) + (d/dY)(dU/dY) + (d/dZ)(dU/dZ)
  !
  !                    + LAMBDA*U = F(X,Y,Z) .
  !
  !    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  !
  !
  !    * * * * * * * *    Parameter Description     * * * * * * * * * *
  !
  !
  !            * * * * * *   On Input    * * * * * *
  !
  !     XS,XF
  !        The range of X, i.e. XS .LE. X .LE. XF .
  !        XS must be less than XF.
  !
  !     L
  !        The number of panels into which the interval (XS,XF) is
  !        subdivided.  Hence, there will be L+1 grid points in the
  !        X-direction given by X(I) = XS+(I-1)DX for I=1,2,...,L+1,
  !        where DX = (XF-XS)/L is the panel width.  L must be at
  !        least 5 .
  !
  !     LBDCND
  !        Indicates the type of boundary conditions at X = XS and X = XF.
  !
  !        = 0  If the solution is periodic in X, i.e.
  !             U(L+I,J,K) = U(I,J,K).
  !        = 1  If the solution is specified at X = XS and X = XF.
  !        = 2  If the solution is specified at X = XS and the derivative
  !             of the solution with respect to X is specified at X = XF.
  !        = 3  If the derivative of the solution with respect to X is
  !             specified at X = XS and X = XF.
  !        = 4  If the derivative of the solution with respect to X is
  !             specified at X = XS and the solution is specified at X=XF.
  !
  !     BDXS
  !        A two-dimensional array that specifies the values of the
  !        derivative of the solution with respect to X at X = XS.
  !        when LBDCND = 3 or 4,
  !
  !             BDXS(J,K) = (d/dX)U(XS,Y(J),Z(K)), J=1,2,...,M+1,
  !                                                K=1,2,...,N+1.
  !
  !        When LBDCND has any other value, BDXS is a dummy variable.
  !        BDXS must be dimensioned at least (M+1)*(N+1).
  !
  !     BDXF
  !        A two-dimensional array that specifies the values of the
  !        derivative of the solution with respect to X at X = XF.
  !        When LBDCND = 2 or 3,
  !
  !             BDXF(J,K) = (d/dX)U(XF,Y(J),Z(K)), J=1,2,...,M+1,
  !                                                K=1,2,...,N+1.
  !
  !        When LBDCND has any other value, BDXF is a dummy variable.
  !        BDXF must be dimensioned at least (M+1)*(N+1).
  !
  !     YS,YF
  !        The range of Y, i.e. YS .LE. Y .LE. YF.
  !        YS must be less than YF.
  !
  !     M
  !        The number of panels into which the interval (YS,YF) is
  !        subdivided.  Hence, there will be M+1 grid points in the
  !        Y-direction given by Y(J) = YS+(J-1)DY for J=1,2,...,M+1,
  !        where DY = (YF-YS)/M is the panel width.  M must be at
  !        least 5 .
  !
  !     MBDCND
  !        Indicates the type of boundary conditions at Y = YS and Y = YF.
  !
  !        = 0  If the solution is periodic in Y, i.e.
  !             U(I,M+J,K) = U(I,J,K).
  !        = 1  If the solution is specified at Y = YS and Y = YF.
  !        = 2  If the solution is specified at Y = YS and the derivative
  !             of the solution with respect to Y is specified at Y = YF.
  !        = 3  If the derivative of the solution with respect to Y is
  !             specified at Y = YS and Y = YF.
  !        = 4  If the derivative of the solution with respect to Y is
  !             specified at Y = YS and the solution is specified at Y=YF.
  !
  !     BDYS
  !        A two-dimensional array that specifies the values of the
  !        derivative of the solution with respect to Y at Y = YS.
  !        When MBDCND = 3 or 4,
  !
  !             BDYS(I,K) = (d/dY)U(X(I),YS,Z(K)), I=1,2,...,L+1,
  !                                                K=1,2,...,N+1.
  !
  !        When MBDCND has any other value, BDYS is a dummy variable.
  !        BDYS must be dimensioned at least (L+1)*(N+1).
  !
  !     BDYF
  !        A two-dimensional array that specifies the values of the
  !        derivative of the solution with respect to Y at Y = YF.
  !        When MBDCND = 2 or 3,
  !
  !             BDYF(I,K) = (d/dY)U(X(I),YF,Z(K)), I=1,2,...,L+1,
  !                                                K=1,2,...,N+1.
  !
  !        When MBDCND has any other value, BDYF is a dummy variable.
  !        BDYF must be dimensioned at least (L+1)*(N+1).
  !
  !     ZS,ZF
  !        The range of Z, i.e. ZS .LE. Z .LE. ZF.
  !        ZS must be less than ZF.
  !
  !     N
  !        The number of panels into which the interval (ZS,ZF) is
  !        subdivided.  Hence, there will be N+1 grid points in the
  !        Z-direction given by Z(K) = ZS+(K-1)DZ for K=1,2,...,N+1,
  !        where DZ = (ZF-ZS)/N is the panel width.  N must be at least 5.
  !
  !     NBDCND
  !        Indicates the type of boundary conditions at Z = ZS and Z = ZF.
  !
  !        = 0  If the solution is periodic in Z, i.e.
  !             U(I,J,N+K) = U(I,J,K).
  !        = 1  If the solution is specified at Z = ZS and Z = ZF.
  !        = 2  If the solution is specified at Z = ZS and the derivative
  !             of the solution with respect to Z is specified at Z = ZF.
  !        = 3  If the derivative of the solution with respect to Z is
  !             specified at Z = ZS and Z = ZF.
  !        = 4  If the derivative of the solution with respect to Z is
  !             specified at Z = ZS and the solution is specified at Z=ZF.
  !
  !     BDZS
  !        A two-dimensional array that specifies the values of the
  !        derivative of the solution with respect to Z at Z = ZS.
  !        When NBDCND = 3 or 4,
  !
  !             BDZS(I,J) = (d/dZ)U(X(I),Y(J),ZS), I=1,2,...,L+1,
  !                                                J=1,2,...,M+1.
  !
  !        When NBDCND has any other value, BDZS is a dummy variable.
  !        BDZS must be dimensioned at least (L+1)*(M+1).
  !
  !     BDZF
  !        A two-dimensional array that specifies the values of the
  !        derivative of the solution with respect to Z at Z = ZF.
  !        When NBDCND = 2 or 3,
  !
  !             BDZF(I,J) = (d/dZ)U(X(I),Y(J),ZF), I=1,2,...,L+1,
  !                                                J=1,2,...,M+1.
  !
  !        When NBDCND has any other value, BDZF is a dummy variable.
  !        BDZF must be dimensioned at least (L+1)*(M+1).
  !
  !     ELMBDA
  !        The constant LAMBDA in the Helmholtz equation. If
  !        LAMBDA .GT. 0, a solution may not exist.  However, HW3CRT will
  !        attempt to find a solution.
  !
  !     F
  !        A three-dimensional array that specifies the values of the
  !        right side of the Helmholtz equation and boundary values (if
  !        any).  For I=2,3,...,L, J=2,3,...,M, and K=2,3,...,N
  !
  !                   F(I,J,K) = F(X(I),Y(J),Z(K)).
  !
  !        On the boundaries F is defined by
  !
  !        LBDCND      F(1,J,K)         F(L+1,J,K)
  !        ------   ---------------   ---------------
  !
  !          0      F(XS,Y(J),Z(K))   F(XS,Y(J),Z(K))
  !          1      U(XS,Y(J),Z(K))   U(XF,Y(J),Z(K))
  !          2      U(XS,Y(J),Z(K))   F(XF,Y(J),Z(K))   J=1,2,...,M+1
  !          3      F(XS,Y(J),Z(K))   F(XF,Y(J),Z(K))   K=1,2,...,N+1
  !          4      F(XS,Y(J),Z(K))   U(XF,Y(J),Z(K))
  !
  !        MBDCND      F(I,1,K)         F(I,M+1,K)
  !        ------   ---------------   ---------------
  !
  !          0      F(X(I),YS,Z(K))   F(X(I),YS,Z(K))
  !          1      U(X(I),YS,Z(K))   U(X(I),YF,Z(K))
  !          2      U(X(I),YS,Z(K))   F(X(I),YF,Z(K))   I=1,2,...,L+1
  !          3      F(X(I),YS,Z(K))   F(X(I),YF,Z(K))   K=1,2,...,N+1
  !          4      F(X(I),YS,Z(K))   U(X(I),YF,Z(K))
  !
  !        NBDCND      F(I,J,1)         F(I,J,N+1)
  !        ------   ---------------   ---------------
  !
  !          0      F(X(I),Y(J),ZS)   F(X(I),Y(J),ZS)
  !          1      U(X(I),Y(J),ZS)   U(X(I),Y(J),ZF)
  !          2      U(X(I),Y(J),ZS)   F(X(I),Y(J),ZF)   I=1,2,...,L+1
  !          3      F(X(I),Y(J),ZS)   F(X(I),Y(J),ZF)   J=1,2,...,M+1
  !          4      F(X(I),Y(J),ZS)   U(X(I),Y(J),ZF)
  !
  !        F must be dimensioned at least (L+1)*(M+1)*(N+1).
  !
  !        NOTE:
  !
  !        If the table calls for both the solution U and the right side F
  !        on a boundary, then the solution must be specified.
  !
  !     LDIMF
  !        The row (or first) dimension of the arrays F,BDYS,BDYF,BDZS,
  !        and BDZF as it appears in the program calling HW3CRT. this
  !        parameter is used to specify the variable dimension of these
  !        arrays.  LDIMF must be at least L+1.
  !
  !     MDIMF
  !        The column (or second) dimension of the array F and the row (or
  !        first) dimension of the arrays BDXS and BDXF as it appears in
  !        the program calling HW3CRT.  This parameter is used to specify
  !        the variable dimension of these arrays.
  !        MDIMF must be at least M+1.
  !
  !     W
  !        A one-dimensional array that must be provided by the user for
  !        work space.  The length of W must be at least 30 + L + M + 5*N
  !        + MAX(L,M,N) + 7*(INT((L+1)/2) + INT((M+1)/2))
  !
  !
  !            * * * * * *   On Output   * * * * * *
  !
  !     F
  !        Contains the solution U(I,J,K) of the finite difference
  !        approximation for the grid point (X(I),Y(J),Z(K)) for
  !        I=1,2,...,L+1, J=1,2,...,M+1, and K=1,2,...,N+1.
  !
  !     PERTRB
  !        If a combination of periodic or derivative boundary conditions
  !        is specified for a Poisson equation (LAMBDA = 0), a solution
  !        may not exist.  PERTRB is a constant, calculated and subtracted
  !        from F, which ensures that a solution exists.  PWSCRT then
  !        computes this solution, which is a least squares solution to
  !        the original approximation.  This solution is not unique and is
  !        unnormalized.  The value of PERTRB should be small compared to
  !        the right side F.  Otherwise, a solution is obtained to an
  !        essentially different problem.  This comparison should always
  !        be made to insure that a meaningful solution has been obtained.
  !
  !     IERROR
  !        An error flag that indicates invalid input parameters.  Except
  !        for numbers 0 and 12, a solution is not attempted.
  !
  !        =  0  No error
  !        =  1  XS .GE. XF
  !        =  2  L .LT. 5
  !        =  3  LBDCND .LT. 0 .OR. LBDCND .GT. 4
  !        =  4  YS .GE. YF
  !        =  5  M .LT. 5
  !        =  6  MBDCND .LT. 0 .OR. MBDCND .GT. 4
  !        =  7  ZS .GE. ZF
  !        =  8  N .LT. 5
  !        =  9  NBDCND .LT. 0 .OR. NBDCND .GT. 4
  !        = 10  LDIMF .LT. L+1
  !        = 11  MDIMF .LT. M+1
  !        = 12  LAMBDA .GT. 0
  !
  !        Since this is the only means of indicating a possibly incorrect
  !        call to HW3CRT, the user should test IERROR after the call.
  !
  ! *Long Description:
  !
  !    * * * * * * *   Program Specifications    * * * * * * * * * * * *
  !
  !     Dimension of   BDXS(MDIMF,N+1),BDXF(MDIMF,N+1),BDYS(LDIMF,N+1),
  !     Arguments      BDYF(LDIMF,N+1),BDZS(LDIMF,M+1),BDZF(LDIMF,M+1),
  !                    F(LDIMF,MDIMF,N+1),W(see argument list)
  !
  !     Latest         December 1, 1978
  !     Revision
  !
  !     Subprograms    HW3CRT,POIS3D,POS3D1,TRIDQ,RFFTI,RFFTF,RFFTF1,
  !     Required       RFFTB,RFFTB1,COSTI,COST,SINTI,SINT,COSQI,COSQF,
  !                    COSQF1,COSQB,COSQB1,SINQI,SINQF,SINQB,CFFTI,
  !                    CFFTI1,CFFTB,CFFTB1,PASSB2,PASSB3,PASSB4,PASSB,
  !                    CFFTF,CFFTF1,PASSF1,PASSF2,PASSF3,PASSF4,PASSF,
  !                    PIMACH
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
  !     History        Written by Roland Sweet at NCAR in July 1977
  !
  !     Algorithm      This subroutine defines the finite difference
  !                    equations, incorporates boundary data, and
  !                    adjusts the right side of singular systems and
  !                    then calls POIS3D to solve the system.
  !
  !     Space          7862(decimal) = 17300(octal) locations on the
  !     Required       NCAR Control Data 7600
  !
  !     Timing and        The execution time T on the NCAR Control Data
  !     Accuracy       7600 for subroutine HW3CRT is roughly proportional
  !                    to L*M*N*(log2(L)+log2(M)+5), but also depends on
  !                    input parameters LBDCND and MBDCND.  Some typical
  !                    values are listed in the table below.
  !                       The solution process employed results in a loss
  !                    of no more than three significant digits for L,M
  !                    and N as large as 32.  More detailed information
  !                    about accuracy can be found in the documentation
  !                    for subroutine POIS3D which is the routine that
  !                    actually solves the finite difference equations.
  !
  !
  !                       L(=M=N)     LBDCND(=MBDCND=NBDCND)      T(MSECS)
  !                       -------     ----------------------      --------
  !
  !                         16                  0                    300
  !                         16                  1                    302
  !                         16                  3                    348
  !                         32                  0                   1925
  !                         32                  1                   1929
  !                         32                  3                   2109
  !
  !     Portability    American National Standards Institute FORTRAN.
  !                    The machine dependent constant PI is defined in
  !                    function PIMACH.
  !
  !     Required       COS,SIN,ATAN
  !     Resident
  !     Routines
  !
  !     Reference      NONE
  !
  !    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  POIS3D
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  HW3CRT
  !
  !
  DIMENSION Bdxs(Mdimf,*), Bdxf(Mdimf,*), Bdys(Ldimf,*), Bdyf(Ldimf,*), &
    Bdzs(Ldimf,*), Bdzf(Ldimf,*), F(Ldimf,Mdimf,*), W(*)
  !***FIRST EXECUTABLE STATEMENT  HW3CRT
  Ierror = 0
  IF ( Xf<=Xs ) Ierror = 1
  IF ( L<5 ) Ierror = 2
  IF ( Lbdcnd<0.OR.Lbdcnd>4 ) Ierror = 3
  IF ( Yf<=Ys ) Ierror = 4
  IF ( M<5 ) Ierror = 5
  IF ( Mbdcnd<0.OR.Mbdcnd>4 ) Ierror = 6
  IF ( Zf<=Zs ) Ierror = 7
  IF ( N<5 ) Ierror = 8
  IF ( Nbdcnd<0.OR.Nbdcnd>4 ) Ierror = 9
  IF ( Ldimf<L+1 ) Ierror = 10
  IF ( Mdimf<M+1 ) Ierror = 11
  IF ( Ierror==0 ) THEN
    dy = (Yf-Ys)/M
    twbydy = 2./dy
    c2 = 1./(dy**2)
    mstart = 1
    mstop = M
    mp1 = M + 1
    mp = Mbdcnd + 1
    SELECT CASE (mp)
      CASE (1)
        GOTO 50
      CASE (4,5)
      CASE DEFAULT
        mstart = 2
    END SELECT
    SELECT CASE (mp)
      CASE (1,2,5)
      CASE DEFAULT
        mstop = mp1
    END SELECT
    50     munk = mstop - mstart + 1
    dz = (Zf-Zs)/N
    twbydz = 2./dz
    np = Nbdcnd + 1
    c3 = 1./(dz**2)
    np1 = N + 1
    nstart = 1
    nstop = N
    SELECT CASE (np)
      CASE (1)
        GOTO 100
      CASE (4,5)
      CASE DEFAULT
        nstart = 2
    END SELECT
    SELECT CASE (np)
      CASE (1,2,5)
      CASE DEFAULT
        nstop = np1
    END SELECT
    100    nunk = nstop - nstart + 1
    lp1 = L + 1
    dx = (Xf-Xs)/L
    c1 = 1./(dx**2)
    twbydx = 2./dx
    lp = Lbdcnd + 1
    lstart = 1
    lstop = L
    !
    !     ENTER BOUNDARY DATA FOR X-BOUNDARIES.
    !
    SELECT CASE (lp)
      CASE (1)
        GOTO 150
      CASE (4,5)
        DO j = mstart, mstop
          DO k = nstart, nstop
            F(1,j,k) = F(1,j,k) + twbydx*Bdxs(j,k)
          ENDDO
        ENDDO
      CASE DEFAULT
        lstart = 2
        DO j = mstart, mstop
          DO k = nstart, nstop
            F(2,j,k) = F(2,j,k) - c1*F(1,j,k)
          ENDDO
        ENDDO
    END SELECT
    SELECT CASE (lp)
      CASE (1)
      CASE (3,4)
        lstop = lp1
        DO j = mstart, mstop
          DO k = nstart, nstop
            F(lp1,j,k) = F(lp1,j,k) - twbydx*Bdxf(j,k)
          ENDDO
        ENDDO
      CASE DEFAULT
        DO j = mstart, mstop
          DO k = nstart, nstop
            F(L,j,k) = F(L,j,k) - c1*F(lp1,j,k)
          ENDDO
        ENDDO
    END SELECT
    150    lunk = lstop - lstart + 1
    !
    !     ENTER BOUNDARY DATA FOR Y-BOUNDARIES.
    !
    SELECT CASE (mp)
      CASE (1)
        GOTO 200
      CASE (4,5)
        DO i = lstart, lstop
          DO k = nstart, nstop
            F(i,1,k) = F(i,1,k) + twbydy*Bdys(i,k)
          ENDDO
        ENDDO
      CASE DEFAULT
        DO i = lstart, lstop
          DO k = nstart, nstop
            F(i,2,k) = F(i,2,k) - c2*F(i,1,k)
          ENDDO
        ENDDO
    END SELECT
    SELECT CASE (mp)
      CASE (1)
      CASE (3,4)
        DO i = lstart, lstop
          DO k = nstart, nstop
            F(i,mp1,k) = F(i,mp1,k) - twbydy*Bdyf(i,k)
          ENDDO
        ENDDO
      CASE DEFAULT
        DO i = lstart, lstop
          DO k = nstart, nstop
            F(i,M,k) = F(i,M,k) - c2*F(i,mp1,k)
          ENDDO
        ENDDO
    END SELECT
    !
    !     ENTER BOUNDARY DATA FOR Z-BOUNDARIES.
    !
    200 CONTINUE
    SELECT CASE (np)
      CASE (1)
        GOTO 250
      CASE (4,5)
        DO i = lstart, lstop
          DO j = mstart, mstop
            F(i,j,1) = F(i,j,1) + twbydz*Bdzs(i,j)
          ENDDO
        ENDDO
      CASE DEFAULT
        DO i = lstart, lstop
          DO j = mstart, mstop
            F(i,j,2) = F(i,j,2) - c3*F(i,j,1)
          ENDDO
        ENDDO
    END SELECT
    SELECT CASE (np)
      CASE (1)
      CASE (3,4)
        DO i = lstart, lstop
          DO j = mstart, mstop
            F(i,j,np1) = F(i,j,np1) - twbydz*Bdzf(i,j)
          ENDDO
        ENDDO
      CASE DEFAULT
        DO i = lstart, lstop
          DO j = mstart, mstop
            F(i,j,N) = F(i,j,N) - c3*F(i,j,np1)
          ENDDO
        ENDDO
    END SELECT
    !
    !     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
    !
    250    iwb = nunk + 1
    iwc = iwb + nunk
    iww = iwc + nunk
    DO k = 1, nunk
      i = iwc + k - 1
      W(k) = c3
      W(i) = c3
      i = iwb + k - 1
      W(i) = -2.*c3 + Elmbda
    ENDDO
    SELECT CASE (np)
      CASE (1,2)
        GOTO 300
      CASE (3)
      CASE DEFAULT
        W(iwc) = 2.*c3
    END SELECT
    SELECT CASE (np)
      CASE (1,2,5)
      CASE DEFAULT
        W(iwb-1) = 2.*c3
    END SELECT
    300    Pertrb = 0.
    !
    !     FOR SINGULAR PROBLEMS ADJUST DATA TO INSURE A SOLUTION WILL EXIST.
    !
    SELECT CASE (lp)
      CASE (2,3,5)
      CASE DEFAULT
        SELECT CASE (mp)
          CASE (2,3,5)
          CASE DEFAULT
            SELECT CASE (np)
              CASE (2,3,5)
              CASE DEFAULT
                IF ( Elmbda<0 ) THEN
                ELSEIF ( Elmbda==0 ) THEN
                  mstpm1 = mstop - 1
                  lstpm1 = lstop - 1
                  nstpm1 = nstop - 1
                  xlp = (2+lp)/3
                  ylp = (2+mp)/3
                  zlp = (2+np)/3
                  s1 = 0.
                  DO k = 2, nstpm1
                    DO j = 2, mstpm1
                      DO i = 2, lstpm1
                        s1 = s1 + F(i,j,k)
                      ENDDO
                      s1 = s1 + (F(1,j,k)+F(lstop,j,k))/xlp
                    ENDDO
                    s2 = 0.
                    DO i = 2, lstpm1
                      s2 = s2 + F(i,1,k) + F(i,mstop,k)
                    ENDDO
                    s2 = (s2+(F(1,1,k)+F(1,mstop,k)+F(lstop,1,k)+F(lstop,mstop,k))&
                      /xlp)/ylp
                    s1 = s1 + s2
                  ENDDO
                  s = (F(1,1,1)+F(lstop,1,1)+F(1,1,nstop)+F(lstop,1,nstop)&
                    +F(1,mstop,1)+F(lstop,mstop,1)+F(1,mstop,nstop)&
                    +F(lstop,mstop,nstop))/(xlp*ylp)
                  DO j = 2, mstpm1
                    DO i = 2, lstpm1
                      s = s + F(i,j,1) + F(i,j,nstop)
                    ENDDO
                  ENDDO
                  s2 = 0.
                  DO i = 2, lstpm1
                    s2 = s2 + F(i,1,1) + F(i,1,nstop) + F(i,mstop,1)&
                      + F(i,mstop,nstop)
                  ENDDO
                  s = s2/ylp + s
                  s2 = 0.
                  DO j = 2, mstpm1
                    s2 = s2 + F(1,j,1) + F(1,j,nstop) + F(lstop,j,1)&
                      + F(lstop,j,nstop)
                  ENDDO
                  s = s2/xlp + s
                  Pertrb = (s/zlp+s1)/((lunk+1.-xlp)*(munk+1.-ylp)*(nunk+1.-zlp))
                  DO i = 1, lunk
                    DO j = 1, munk
                      DO k = 1, nunk
                        F(i,j,k) = F(i,j,k) - Pertrb
                      ENDDO
                    ENDDO
                  ENDDO
                ELSE
                  Ierror = 12
                ENDIF
            END SELECT
        END SELECT
    END SELECT
    nperod = 0
    IF ( Nbdcnd/=0 ) THEN
      nperod = 1
      W(1) = 0.
      W(iww-1) = 0.
    ENDIF
    CALL POIS3D(Lbdcnd,lunk,c1,Mbdcnd,munk,c2,nperod,nunk,W,W(iwb),W(iwc),&
      Ldimf,Mdimf,F(lstart,mstart,nstart),ir,W(iww))
    !
    !     FILL IN SIDES FOR PERIODIC BOUNDARY CONDITIONS.
    !
    IF ( lp==1 ) THEN
      IF ( mp==1 ) THEN
        DO k = nstart, nstop
          F(1,mp1,k) = F(1,1,k)
        ENDDO
        mstop = mp1
      ENDIF
      IF ( np==1 ) THEN
        DO j = mstart, mstop
          F(1,j,np1) = F(1,j,1)
        ENDDO
        nstop = np1
      ENDIF
      DO j = mstart, mstop
        DO k = nstart, nstop
          F(lp1,j,k) = F(1,j,k)
        ENDDO
      ENDDO
    ENDIF
    IF ( mp==1 ) THEN
      IF ( np==1 ) THEN
        DO i = lstart, lstop
          F(i,1,np1) = F(i,1,1)
        ENDDO
        nstop = np1
      ENDIF
      DO i = lstart, lstop
        DO k = nstart, nstop
          F(i,mp1,k) = F(i,1,k)
        ENDDO
      ENDDO
    ENDIF
    IF ( np==1 ) THEN
      DO i = lstart, lstop
        DO j = mstart, mstop
          F(i,j,np1) = F(i,j,1)
        ENDDO
      ENDDO
    ENDIF
  ENDIF
END SUBROUTINE HW3CRT
