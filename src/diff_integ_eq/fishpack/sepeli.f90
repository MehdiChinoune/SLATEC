!** SEPELI
SUBROUTINE SEPELI(Intl,Iorder,A,B,M,Mbdcnd,Bda,Alpha,Bdb,Beta,C,D,N,&
    Nbdcnd,Bdc,Gama,Bdd,Xnu,COFX,COFY,Grhs,Usol,Idmn,W,Pertrb,Ierror)
  IMPLICIT NONE
  !>
  !***
  !  Discretize and solve a second and, optionally, a fourth
  !            order finite difference approximation on a uniform grid to
  !            the general separable elliptic partial differential
  !            equation on a rectangle with any combination of periodic or
  !            mixed boundary conditions.
  !***
  ! **Library:**   SLATEC (FISHPACK)
  !***
  ! **Category:**  I2B1A2
  !***
  ! **Type:**      SINGLE PRECISION (SEPELI-S)
  !***
  ! **Keywords:**  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE, SEPARABLE
  !***
  ! **Author:**  Adams, J., (NCAR)
  !           Swarztrauber, P. N., (NCAR)
  !           Sweet, R., (NCAR)
  !***
  ! **Description:**
  !
  ! Dimension of           BDA(N+1), BDB(N+1), BDC(M+1), BDD(M+1),
  ! Arguments              USOL(IDMN,N+1), GRHS(IDMN,N+1),
  !                        W (see argument list)
  !
  ! Latest Revision        March 1977
  !
  ! Purpose                SEPELI solves for either the second-order
  !                        finite difference approximation or a
  !                        fourth-order approximation to a separable
  !                        elliptic equation.
  !
  !                                    2    2
  !                             AF(X)*d U/dX + BF(X)*dU/dX  + CF(X)*U +
  !                                    2    2
  !                             DF(Y)*d U/dY  + EF(Y)*dU/dY + FF(Y)*U
  !
  !                             = G(X,Y)
  !
  !                        on a rectangle (X greater than or equal to A
  !                        and less than or equal to B; Y greater than
  !                        or equal to C and less than or equal to D).
  !                        Any combination of periodic or mixed boundary
  !                        conditions is allowed.
  !
  ! Purpose                The possible boundary conditions are:
  !                        in the X-direction:
  !                         (0) Periodic, U(X+B-A,Y)=U(X,Y) for all Y,X
  !                         (1) U(A,Y), U(B,Y) are specified for all Y
  !                         (2) U(A,Y), dU(B,Y)/dX+BETA*U(B,Y) are
  !                             specified for all Y
  !                         (3) dU(A,Y)/dX+ALPHA*U(A,Y),dU(B,Y)/dX+
  !                             BETA*U(B,Y) are specified for all Y
  !                         (4) dU(A,Y)/dX+ALPHA*U(A,Y),U(B,Y) are
  !                             specified for all Y
  !
  !                        in the Y-direction:
  !                         (0) Periodic, U(X,Y+D-C)=U(X,Y) for all X,Y
  !                         (1) U(X,C),U(X,D) are specified for all X
  !                         (2) U(X,C),dU(X,D)/dY+XNU*U(X,D) are specified
  !                             for all X
  !                         (3) dU(X,C)/dY+GAMA*U(X,C),dU(X,D)/dY+
  !                             XNU*U(X,D) are specified for all X
  !                         (4) dU(X,C)/dY+GAMA*U(X,C),U(X,D) are
  !                             specified for all X
  !
  ! Arguments
  !
  ! On Input               INTL
  !                          = 0 On initial entry to SEPELI or if any of
  !                              the arguments C, D, N, NBDCND, COFY are
  !                              changed from a previous call
  !                          = 1 If C, D, N, NBDCND, COFY are unchanged
  !                              from the previous call.
  !
  !                        IORDER
  !                          = 2 If a second-order approximation is sought
  !                          = 4 If a fourth-order approximation is sought
  !
  !                        A,B
  !                          The range of the X-independent variable;
  !                          i.e., X is greater than or equal to A and
  !                          less than or equal to B.  A must be less than
  !                          B.
  !
  !                        M
  !                          The number of panels into which the interval
  !                          [A,B] is subdivided.  Hence, there will be
  !                          M+1 grid points in the X-direction given by
  !                          XI=A+(I-1)*DLX for I=1,2,...,M+1 where
  !                          DLX=(B-A)/M is the panel width.  M must be
  !                          less than IDMN and greater than 5.
  !
  !                        MBDCND
  !                          Indicates the type of boundary condition at
  !                          X=A and X=B
  !                          = 0 If the solution is periodic in X; i.e.,
  !                              U(X+B-A,Y)=U(X,Y) for all Y,X
  !                          = 1 If the solution is specified at X=A and
  !                              X=B; i.e., U(A,Y) and U(B,Y) are
  !                              specified for all Y
  !                          = 2 If the solution is specified at X=A and
  !                              the boundary condition is mixed at X=B;
  !                              i.e., U(A,Y) and dU(B,Y)/dX+BETA*U(B,Y)
  !                              are specified for all Y
  !                          = 3 If the boundary conditions at X=A and X=B
  !                              are mixed; i.e., dU(A,Y)/dX+ALPHA*U(A,Y)
  !                              and dU(B,Y)/dX+BETA*U(B,Y) are specified
  !                              for all Y
  !                          = 4 If the boundary condition at X=A is mixed
  !                              and the solution is specified at X=B;
  !                              i.e., dU(A,Y)/dX+ALPHA*U(A,Y) and U(B,Y)
  !                              are specified for all Y
  !
  !                        BDA
  !                          A one-dimensional array of length N+1 that
  !                          specifies the values of dU(A,Y)/dX+
  !                          ALPHA*U(A,Y) at X=A, when MBDCND=3 or 4.
  !                               BDA(J) = dU(A,YJ)/dX+ALPHA*U(A,YJ);
  !                               J=1,2,...,N+1
  !                          when MBDCND has any other value, BDA is a
  !                          dummy parameter.
  !
  ! On Input               ALPHA
  !                          The scalar multiplying the solution in case
  !                          of a mixed boundary condition at X=A (see
  !                          argument BDA).  If MBDCND = 3,4 then ALPHA is
  !                          a dummy parameter.
  !
  !                        BDB
  !                          A one-dimensional array of length N+1 that
  !                          specifies the values of dU(B,Y)/dX+
  !                          BETA*U(B,Y) at X=B.  When MBDCND=2 or 3
  !                               BDB(J) = dU(B,YJ)/dX+BETA*U(B,YJ);
  !                               J=1,2,...,N+1
  !                          When MBDCND has any other value, BDB is a
  !                          dummy parameter.
  !
  !                        BETA
  !                          The scalar multiplying the solution in case
  !                          of a mixed boundary condition at X=B (see
  !                          argument BDB).  If MBDCND=2,3 then BETA is a
  !                          dummy parameter.
  !
  !                        C,D
  !                          The range of the Y-independent variable;
  !                          i.e., Y is greater than or equal to C and
  !                          less than or equal to D.  C must be less than
  !                          D.
  !
  !                        N
  !                          The number of panels into which the interval
  !                          [C,D] is subdivided.  Hence, there will be
  !                          N+1 grid points in the Y-direction given by
  !                          YJ=C+(J-1)*DLY for J=1,2,...,N+1 where
  !                          DLY=(D-C)/N is the panel width.  In addition,
  !                          N must be greater than 4.
  !
  !                        NBDCND
  !                          Indicates the types of boundary conditions at
  !                          Y=C and Y=D
  !                          = 0 If the solution is periodic in Y; i.e.,
  !                              U(X,Y+D-C)=U(X,Y) for all X,Y
  !                          = 1 If the solution is specified at Y=C and
  !                              Y = D, i.e., U(X,C) and U(X,D) are
  !                              specified for all X
  !                          = 2 If the solution is specified at Y=C and
  !                              the boundary condition is mixed at Y=D;
  !                              i.e., U(X,C) and dU(X,D)/dY+XNU*U(X,D)
  !                              are specified for all X
  !                          = 3 If the boundary conditions are mixed at
  !                              Y=C and Y=D; i.e., dU(X,D)/dY+GAMA*U(X,C)
  !                              and dU(X,D)/dY+XNU*U(X,D) are specified
  !                              for all X
  !                          = 4 If the boundary condition is mixed at Y=C
  !                              and the solution is specified at Y=D;
  !                              i.e. dU(X,C)/dY+GAMA*U(X,C) and U(X,D)
  !                              are specified for all X
  !
  !                        BDC
  !                          A one-dimensional array of length M+1 that
  !                          specifies the value of dU(X,C)/dY+GAMA*U(X,C)
  !                          at Y=C.  When NBDCND=3 or 4
  !                             BDC(I) = dU(XI,C)/dY + GAMA*U(XI,C);
  !                             I=1,2,...,M+1.
  !                          When NBDCND has any other value, BDC is a
  !                          dummy parameter.
  !
  !                        GAMA
  !                          The scalar multiplying the solution in case
  !                          of a mixed boundary condition at Y=C (see
  !                          argument BDC).  If NBDCND=3,4 then GAMA is a
  !                          dummy parameter.
  !
  !                        BDD
  !                          A one-dimensional array of length M+1 that
  !                          specifies the value of dU(X,D)/dY +
  !                          XNU*U(X,D) at Y=C.  When NBDCND=2 or 3
  !                            BDD(I) = dU(XI,D)/dY + XNU*U(XI,D);
  !                            I=1,2,...,M+1.
  !                          When NBDCND has any other value, BDD is a
  !                          dummy parameter.
  !
  !                        XNU
  !                          The scalar multiplying the solution in case
  !                          of a mixed boundary condition at Y=D (see
  !                          argument BDD).  If NBDCND=2 or 3 then XNU is
  !                          a dummy parameter.
  !
  !                        COFX
  !                          A user-supplied subprogram with
  !                          parameters X, AFUN, BFUN, CFUN which
  !                          returns the values of the X-dependent
  !                          coefficients AF(X), BF(X), CF(X) in
  !                          the elliptic equation at X.
  !
  !                        COFY
  !                          A user-supplied subprogram with
  !                          parameters Y, DFUN, EFUN, FFUN which
  !                          returns the values of the Y-dependent
  !                          coefficients DF(Y), EF(Y), FF(Y) in
  !                          the elliptic equation at Y.
  !
  !                        NOTE:  COFX and COFY must be declared external
  !                        in the calling routine.  The values returned in
  !                        AFUN and DFUN must satisfy AFUN*DFUN greater
  !                        than 0 for A less than X less than B,
  !                        C less than Y less than D (see IERROR=10).
  !                        The coefficients provided may lead to a matrix
  !                        equation which is not diagonally dominant in
  !                        which case solution may fail (see IERROR=4).
  !
  !                        GRHS
  !                          A two-dimensional array that specifies the
  !                          values of the right-hand side of the elliptic
  !                          equation; i.e., GRHS(I,J)=G(XI,YI), for
  !                          I=2,...,M; J=2,...,N.  At the boundaries,
  !                          GRHS is defined by
  !
  !                          MBDCND   GRHS(1,J)   GRHS(M+1,J)
  !                          ------   ---------   -----------
  !                            0      G(A,YJ)     G(B,YJ)
  !                            1         *           *
  !                            2         *        G(B,YJ)  J=1,2,...,N+1
  !                            3      G(A,YJ)     G(B,YJ)
  !                            4      G(A,YJ)        *
  !
  !                          NBDCND   GRHS(I,1)   GRHS(I,N+1)
  !                          ------   ---------   -----------
  !                            0      G(XI,C)     G(XI,D)
  !                            1         *           *
  !                            2         *        G(XI,D)  I=1,2,...,M+1
  !                            3      G(XI,C)     G(XI,D)
  !                            4      G(XI,C)        *
  !
  !                          where * means these quantities are not used.
  !                          GRHS should be dimensioned IDMN by at least
  !                          N+1 in the calling routine.
  !
  !                        USOL
  !                          A two-dimensional array that specifies the
  !                          values of the solution along the boundaries.
  !                          At the boundaries, USOL is defined by
  !
  !                          MBDCND   USOL(1,J)   USOL(M+1,J)
  !                          ------   ---------   -----------
  !                            0         *           *
  !                            1      U(A,YJ)     U(B,YJ)
  !                            2      U(A,YJ)        *     J=1,2,...,N+1
  !                            3         *           *
  !                            4         *        U(B,YJ)
  !
  !                          NBDCND   USOL(I,1)   USOL(I,N+1)
  !                          ------   ---------   -----------
  !                            0         *           *
  !                            1      U(XI,C)     U(XI,D)
  !                            2      U(XI,C)        *     I=1,2,...,M+1
  !                            3         *           *
  !                            4         *        U(XI,D)
  !
  !                          where * means the quantities are not used in
  !                          the solution.
  !
  !                          If IORDER=2, the user may equivalence GRHS
  !                          and USOL to save space.  Note that in this
  !                          case the tables specifying the boundaries of
  !                          the GRHS and USOL arrays determine the
  !                          boundaries uniquely except at the corners.
  !                          If the tables call for both G(X,Y) and
  !                          U(X,Y) at a corner then the solution must be
  !                          chosen.  For example, if MBDCND=2 and
  !                          NBDCND=4, then U(A,C), U(A,D), U(B,D) must be
  !                          chosen at the corners in addition to G(B,C).
  !
  !                          If IORDER=4, then the two arrays, USOL and
  !                          GRHS, must be distinct.
  !
  !                          USOL should be dimensioned IDMN by at least
  !                          N+1 in the calling routine.
  !
  !                        IDMN
  !                          The row (or first) dimension of the arrays
  !                          GRHS and USOL as it appears in the program
  !                          calling SEPELI.  This parameter is used to
  !                          specify the variable dimension of GRHS and
  !                          USOL.  IDMN must be at least 7 and greater
  !                          than or equal to M+1.
  !
  !                        W
  !                          A one-dimensional array that must be provided
  !                          by the user for work space.  Let
  !                          K=INT(log2(N+1))+1 and set  L=2**(K+1).
  !                          then (K-2)*L+K+10*N+12*M+27 will suffice
  !                          as a length of W.  THE actual length of W in
  !                          the calling routine must be set in W(1) (see
  !                          IERROR=11).
  !
  ! On Output              USOL
  !                          Contains the approximate solution to the
  !                          elliptic equation.  USOL(I,J) is the
  !                          approximation to U(XI,YJ) for I=1,2...,M+1
  !                          and J=1,2,...,N+1.  The approximation has
  !                          error O(DLX**2+DLY**2) if called with
  !                          IORDER=2 and O(DLX**4+DLY**4) if called with
  !                          IORDER=4.
  !
  !                        W
  !                          Contains intermediate values that must not be
  !                          destroyed if SEPELI is called again with
  !                          INTL=1.  In addition W(1) contains the exact
  !                          minimal length (in floating point) required
  !                          for the work space (see IERROR=11).
  !
  !                        PERTRB
  !                          If a combination of periodic or derivative
  !                          boundary conditions (i.e., ALPHA=BETA=0 if
  !                          MBDCND=3; GAMA=XNU=0 if NBDCND=3) is
  !                          specified and if the coefficients of U(X,Y)
  !                          in the separable elliptic equation are zero
  !                          (i.e., CF(X)=0 for X greater than or equal to
  !                          A and less than or equal to B; FF(Y)=0 for
  !                          Y greater than or equal to C and less than
  !                          or equal to D) then a solution may not exist.
  !                          PERTRB is a constant calculated and
  !                          subtracted from the right-hand side of the
  !                          matrix equations generated by SEPELI which
  !                          insures that a solution exists.  SEPELI then
  !                          computes this solution which is a weighted
  !                          minimal least squares solution to the
  !                          original problem.
  !
  !                        IERROR
  !                          An error flag that indicates invalid input
  !                          parameters or failure to find a solution
  !                          = 0 No error
  !                          = 1 If A greater than B or C greater than D
  !                          = 2 If MBDCND less than 0 or MBDCND greater
  !                              than 4
  !                          = 3 If NBDCND less than 0 or NBDCND greater
  !                              than 4
  !                          = 4 If attempt to find a solution fails.
  !                              (the linear system generated is not
  !                              diagonally dominant.)
  !                          = 5 If IDMN is too small (see discussion of
  !                              IDMN)
  !                          = 6 If M is too small or too large (see
  !                              discussion of M)
  !                          = 7 If N is too small (see discussion of N)
  !                          = 8 If IORDER is not 2 or 4
  !                          = 9 If INTL is not 0 or 1
  !                          = 10 If AFUN*DFUN less than or equal to 0 for
  !                               some interior mesh point (XI,YJ)
  !                          = 11 If the work space length input in W(1)
  !                               is less than the exact minimal work
  !                               space length required output in W(1).
  !
  !                          NOTE (concerning IERROR=4):  for the
  !                          coefficients input through COFX, COFY, the
  !                          discretization may lead to a block
  !                          tridiagonal linear system which is not
  !                          diagonally dominant (for example, this
  !                          happens if CFUN=0 and BFUN/(2.*DLX) greater
  !                          than AFUN/DLX**2).  In this case solution may
  !                          fail.  This cannot happen in the limit as
  !                          DLX, DLY approach zero.  Hence, the condition
  !                          may be remedied by taking larger values for M
  !                          or N.
  !
  ! Entry Points           SEPELI, SPELIP, CHKPRM, CHKSNG, ORTHOG, MINSOL,
  !                        TRISP, DEFER, DX, DY, BLKTRI, BLKTR1, INDXB,
  !                        INDXA, INDXC, PROD, PRODP, CPROD, CPRODP,
  !                        PPADD, PSGF, BSRH, PPSGF, PPSPF, COMPB,
  !                        TRUN1, STOR1, TQLRAT
  !
  ! Special Conditions     NONE
  !
  ! Common Blocks          SPLP, CBLKT
  !
  ! I/O                    NONE
  !
  ! Precision              Single
  !
  ! Specialist             John C. Adams, NCAR, Boulder, Colorado  80307
  !
  ! Language               FORTRAN
  !
  ! History                Developed at NCAR during 1975-76.
  !
  ! Algorithm              SEPELI automatically discretizes the separable
  !                        elliptic equation which is then solved by a
  !                        generalized cyclic reduction algorithm in the
  !                        subroutine, BLKTRI.  The fourth-order solution
  !                        is obtained using 'Deferred Corrections' which
  !                        is described and referenced in sections,
  !                        references and method.
  !
  ! Space Required         14654 (octal) = 6572 (decimal)
  !
  ! Accuracy and Timing    The following computational results were
  !                        obtained by solving the sample problem at the
  !                        end of this write-up on the Control Data 7600.
  !                        The op count is proportional to M*N*log2(N).
  !                        In contrast to the other routines in this
  !                        chapter, accuracy is tested by computing and
  !                        tabulating second- and fourth-order
  !                        discretization errors.  Below is a table
  !                        containing computational results.  The times
  !                        given do not include initialization (i.e.,
  !                        times are for INTL=1).  Note that the
  !                        fourth-order accuracy is not realized until the
  !                        mesh is sufficiently refined.
  !
  !              Second-order    Fourth-order   Second-order  Fourth-order
  !    M    N   Execution Time  Execution Time    Error         Error
  !               (M SEC)         (M SEC)
  !     6    6         6              14          6.8E-1        1.2E0
  !    14   14        23              58          1.4E-1        1.8E-1
  !    30   30       100             247          3.2E-2        9.7E-3
  !    62   62       445           1,091          7.5E-3        3.0E-4
  !   126  126     2,002           4,772          1.8E-3        3.5E-6
  !
  ! Portability            There are no machine-dependent constants.
  !
  ! Required Resident      SQRT, ABS, LOG
  ! Routines
  !
  ! References             Keller, H.B., 'Numerical Methods for Two-point
  !                          Boundary-value Problems', Blaisdel (1968),
  !                          Waltham, Mass.
  !
  !                        Swarztrauber, P., and R. Sweet (1975):
  !                          'Efficient FORTRAN Subprograms for The
  !                          Solution of Elliptic Partial Differential
  !                          Equations'.  NCAR Technical Note
  !                          NCAR-TN/IA-109, pp. 135-137.
  !
  !***
  ! **References:**  H. B. Keller, Numerical Methods for Two-point
  !                 Boundary-value Problems, Blaisdel, Waltham, Mass.,
  !                 1968.
  !               P. N. Swarztrauber and R. Sweet, Efficient Fortran
  !                 subprograms for the solution of elliptic equations,
  !                 NCAR TN/IA-109, July 1975, 138 pp.
  !***
  ! **Routines called:**  CHKPRM, SPELIP

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER i1, i10, i11, i12, i13, i2, i3, i4, i5, i6, i7, i8, &
    i9, Idmn, Ierror, Intl, Iorder, k, l, length
  REAL A, Alpha, B, Bda(*), Bdb(*), Bdc(*), Bdd(*), Beta, C, D, &
    Gama, Grhs(Idmn,*), Pertrb, Usol(Idmn,*), W(*), Xnu
  INTEGER linput, ll, logb2n, loutpt, M, Mbdcnd, N, Nbdcnd
  EXTERNAL :: COFX, COFY
  !* FIRST EXECUTABLE STATEMENT  SEPELI
  CALL CHKPRM(Intl,Iorder,A,B,M,Mbdcnd,C,D,N,Nbdcnd,COFX,COFY,Idmn,Ierror)
  IF ( Ierror/=0 ) RETURN
  !
  !     COMPUTE MINIMUM WORK SPACE AND CHECK WORK SPACE LENGTH INPUT
  !
  l = N + 1
  IF ( Nbdcnd==0 ) l = N
  logb2n = INT(LOG(l+0.5)/LOG(2.0)) + 1
  ll = 2**(logb2n+1)
  k = M + 1
  l = N + 1
  length = (logb2n-2)*ll + logb2n + MAX(2*l,6*k) + 5
  IF ( Nbdcnd==0 ) length = length + 2*l
  Ierror = 11
  linput = INT(W(1)+0.5)
  loutpt = length + 6*(k+l) + 1
  W(1) = loutpt
  IF ( loutpt>linput ) RETURN
  Ierror = 0
  !
  !     SET WORK SPACE INDICES
  !
  i1 = length + 2
  i2 = i1 + l
  i3 = i2 + l
  i4 = i3 + l
  i5 = i4 + l
  i6 = i5 + l
  i7 = i6 + l
  i8 = i7 + k
  i9 = i8 + k
  i10 = i9 + k
  i11 = i10 + k
  i12 = i11 + k
  i13 = 2
  CALL SPELIP(Intl,Iorder,A,B,M,Mbdcnd,Bda,Alpha,Bdb,Beta,C,D,N,Nbdcnd,Bdc,&
    Gama,Bdd,Xnu,COFX,COFY,W(i1),W(i2),W(i3),W(i4),W(i5),W(i6),&
    W(i7),W(i8),W(i9),W(i10),W(i11),W(i12),Grhs,Usol,Idmn,W(i13),Pertrb,Ierror)
END SUBROUTINE SEPELI
