!** SEPX4
SUBROUTINE SEPX4(Iorder,A,B,M,Mbdcnd,Bda,Alpha,Bdb,Beta,C,D,N,Nbdcnd,Bdc,&
    Bdd,COFX,Grhs,Usol,Idmn,W,Pertrb,Ierror)
  !> Solve for either the second or fourth order finite difference approximation to
  !  the solution of a separable elliptic partial differential equation on a rectangle.
  !  Any combination of periodic or mixed boundary conditions is allowed.
  !***
  ! **Library:**   SLATEC (FISHPACK)
  !***
  ! **Category:**  I2B1A2
  !***
  ! **Type:**      SINGLE PRECISION (SEPX4-S)
  !***
  ! **Keywords:**  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE, SEPARABLE
  !***
  ! **Author:**  Adams, J., (NCAR)
  !           Swarztrauber, P. N., (NCAR)
  !           Sweet, R., (NCAR)
  !***
  ! **Description:**
  !
  ! Purpose                SEPX4 solves for either the second-order
  !                        finite difference approximation or a
  !                        fourth-order approximation  to the
  !                        solution of a separable elliptic equation
  !                             AF(X)*UXX+BF(X)*UX+CF(X)*U+UYY = G(X,Y)
  !
  !                        on a rectangle (X greater than or equal to A
  !                        and less than or equal to B; Y greater than
  !                        or equal to C and less than or equal to D).
  !                        Any combination of periodic or mixed boundary
  !                        conditions is allowed.
  !                        If boundary conditions in the X direction
  !                        are periodic (see MBDCND=0 below) then the
  !                        coefficients must satisfy
  !                        AF(X)=C1,BF(X)=0,CF(X)=C2 for all X.
  !                        Here C1,C2 are constants, C1>0.
  !
  !                        The possible boundary conditions are
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
  !                        In the Y-direction:
  !                         (0) Periodic, U(X,Y+D-C)=U(X,Y) for all X,Y
  !                         (1) U(X,C),U(X,D) are specified for all X
  !                         (2) U(X,C),dU(X,D)/dY are specified for all X
  !                         (3) dU(X,C)/DY,dU(X,D)/dY are specified for
  !                            all X
  !                        (4) dU(X,C)/DY,U(X,D) are specified for all X
  !
  ! Usage                  Call SEPX4(IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,
  !                                  BETA,C,D,N,NBDCND,BDC,BDD,COFX,
  !                                  GRHS,USOL,IDMN,W,PERTRB,IERROR)
  !
  ! Arguments
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
  !                          When MBDCND has any other value, BDA is a
  !                          dummy parameter.
  !
  ! On Input               ALPHA
  !                          The scalar multiplying the solution in case
  !                          of a mixed boundary condition AT X=A (see
  !                          argument BDA).  If MBDCND = 3,4 then ALPHA is
  !                          a dummy parameter.
  !
  !                        BDB
  !                          A one-dimensional array of length N+1 that
  !                          specifies the values of dU(B,Y)/dX+
  !                          BETA*U(B,Y) at X=B.  when MBDCND=2 or 3
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
  !                              i.e., dU(X,C)/dY and U(X,D)
  !                              are specified for all X
  !                          = 3 If the boundary conditions are mixed at
  !                              Y= C and Y=D i.e., dU(X,D)/DY
  !                              and dU(X,D)/dY are specified
  !                              for all X
  !                          = 4 If the boundary condition is mixed at Y=C
  !                              and the solution is specified at Y=D;
  !                              i.e. dU(X,C)/dY+GAMA*U(X,C) and U(X,D)
  !                              are specified for all X
  !
  !                        BDC
  !                          A one-dimensional array of length M+1 that
  !                          specifies the value dU(X,C)/DY
  !                          at Y=C.  When NBDCND=3 or 4
  !                            BDC(I) = dU(XI,C)/DY
  !                             I=1,2,...,M+1.
  !                          When NBDCND has any other value, BDC is a
  !                          dummy parameter.
  !
  !
  !                        BDD
  !                          A one-dimensional array of length M+1 that
  !                          specifies the value of dU(X,D)/DY
  !                          at Y=D.  When NBDCND=2 or 3
  !                            BDD(I)=dU(XI,D)/DY
  !                             I=1,2,...,M+1.
  !                          When NBDCND has any other value, BDD is a
  !                          dummy parameter.
  !
  !
  !                        COFX
  !                          A user-supplied subprogram with
  !                          parameters X, AFUN, BFUN, CFUN which
  !                          returns the values of the X-dependent
  !                          coefficients AF(X), BF(X), CF(X) in
  !                          the elliptic equation at X.
  !                          If boundary conditions in the X direction
  !                          are periodic then the coefficients
  !                          must satisfy AF(X)=C1,BF(X)=0,CF(X)=C2 for
  !                          all X.  Here C1>0 and C2 are constants.
  !
  !                          Note that COFX must be declared external
  !                          in the calling routine.
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
  !                          calling SEPX4.  This parameter is used to
  !                          specify the variable dimension of GRHS and
  !                          USOL.  IDMN must be at least 7 and greater
  !                          than or equal to M+1.
  !
  !                        W
  !                          A one-dimensional array that must be provided
  !                          by the user for work space.
  !                          10*N+(16+INT(log2(N)))*(M+1)+23 will suffice
  !                          as a length for W.  The actual length of
  !                          W in the calling routine must be set in W(1)
  !                          (see IERROR=11).
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
  !                          W(1) contains the exact minimal length (in
  !                          floating point) required for the work space
  !                          (see IERROR=11).
  !
  !                        PERTRB
  !                          If a combination of periodic or derivative
  !                          boundary conditions (i.e., ALPHA=BETA=0 if
  !                          MBDCND=3) is specified and if CF(X)=0 for all
  !                          X, then a solution to the discretized matrix
  !                          equation may not exist (reflecting the non-
  !                          uniqueness of solutions to the PDE).  PERTRB
  !                          is a constant calculated and subtracted from
  !                          the right hand side of the matrix equation
  !                          insuring the existence of a solution.
  !                          SEPX4 computes this solution which is a
  !                          weighted minimal least squares solution to
  !                          the original problem.  If singularity is
  !                          not detected PERTRB=0.0 is returned by
  !                          SEPX4.
  !
  !                        IERROR
  !                          An error flag that indicates invalid input
  !                          parameters or failure to find a solution
  !                          = 0  No error
  !                          = 1  If A greater than B or C greater than D
  !                          = 2  If MBDCND less than 0 or MBDCND greater
  !                               than 4
  !                          = 3  If NBDCND less than 0 or NBDCND greater
  !                               than 4
  !                          = 4  If attempt to find a solution fails.
  !                               (the linear system generated is not
  !                               diagonally dominant.)
  !                          = 5  If IDMN is too small (see discussion of
  !                               IDMN)
  !                          = 6  If M is too small or too large (see
  !                               discussion of M)
  !                          = 7  If N is too small (see discussion of N)
  !                          = 8  If IORDER is not 2 or 4
  !                          = 10 If AFUN is less than or equal to zero
  !                               for some interior mesh point XI
  !                          = 11 If the work space length input in W(1)
  !                               is less than the exact minimal work
  !                               space length required output in W(1).
  !                          = 12 If MBDCND=0 and AF(X)=CF(X)=constant
  !                               or BF(X)=0 for all X is not true.
  !
  !- Long Description:
  !
  ! Dimension of           BDA(N+1), BDB(N+1), BDC(M+1), BDD(M+1),
  ! Arguments              USOL(IDMN,N+1), GRHS(IDMN,N+1),
  !                        W (see argument list)
  !
  ! Latest Revision        October 1980
  !
  ! Special Conditions     NONE
  !
  ! Common Blocks          SPL4
  !
  ! I/O                    NONE
  !
  ! Precision              Single
  !
  ! Required Library       NONE
  ! Files
  !
  ! Specialist             John C. Adams, NCAR, Boulder, Colorado  80307
  !
  ! Language               FORTRAN
  !
  !
  ! Entry Points           SEPX4,SPELI4,CHKPR4,CHKSN4,ORTHO4,MINSO4,TRIS4,
  !                        DEFE4,DX4,DY4
  !
  ! History                SEPX4 was developed by modifying the ULIB
  !                        routine SEPELI during October 1978.
  !                        It should be used instead of SEPELI whenever
  !                        possible.  The increase in speed is at least
  !                        a factor of three.
  !
  ! Algorithm              SEPX4 automatically discretizes the separable
  !                        elliptic equation which is then solved by a
  !                        generalized cyclic reduction algorithm in the
  !                        subroutine POIS.  The fourth order solution
  !                        is obtained using the technique of
  !                        deferred corrections referenced below.
  !
  !
  ! References             Keller, H.B., 'Numerical Methods for Two-point
  !                          Boundary-value Problems', Blaisdel (1968),
  !                          Waltham, Mass.
  !
  !                        Swarztrauber, P., and R. Sweet (1975):
  !                          'Efficient FORTRAN Subprograms For The
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
  ! **Routines called:**  CHKPR4, SPELI4

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920122  Minor corrections and modifications to prologue.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTERFACE
    PURE SUBROUTINE COFX(X,A,B,C)
      IMPORT SP
      REAL(SP), INTENT(IN) :: X
      REAL(SP), INTENT(OUT) :: A, B, C
    END SUBROUTINE COFX
  END INTERFACE
  INTEGER, INTENT(IN) :: Idmn, Iorder, M, Mbdcnd, N, Nbdcnd
  INTEGER, INTENT(OUT) :: Ierror
  REAL(SP), INTENT(IN) :: A, Alpha, B, Beta, C, D
  REAL(SP), INTENT(OUT) :: Pertrb
  REAL(SP), INTENT(IN) :: Bda(N+1), Bdb(N+1), Bdc(M+1), Bdd(M+1)
  REAL(SP), INTENT(INOUT) :: Grhs(Idmn,N), W(:), Usol(Idmn,N+1)
  !
  INTEGER :: i1, i10, i11, i12, i13, i2, i3, i4, i5, i6, i7, i8, &
    i9, k, l, length, linput, log2n, loutpt
  !* FIRST EXECUTABLE STATEMENT  SEPX4
  CALL CHKPR4(Iorder,A,B,M,Mbdcnd,C,D,N,Nbdcnd,COFX,Idmn,Ierror)
  IF( Ierror/=0 ) RETURN
  !
  !     COMPUTE MINIMUM WORK SPACE AND CHECK WORK SPACE LENGTH INPUT
  !
  l = N + 1
  IF( Nbdcnd==0 ) l = N
  k = M + 1
  l = N + 1
  !     ESTIMATE LOG BASE 2 OF N
  log2n = INT(LOG(REAL(N+1))/LOG(2._SP)+0.5_SP)
  length = 4*(N+1) + (10+log2n)*(M+1)
  Ierror = 11
  linput = INT(W(1)+0.5_SP)
  loutpt = length + 6*(k+l) + 1
  W(1) = loutpt
  IF( loutpt>linput ) RETURN
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
  CALL SPELI4(Iorder,A,B,M,Mbdcnd,Bda,Alpha,Bdb,Beta,C,D,N,Nbdcnd,Bdc,Bdd,&
    COFX,W(i1:i2-1),W(i2:i3-1),W(i3:i4-1),W(i4:i5-1),W(i5:i6-1),W(i6:i7-1),W(i7:i8-1),&
    W(i8:i9-1),W(i9:i10-1),W(i10:i11-1),W(i11:i12-1),W(i12:i13-1),Grhs,Usol,Idmn,&
    W(i13:),Pertrb,Ierror)
  !
END SUBROUTINE SEPX4