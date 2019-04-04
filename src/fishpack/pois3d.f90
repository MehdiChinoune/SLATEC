!** POIS3D
SUBROUTINE POIS3D(Lperod,L,C1,Mperod,M,C2,Nperod,N,A,B,C,Ldimf,Mdimf,F,Ierror,W)
  IMPLICIT NONE
  !>
  !***
  !  Solve a three-dimensional block tridiagonal linear system
  !            which arises from a finite difference approximation to a
  !            three-dimensional Poisson equation using the Fourier
  !            transform package FFTPAK written by Paul Swarztrauber.
  !***
  ! **Library:**   SLATEC (FISHPACK)
  !***
  ! **Category:**  I2B4B
  !***
  ! **Type:**      SINGLE PRECISION (POIS3D-S)
  !***
  ! **Keywords:**  ELLIPTIC PDE, FISHPACK, HELMHOLTZ, POISSON
  !***
  ! **Author:**  Adams, J., (NCAR)
  !           Swarztrauber, P. N., (NCAR)
  !           Sweet, R., (NCAR)
  !***
  ! **Description:**
  !
  !     Subroutine POIS3D solves the linear system of equations
  !
  !       C1*(X(I-1,J,K)-2.*X(I,J,K)+X(I+1,J,K))
  !     + C2*(X(I,J-1,K)-2.*X(I,J,K)+X(I,J+1,K))
  !     + A(K)*X(I,J,K-1)+B(K)*X(I,J,K)+C(K)*X(I,J,K+1) = F(I,J,K)
  !
  !     for  I=1,2,...,L, J=1,2,...,M, and K=1,2,...,N .
  !
  !     The indices K-1 and K+1 are evaluated modulo N, i.e.
  !     X(I,J,0) = X(I,J,N) and X(I,J,N+1) = X(I,J,1). The unknowns
  !     X(0,J,K), X(L+1,J,K), X(I,0,K), and X(I,M+1,K) are assumed to take
  !     on certain prescribed values described below.
  !
  !    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  !
  !
  !    * * * * * * * *    Parameter Description     * * * * * * * * * *
  !
  !
  !            * * * * * *   On Input    * * * * * *
  !
  !     LPEROD   Indicates the values that X(0,J,K) and X(L+1,J,K) are
  !              assumed to have.
  !
  !              = 0  If X(0,J,K) = X(L,J,K) and X(L+1,J,K) = X(1,J,K).
  !              = 1  If X(0,J,K) = X(L+1,J,K) = 0.
  !              = 2  If X(0,J,K) = 0  and X(L+1,J,K) = X(L-1,J,K).
  !              = 3  If X(0,J,K) = X(2,J,K) and X(L+1,J,K) = X(L-1,J,K).
  !              = 4  If X(0,J,K) = X(2,J,K) and X(L+1,J,K) = 0.
  !
  !     L        The number of unknowns in the I-direction. L must be at
  !              least 3.
  !
  !     C1       The real constant that appears in the above equation.
  !
  !     MPEROD   Indicates the values that X(I,0,K) and X(I,M+1,K) are
  !              assumed to have.
  !
  !              = 0  If X(I,0,K) = X(I,M,K) and X(I,M+1,K) = X(I,1,K).
  !              = 1  If X(I,0,K) = X(I,M+1,K) = 0.
  !              = 2  If X(I,0,K) = 0 and X(I,M+1,K) = X(I,M-1,K).
  !              = 3  If X(I,0,K) = X(I,2,K) and X(I,M+1,K) = X(I,M-1,K).
  !              = 4  If X(I,0,K) = X(I,2,K) and X(I,M+1,K) = 0.
  !
  !     M        The number of unknowns in the J-direction. M must be at
  !              least 3.
  !
  !     C2       The real constant which appears in the above equation.
  !
  !     NPEROD   = 0  If A(1) and C(N) are not zero.
  !              = 1  If A(1) = C(N) = 0.
  !
  !     N        The number of unknowns in the K-direction. N must be at
  !              least 3.
  !
  !
  !     A,B,C    One-dimensional arrays of length N that specify the
  !              coefficients in the linear equations given above.
  !
  !              If NPEROD = 0 the array elements must not depend upon the
  !              index K, but must be constant.  Specifically, the
  !              subroutine checks the following condition
  !
  !                          A(K) = C(1)
  !                          C(K) = C(1)
  !                          B(K) = B(1)
  !
  !                  for K=1,2,...,N.
  !
  !     LDIMF    The row (or first) dimension of the three-dimensional
  !              array F as it appears in the program calling POIS3D.
  !              This parameter is used to specify the variable dimension
  !              of F.  LDIMF must be at least L.
  !
  !     MDIMF    The column (or second) dimension of the three-dimensional
  !              array F as it appears in the program calling POIS3D.
  !              This parameter is used to specify the variable dimension
  !              of F.  MDIMF must be at least M.
  !
  !     F        A three-dimensional array that specifies the values of
  !              the right side of the linear system of equations given
  !              above.  F must be dimensioned at least L x M x N.
  !
  !     W        A one-dimensional array that must be provided by the
  !              user for work space.  The length of W must be at least
  !              30 + L + M + 2*N + MAX(L,M,N) +
  !              7*(INT((L+1)/2) + INT((M+1)/2)).
  !
  !
  !            * * * * * *   On Output   * * * * * *
  !
  !     F        Contains the solution X.
  !
  !     IERROR   An error flag that indicates invalid input parameters.
  !              Except for number zero, a solution is not attempted.
  !              = 0  No error
  !              = 1  If LPEROD .LT. 0 or .GT. 4
  !              = 2  If L .LT. 3
  !              = 3  If MPEROD .LT. 0 or .GT. 4
  !              = 4  If M .LT. 3
  !              = 5  If NPEROD .LT. 0 or .GT. 1
  !              = 6  If N .LT. 3
  !              = 7  If LDIMF .LT. L
  !              = 8  If MDIMF .LT. M
  !              = 9  If A(K) .NE. C(1) or C(K) .NE. C(1) or B(I) .NE.B(1)
  !                      for some K=1,2,...,N.
  !              = 10 If NPEROD = 1 and A(1) .NE. 0 or C(N) .NE. 0
  !
  !              Since this is the only means of indicating a possibly
  !              incorrect call to POIS3D, the user should test IERROR
  !              after the call.
  !
  !- Long Description:
  !
  !    * * * * * * *   Program Specifications    * * * * * * * * * * * *
  !
  !     Dimension of   A(N),B(N),C(N),F(LDIMF,MDIMF,N),
  !     Arguments      W(see argument list)
  !
  !     Latest         December 1, 1978
  !     Revision
  !
  !     Subprograms    POIS3D,POS3D1,TRIDQ,RFFTI,RFFTF,RFFTF1,RFFTB,
  !     Required       RFFTB1,COSTI,COST,SINTI,SINT,COSQI,COSQF,COSQF1
  !                    COSQB,COSQB1,SINQI,SINQF,SINQB,CFFTI,CFFTI1,
  !                    CFFTB,CFFTB1,PASSB2,PASSB3,PASSB4,PASSB,CFFTF,
  !                    CFFTF1,PASSF1,PASSF2,PASSF3,PASSF4,PASSF,PIMACH,
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
  !     Algorithm      This subroutine solves three-dimensional block
  !                    tridiagonal linear systems arising from finite
  !                    difference approximations to three-dimensional
  !                    Poisson equations using the Fourier transform
  !                    package FFTPAK written by Paul Swarztrauber.
  !
  !     Space          6561(decimal) = 14641(octal) locations on the
  !     Required       NCAR Control Data 7600
  !
  !     Timing and        The execution time T on the NCAR Control Data
  !     Accuracy       7600 for subroutine POIS3D is roughly proportional
  !                    to L*M*N*(log2(L)+log2(M)+5), but also depends on
  !                    input parameters LPEROD and MPEROD.  Some typical
  !                    values are listed in the table below when NPEROD=0.
  !                       To measure the accuracy of the algorithm a
  !                    uniform random number generator was used to create
  !                    a solution array X for the system given in the
  !                    'PURPOSE' with
  !
  !                       A(K) = C(K) = -0.5*B(K) = 1,       K=1,2,...,N
  !
  !                    and, when NPEROD = 1
  !
  !                       A(1) = C(N) = 0
  !                       A(N) = C(1) = 2.
  !
  !                    The solution X was substituted into the given sys-
  !                    tem and, using double precision, a right side Y was
  !                    computed.  Using this array Y subroutine POIS3D was
  !                    called to produce an approximate solution Z.  Then
  !                    the relative error, defined as
  !
  !                    E = MAX(ABS(Z(I,J,K)-X(I,J,K)))/MAX(ABS(X(I,J,K)))
  !
  !                    where the two maxima are taken over I=1,2,...,L,
  !                    J=1,2,...,M and K=1,2,...,N, was computed.  The
  !                    value of E is given in the table below for some
  !                    typical values of L,M and N.
  !
  !
  !                       L(=M=N)   LPEROD    MPEROD    T(MSECS)    E
  !                       ------    ------    ------    --------  ------
  !
  !                         16        0         0         272     1.E-13
  !                         15        1         1         287     4.E-13
  !                         17        3         3         338     2.E-13
  !                         32        0         0        1755     2.E-13
  !                         31        1         1        1894     2.E-12
  !                         33        3         3        2042     7.E-13
  !
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
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  POS3D1

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)

  INTEGER i, Ierror, iwbb, iwd, iwt, iwx, iwy, iwyrt, j, k, L, &
    Ldimf, lp, Lperod, M, Mdimf, mp, Mperod, N, nh
  REAL A(*), B(*), C(*), C1, C2, F(Ldimf,Mdimf,*), save(6), W(*)
  INTEGER nhm1, nhmk, nhpk, nodd, np, Nperod
  !* FIRST EXECUTABLE STATEMENT  POIS3D
  lp = Lperod + 1
  mp = Mperod + 1
  np = Nperod + 1
  !
  !     CHECK FOR INVALID INPUT.
  !
  Ierror = 0
  IF ( lp<1.OR.lp>5 ) Ierror = 1
  IF ( L<3 ) Ierror = 2
  IF ( mp<1.OR.mp>5 ) Ierror = 3
  IF ( M<3 ) Ierror = 4
  IF ( np<1.OR.np>2 ) Ierror = 5
  IF ( N<3 ) Ierror = 6
  IF ( Ldimf<L ) Ierror = 7
  IF ( Mdimf<M ) Ierror = 8
  IF ( np/=1 ) GOTO 200
  DO k = 1, N
    IF ( A(k)/=C(1) ) GOTO 100
    IF ( C(k)/=C(1) ) GOTO 100
    IF ( B(k)/=B(1) ) GOTO 100
  END DO
  GOTO 300
  100  Ierror = 9
  200 CONTINUE
  IF ( Nperod==1.AND.(A(1)/=0..OR.C(N)/=0.) ) Ierror = 10
  300 CONTINUE
  IF ( Ierror==0 ) THEN
    iwyrt = L + 1
    iwt = iwyrt + M
    iwd = iwt + MAX(L,M,N) + 1
    iwbb = iwd + N
    iwx = iwbb + N
    iwy = iwx + 7*((L+1)/2) + 15
    IF ( np/=2 ) THEN
      !
      !     REORDER UNKNOWNS WHEN NPEROD = 0.
      !
      nh = (N+1)/2
      nhm1 = nh - 1
      nodd = 1
      IF ( 2*nh==N ) nodd = 2
      DO i = 1, L
        DO j = 1, M
          DO k = 1, nhm1
            nhpk = nh + k
            nhmk = nh - k
            W(k) = F(i,j,nhmk) - F(i,j,nhpk)
            W(nhpk) = F(i,j,nhmk) + F(i,j,nhpk)
          END DO
          W(nh) = 2.*F(i,j,nh)
          IF ( nodd/=1 ) W(N) = 2.*F(i,j,N)
          DO k = 1, N
            F(i,j,k) = W(k)
          END DO
        END DO
      END DO
      save(1) = C(nhm1)
      save(2) = A(nh)
      save(3) = C(nh)
      save(4) = B(nhm1)
      save(5) = B(N)
      save(6) = A(N)
      C(nhm1) = 0.
      A(nh) = 0.
      C(nh) = 2.*C(nh)
      IF ( nodd==2 ) THEN
        A(N) = C(nh)
      ELSE
        B(nhm1) = B(nhm1) - A(nh-1)
        B(N) = B(N) + A(N)
      END IF
    END IF
    CALL POS3D1(lp,L,mp,M,N,A,B,C,Ldimf,Mdimf,F,W,W(iwyrt),W(iwt),W(iwd),&
      W(iwx),W(iwy),C1,C2,W(iwbb))
    IF ( np/=2 ) THEN
      DO i = 1, L
        DO j = 1, M
          DO k = 1, nhm1
            nhmk = nh - k
            nhpk = nh + k
            W(nhmk) = .5*(F(i,j,nhpk)+F(i,j,k))
            W(nhpk) = .5*(F(i,j,nhpk)-F(i,j,k))
          END DO
          W(nh) = .5*F(i,j,nh)
          IF ( nodd/=1 ) W(N) = .5*F(i,j,N)
          DO k = 1, N
            F(i,j,k) = W(k)
          END DO
        END DO
      END DO
      C(nhm1) = save(1)
      A(nh) = save(2)
      C(nh) = save(3)
      B(nhm1) = save(4)
      B(N) = save(5)
      A(N) = save(6)
    END IF
  END IF
END SUBROUTINE POIS3D
