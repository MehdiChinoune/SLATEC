!DECK CMGNBN
SUBROUTINE CMGNBN(Nperod,N,Mperod,M,A,B,C,Idimy,Y,Ierror,W)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CMGNBN
  !***PURPOSE  Solve a complex block tridiagonal linear system of
  !            equations by a cyclic reduction algorithm.
  !***LIBRARY   SLATEC (FISHPACK)
  !***CATEGORY  I2B4B
  !***TYPE      COMPLEX (GENBUN-S, CMGNBN-C)
  !***KEYWORDS  CYCLIC REDUCTION, ELLIPTIC PDE, FISHPACK,
  !             TRIDIAGONAL LINEAR SYSTEM
  !***AUTHOR  Adams, J., (NCAR)
  !           Swarztrauber, P. N., (NCAR)
  !           Sweet, R., (NCAR)
  !***DESCRIPTION
  !
  !     Subroutine CMGNBN solves the complex linear system of equations
  !
  !          A(I)*X(I-1,J) + B(I)*X(I,J) + C(I)*X(I+1,J)
  !
  !          + X(I,J-1) - 2.*X(I,J) + X(I,J+1) = Y(I,J)
  !
  !               For I = 1,2,...,M  and  J = 1,2,...,N.
  !
  !     The indices I+1 and I-1 are evaluated modulo M, i.e.,
  !     X(0,J) = X(M,J) and X(M+1,J) = X(1,J), and X(I,0) may be equal to
  !     0, X(I,2), or X(I,N) and X(I,N+1) may be equal to 0, X(I,N-1), or
  !     X(I,1) depending on an input parameter.
  !
  !
  !     * * * * * * * *    Parameter Description     * * * * * * * * * *
  !
  !             * * * * * *   On Input    * * * * * *
  !
  !     NPEROD
  !       Indicates the values that X(I,0) and X(I,N+1) are assumed to
  !       have.
  !
  !       = 0  If X(I,0) = X(I,N) and X(I,N+1) = X(I,1).
  !       = 1  If X(I,0) = X(I,N+1) = 0  .
  !       = 2  If X(I,0) = 0 and X(I,N+1) = X(I,N-1).
  !       = 3  If X(I,0) = X(I,2) and X(I,N+1) = X(I,N-1).
  !       = 4  If X(I,0) = X(I,2) and X(I,N+1) = 0.
  !
  !     N
  !       The number of unknowns in the J-direction.  N must be greater
  !       than 2.
  !
  !     MPEROD
  !       = 0 If A(1) and C(M) are not zero
  !       = 1 If A(1) = C(M) = 0
  !
  !     M
  !       The number of unknowns in the I-direction.  N must be greater
  !       than 2.
  !
  !     A,B,C
  !       One-dimensional complex arrays of length M that specify the
  !       coefficients in the linear equations given above.  If MPEROD = 0
  !       the array elements must not depend upon the index I, but must be
  !       constant.  Specifically, the subroutine checks the following
  !       condition
  !
  !             A(I) = C(1)
  !             C(I) = C(1)
  !             B(I) = B(1)
  !
  !       For I=1,2,...,M.
  !
  !     IDIMY
  !       The row (or first) dimension of the two-dimensional array Y as
  !       it appears in the program calling CMGNBN.  This parameter is
  !       used to specify the variable dimension of Y.  IDIMY must be at
  !       least M.
  !
  !     Y
  !       A two-dimensional complex array that specifies the values of the
  !       right side of the linear system of equations given above.  Y
  !       must be dimensioned at least M*N.
  !
  !     W
  !       A one-dimensional complex array that must be provided by the
  !       user for work space.  W may require up to 4*N +
  !       (10 + INT(log2(N)))*M LOCATIONS.  The actual number of locations
  !       used is computed by CMGNBN and is returned in location W(1).
  !
  !
  !             * * * * * *   On Output     * * * * * *
  !
  !     Y
  !       Contains the solution X.
  !
  !     IERROR
  !       An error flag which indicates invalid input parameters.  Except
  !       for number zero, a solution is not attempted.
  !
  !       = 0  No error.
  !       = 1  M .LE. 2
  !       = 2  N .LE. 2
  !       = 3  IDIMY .LT. M
  !       = 4  NPEROD .LT. 0 or NPEROD .GT. 4
  !       = 5  MPEROD .LT. 0 or MPEROD .GT. 1
  !       = 6  A(I) .NE. C(1) or C(I) .NE. C(1) or B(I) .NE. B(1) for
  !            some I=1,2,...,M.
  !       = 7  A(1) .NE. 0 or C(M) .NE. 0 and MPEROD = 1
  !
  !     W
  !       W(1) contains the required length of W.
  !
  ! *Long Description:
  !
  !     * * * * * * *   Program Specifications    * * * * * * * * * * * *
  !
  !     Dimension of   A(M),B(M),C(M),Y(IDIMY,N),W(see parameter list)
  !     Arguments
  !
  !     Latest         June 1979
  !     Revision
  !
  !     Subprograms    CMGNBN,CMPOSD,CMPOSN,CMPOSP,CMPCSG,CMPMRG,
  !     Required       CMPTRX,CMPTR3,PIMACH
  !
  !     Special        None
  !     Conditions
  !
  !     Common         None
  !     Blocks
  !
  !     I/O            None
  !
  !     Precision      Single
  !
  !     Specialist     Roland Sweet
  !
  !     Language       FORTRAN
  !
  !     History        Written by Roland Sweet at NCAR in June, 1977
  !
  !     Algorithm      The linear system is solved by a cyclic reduction
  !                    algorithm described in the reference.
  !
  !     Space          4944(DECIMAL) = 11520(octal) locations on the NCAR
  !     Required       Control Data 7600
  !
  !     Timing and      The execution time T on the NCAR Control Data
  !     Accuracy       7600 for subroutine CMGNBN is roughly proportional
  !                    to M*N*log2(N), but also depends on the input
  !                    parameter NPEROD.  Some typical values are listed
  !                    in the table below.
  !                       To measure the accuracy of the algorithm a
  !                    uniform random number generator was used to create
  !                    a solution array X for the system given in the
  !                    'PURPOSE' with
  !
  !                       A(I) = C(I) = -0.5*B(I) = 1,       I=1,2,...,M
  !
  !                    and, when MPEROD = 1
  !
  !                       A(1) = C(M) = 0
  !                       A(M) = C(1) = 2.
  !
  !                    The solution X was substituted into the given sys-
  !                    tem and a right side Y was computed.  Using this
  !                    array Y subroutine CMGNBN was called to produce an
  !                    approximate solution Z.  Then the relative error,
  !                    defined as
  !
  !                       E = MAX(ABS(Z(I,J)-X(I,J)))/MAX(ABS(X(I,J)))
  !
  !                    where the two maxima are taken over all I=1,2,...,M
  !                    and J=1,2,...,N, was computed.  The value of E is
  !                    given in the table below for some typical values of
  !                    M and N.
  !
  !
  !                       M (=N)    MPEROD    NPEROD    T(MSECS)    E
  !                       ------    ------    ------    --------  ------
  !
  !                         31        0         0          77     1.E-12
  !                         31        1         1          45     4.E-13
  !                         31        1         3          91     2.E-12
  !                         32        0         0          59     7.E-14
  !                         32        1         1          65     5.E-13
  !                         32        1         3          97     2.E-13
  !                         33        0         0          80     6.E-13
  !                         33        1         1          67     5.E-13
  !                         33        1         3          76     3.E-12
  !                         63        0         0         350     5.E-12
  !                         63        1         1         215     6.E-13
  !                         63        1         3         412     1.E-11
  !                         64        0         0         264     1.E-13
  !                         64        1         1         287     3.E-12
  !                         64        1         3         421     3.E-13
  !                         65        0         0         338     2.E-12
  !                         65        1         1         292     5.E-13
  !                         65        1         3         329     1.E-11
  !
  !     Portability    American National Standards Institute Fortran.
  !                    The machine dependent constant PI is defined in
  !                    function PIMACH.
  !
  !     Required       COS
  !     Resident
  !     Routines
  !
  !     Reference      Sweet, R., 'A Cyclic Reduction Algorithm for
  !                    Solving Block Tridiagonal Systems Of Arbitrary
  !                    Dimensions,' SIAM J. on Numer. Anal.,
  !                    14(SEPT., 1977), PP. 706-720.
  !
  !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  !
  !***REFERENCES  R. Sweet, A cyclic reduction algorithm for solving
  !                 block tridiagonal systems of arbitrary dimensions,
  !                 SIAM Journal on Numerical Analysis 14, (September
  !                 1977), pp. 706-720.
  !***ROUTINES CALLED  CMPOSD, CMPOSN, CMPOSP
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  CMGNBN
  INTEGER i, Idimy, Ierror, ipstor, irev, iwb2, iwb3, iwba, iwbb, &
    iwbc, iwd, iwp, iwtcos, iww1, iww2, iww3, j, k, M, mh
  INTEGER mhm1, mhmi, mhpi, modd, mp, Mperod, mskip, N, nby2, np, Nperod
  COMPLEX A, B, C, Y, W, a1
  DIMENSION Y(Idimy,*)
  DIMENSION W(*), B(*), A(*), C(*)
  !***FIRST EXECUTABLE STATEMENT  CMGNBN
  Ierror = 0
  IF ( M<=2 ) Ierror = 1
  IF ( N<=2 ) Ierror = 2
  IF ( Idimy<M ) Ierror = 3
  IF ( Nperod<0.OR.Nperod>4 ) Ierror = 4
  IF ( Mperod<0.OR.Mperod>1 ) Ierror = 5
  IF ( Mperod==1 ) THEN
    IF ( ABS(A(1))/=0..AND.ABS(C(M))/=0. ) Ierror = 7
  ELSE
    DO i = 2, M
      IF ( ABS(A(i)-C(1))/=0. ) GOTO 100
      IF ( ABS(C(i)-C(1))/=0. ) GOTO 100
      IF ( ABS(B(i)-B(1))/=0. ) GOTO 100
    ENDDO
  ENDIF
  GOTO 200
  100  Ierror = 6
  200 CONTINUE
  IF ( Ierror/=0 ) RETURN
  iwba = M + 1
  iwbb = iwba + M
  iwbc = iwbb + M
  iwb2 = iwbc + M
  iwb3 = iwb2 + M
  iww1 = iwb3 + M
  iww2 = iww1 + M
  iww3 = iww2 + M
  iwd = iww3 + M
  iwtcos = iwd + M
  iwp = iwtcos + 4*N
  DO i = 1, M
    k = iwba + i - 1
    W(k) = -A(i)
    k = iwbc + i - 1
    W(k) = -C(i)
    k = iwbb + i - 1
    W(k) = 2. - B(i)
    DO j = 1, N
      Y(i,j) = -Y(i,j)
    ENDDO
  ENDDO
  mp = Mperod + 1
  np = Nperod + 1
  IF ( mp==1 ) GOTO 600
  300 CONTINUE
  SELECT CASE (np)
    CASE (2)
      CALL CMPOSD(M,N,1,W(iwba),W(iwbb),W(iwbc),Y,Idimy,W,W(iww1),W(iwd),&
        W(iwtcos),W(iwp))
    CASE (3)
      CALL CMPOSN(M,N,1,2,W(iwba),W(iwbb),W(iwbc),Y,Idimy,W,W(iwb2),W(iwb3),&
        W(iww1),W(iww2),W(iww3),W(iwd),W(iwtcos),W(iwp))
    CASE (4)
      CALL CMPOSN(M,N,1,1,W(iwba),W(iwbb),W(iwbc),Y,Idimy,W,W(iwb2),W(iwb3),&
        W(iww1),W(iww2),W(iww3),W(iwd),W(iwtcos),W(iwp))
    CASE (5)
      !
      !     REVERSE COLUMNS WHEN NPEROD = 4
      !
      irev = 1
      nby2 = N/2
      GOTO 700
    CASE DEFAULT
      CALL CMPOSP(M,N,W(iwba),W(iwbb),W(iwbc),Y,Idimy,W,W(iwb2),W(iwb3),&
        W(iww1),W(iww2),W(iww3),W(iwd),W(iwtcos),W(iwp))
  END SELECT
  400  ipstor = INT(W(iww1))
  irev = 2
  IF ( Nperod==4 ) GOTO 700
  500 CONTINUE
  IF ( mp==1 ) GOTO 800
  IF ( mp==2 ) GOTO 900
  !
  !     REORDER UNKNOWNS WHEN MP =0
  !
  600  mh = (M+1)/2
  mhm1 = mh - 1
  modd = 1
  IF ( mh*2==M ) modd = 2
  DO j = 1, N
    DO i = 1, mhm1
      mhpi = mh + i
      mhmi = mh - i
      W(i) = Y(mhmi,j) - Y(mhpi,j)
      W(mhpi) = Y(mhmi,j) + Y(mhpi,j)
    ENDDO
    W(mh) = 2.*Y(mh,j)
    IF ( modd/=1 ) W(M) = 2.*Y(M,j)
    DO i = 1, M
      Y(i,j) = W(i)
    ENDDO
  ENDDO
  k = iwbc + mhm1 - 1
  i = iwba + mhm1
  W(k) = (0.,0.)
  W(i) = (0.,0.)
  W(k+1) = 2.*W(k+1)
  IF ( modd==2 ) THEN
    W(iwbb-1) = W(k+1)
  ELSE
    k = iwbb + mhm1 - 1
    W(k) = W(k) - W(i-1)
    W(iwbc-1) = W(iwbc-1) + W(iwbb-1)
  ENDIF
  GOTO 300
  700 CONTINUE
  DO j = 1, nby2
    mskip = N + 1 - j
    DO i = 1, M
      a1 = Y(i,j)
      Y(i,j) = Y(i,mskip)
      Y(i,mskip) = a1
    ENDDO
  ENDDO
  IF ( irev==1 ) THEN
    CALL CMPOSN(M,N,1,2,W(iwba),W(iwbb),W(iwbc),Y,Idimy,W,W(iwb2),W(iwb3),&
      W(iww1),W(iww2),W(iww3),W(iwd),W(iwtcos),W(iwp))
    GOTO 400
  ELSEIF ( irev==2 ) THEN
    GOTO 500
  ENDIF
  800 CONTINUE
  DO j = 1, N
    DO i = 1, mhm1
      mhmi = mh - i
      mhpi = mh + i
      W(mhmi) = .5*(Y(mhpi,j)+Y(i,j))
      W(mhpi) = .5*(Y(mhpi,j)-Y(i,j))
    ENDDO
    W(mh) = .5*Y(mh,j)
    IF ( modd/=1 ) W(M) = .5*Y(M,j)
    DO i = 1, M
      Y(i,j) = W(i)
    ENDDO
  ENDDO
  !
  !     RETURN STORAGE REQUIREMENTS FOR W ARRAY.
  !
  900  W(1) = CMPLX(REAL(ipstor+iwp-1),0.)
END SUBROUTINE CMGNBN
