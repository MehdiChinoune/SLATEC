!*==GENBUN.f90  processed by SPAG 6.72Dc at 10:55 on  6 Feb 2019
!DECK GENBUN
SUBROUTINE GENBUN(Nperod,N,Mperod,M,A,B,C,Idimy,Y,Ierror,W)
  IMPLICIT NONE
  !*--GENBUN5
  !*** Start of declarations inserted by SPAG
  REAL A , a1 , B , C , W , Y
  INTEGER i , Idimy , Ierror , ipstor , irev , iwb2 , iwb3 , iwba , iwbb , &
    iwbc , iwd , iwp , iwtcos , iww1 , iww2 , iww3 , j , k , M , mh
  INTEGER mhm1 , mhmi , mhpi , modd , mp , mp1 , Mperod , mskip , N , nby2 , &
    np , Nperod
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  GENBUN
  !***PURPOSE  Solve by a cyclic reduction algorithm the linear system
  !            of equations that results from a finite difference
  !            approximation to certain 2-d elliptic PDE's on a centered
  !            grid .
  !***LIBRARY   SLATEC (FISHPACK)
  !***CATEGORY  I2B4B
  !***TYPE      SINGLE PRECISION (GENBUN-S, CMGNBN-C)
  !***KEYWORDS  ELLIPTIC, FISHPACK, PDE, TRIDIAGONAL
  !***AUTHOR  Adams, J., (NCAR)
  !           Swarztrauber, P. N., (NCAR)
  !           Sweet, R., (NCAR)
  !***DESCRIPTION
  !
  !     Subroutine GENBUN solves the linear system of equations
  !
  !          A(I)*X(I-1,J) + B(I)*X(I,J) + C(I)*X(I+1,J)
  !
  !          + X(I,J-1) - 2.*X(I,J) + X(I,J+1) = Y(I,J)
  !
  !               for I = 1,2,...,M  and  J = 1,2,...,N.
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
  !       = 0 if A(1) and C(M) are not zero.
  !       = 1 if A(1) = C(M) = 0.
  !
  !     M
  !       The number of unknowns in the I-direction.  M must be greater
  !       than 2.
  !
  !     A,B,C
  !       One-dimensional arrays of length M that specify the
  !       coefficients in the linear equations given above.  If MPEROD = 0
  !       the array elements must not depend upon the index I, but must be
  !       constant.  Specifically, the subroutine checks the following
  !       condition
  !
  !             A(I) = C(1)
  !             C(I) = C(1)
  !             B(I) = B(1)
  !
  !       for I=1,2,...,M.
  !
  !     IDIMY
  !       The row (or first) dimension of the two-dimensional array Y as
  !       it appears in the program calling GENBUN.  This parameter is
  !       used to specify the variable dimension of Y.  IDIMY must be at
  !       least M.
  !
  !     Y
  !       A two-dimensional array that specifies the values of the right
  !       side of the linear system of equations given above.  Y must be
  !       dimensioned at least M*N.
  !
  !     W
  !       A one-dimensional array that must be provided by the user for
  !       work space.  W may require up to 4*N + (10 + INT(log2(N)))*M
  !       locations.  The actual number of locations used is computed by
  !       GENBUN and is returned in location W(1).
  !
  !
  !             * * * * * *   On Output     * * * * * *
  !
  !     Y
  !       Contains the solution X.
  !
  !     IERROR
  !       An error flag that indicates invalid input parameters.  Except
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
  !     Latest         June 1, 1976
  !     Revision
  !
  !     Subprograms    GENBUN,POISD2,POISN2,POISP2,COSGEN,MERGE,TRIX,TRI3,
  !     Required       PIMACH
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
  !     History        Standardized April 1, 1973
  !                    Revised August 20,1973
  !                    Revised January 1, 1976
  !
  !     Algorithm      The linear system is solved by a cyclic reduction
  !                    algorithm described in the reference.
  !
  !     Space          4944(decimal) = 11520(octal) locations on the NCAR
  !     Required       Control Data 7600.
  !
  !     Timing and        The execution time T on the NCAR Control Data
  !     Accuracy       7600 for subroutine GENBUN is roughly proportional
  !                    to M*N*log2(N), but also depends on the input
  !                    parameter NPEROD.  Some typical values are listed
  !                    in the table below.  More comprehensive timing
  !                    charts may be found in the reference.
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
  !                    tem and, using double precision, a right side Y was
  !                    computed.  Using this array Y subroutine GENBUN was
  !                    called to produce an approximate solution Z.  Then
  !                    the relative error, defined as
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
  !                         31        0         0          36     6.E-14
  !                         31        1         1          21     4.E-13
  !                         31        1         3          41     3.E-13
  !                         32        0         0          29     9.E-14
  !                         32        1         1          32     3.E-13
  !                         32        1         3          48     1.E-13
  !                         33        0         0          36     9.E-14
  !                         33        1         1          30     4.E-13
  !                         33        1         3          34     1.E-13
  !                         63        0         0         150     1.E-13
  !                         63        1         1          91     1.E-12
  !                         63        1         3         173     2.E-13
  !                         64        0         0         122     1.E-13
  !                         64        1         1         128     1.E-12
  !                         64        1         3         199     6.E-13
  !                         65        0         0         143     2.E-13
  !                         65        1         1         120     1.E-12
  !                         65        1         3         138     4.E-13
  !
  !     Portability    American National Standards Institute Fortran.
  !                    The machine dependent constant PI is defined in
  !                    function PIMACH.
  !
  !     Required       COS
  !     Resident
  !     Routines
  !
  !     Reference      Sweet, R., 'A Cyclic Reduction Algorithm For
  !                    Solving Block Tridiagonal Systems Of Arbitrary
  !                    Dimensions,' SIAM J. on Numer. Anal.,
  !                    14(Sept., 1977), PP. 706-720.
  !
  !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  !
  !***REFERENCES  R. Sweet, A cyclic reduction algorithm for solving
  !                 block tridiagonal systems of arbitrary dimensions,
  !                 SIAM Journal on Numerical Analysis 14, (September
  !                 1977), pp. 706-720.
  !***ROUTINES CALLED  POISD2, POISN2, POISP2
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  GENBUN
  !
  !
  DIMENSION Y(Idimy,*)
  DIMENSION W(*) , B(*) , A(*) , C(*)
  !***FIRST EXECUTABLE STATEMENT  GENBUN
  Ierror = 0
  IF ( M<=2 ) Ierror = 1
  IF ( N<=2 ) Ierror = 2
  IF ( Idimy<M ) Ierror = 3
  IF ( Nperod<0.OR.Nperod>4 ) Ierror = 4
  IF ( Mperod<0.OR.Mperod>1 ) Ierror = 5
  IF ( Mperod==1 ) THEN
    IF ( A(1)/=0..OR.C(M)/=0. ) Ierror = 7
  ELSE
    DO i = 2 , M
      IF ( A(i)/=C(1) ) GOTO 100
      IF ( C(i)/=C(1) ) GOTO 100
      IF ( B(i)/=B(1) ) GOTO 100
    ENDDO
  ENDIF
  GOTO 200
  100  Ierror = 6
  200  IF ( Ierror/=0 ) RETURN
  mp1 = M + 1
  iwba = mp1
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
  DO i = 1 , M
    k = iwba + i - 1
    W(k) = -A(i)
    k = iwbc + i - 1
    W(k) = -C(i)
    k = iwbb + i - 1
    W(k) = 2. - B(i)
    DO j = 1 , N
      Y(i,j) = -Y(i,j)
    ENDDO
  ENDDO
  mp = Mperod + 1
  np = Nperod + 1
  IF ( mp==1 ) GOTO 600
  300  SELECT CASE (np)
CASE (2)
  CALL POISD2(M,N,1,W(iwba),W(iwbb),W(iwbc),Y,Idimy,W,W(iww1),W(iwd),&
    W(iwtcos),W(iwp))
CASE (3)
  CALL POISN2(M,N,1,2,W(iwba),W(iwbb),W(iwbc),Y,Idimy,W,W(iwb2),W(iwb3),&
    W(iww1),W(iww2),W(iww3),W(iwd),W(iwtcos),W(iwp))
CASE (4)
  CALL POISN2(M,N,1,1,W(iwba),W(iwbb),W(iwbc),Y,Idimy,W,W(iwb2),W(iwb3),&
    W(iww1),W(iww2),W(iww3),W(iwd),W(iwtcos),W(iwp))
CASE (5)
  !
  !     REVERSE COLUMNS WHEN NPEROD = 4.
  !
  irev = 1
  nby2 = N/2
  GOTO 700
CASE DEFAULT
  CALL POISP2(M,N,W(iwba),W(iwbb),W(iwbc),Y,Idimy,W,W(iwb2),W(iwb3),&
    W(iww1),W(iww2),W(iww3),W(iwd),W(iwtcos),W(iwp))
END SELECT
400  ipstor = W(iww1)
irev = 2
IF ( Nperod==4 ) GOTO 700
500  IF ( mp==1 ) GOTO 800
IF ( mp==2 ) GOTO 900
!
!     REORDER UNKNOWNS WHEN MP =0
!
600  mh = (M+1)/2
mhm1 = mh - 1
modd = 1
IF ( mh*2==M ) modd = 2
DO j = 1 , N
DO i = 1 , mhm1
  mhpi = mh + i
  mhmi = mh - i
  W(i) = Y(mhmi,j) - Y(mhpi,j)
  W(mhpi) = Y(mhmi,j) + Y(mhpi,j)
ENDDO
W(mh) = 2.*Y(mh,j)
IF ( modd/=1 ) W(M) = 2.*Y(M,j)
DO i = 1 , M
  Y(i,j) = W(i)
ENDDO
ENDDO
k = iwbc + mhm1 - 1
i = iwba + mhm1
W(k) = 0.
W(i) = 0.
W(k+1) = 2.*W(k+1)
IF ( modd==2 ) THEN
W(iwbb-1) = W(k+1)
ELSE
k = iwbb + mhm1 - 1
W(k) = W(k) - W(i-1)
W(iwbc-1) = W(iwbc-1) + W(iwbb-1)
ENDIF
GOTO 300
700  DO j = 1 , nby2
mskip = N + 1 - j
DO i = 1 , M
a1 = Y(i,j)
Y(i,j) = Y(i,mskip)
Y(i,mskip) = a1
ENDDO
ENDDO
IF ( irev==1 ) THEN
CALL POISN2(M,N,1,2,W(iwba),W(iwbb),W(iwbc),Y,Idimy,W,W(iwb2),W(iwb3),&
W(iww1),W(iww2),W(iww3),W(iwd),W(iwtcos),W(iwp))
GOTO 400
ELSEIF ( irev==2 ) THEN
GOTO 500
ENDIF
800  DO j = 1 , N
DO i = 1 , mhm1
mhmi = mh - i
mhpi = mh + i
W(mhmi) = .5*(Y(mhpi,j)+Y(i,j))
W(mhpi) = .5*(Y(mhpi,j)-Y(i,j))
ENDDO
W(mh) = .5*Y(mh,j)
IF ( modd/=1 ) W(M) = .5*Y(M,j)
DO i = 1 , M
Y(i,j) = W(i)
ENDDO
ENDDO
!
!     RETURN STORAGE REQUIREMENTS FOR W ARRAY.
!
900  W(1) = ipstor + iwp - 1
END SUBROUTINE GENBUN
