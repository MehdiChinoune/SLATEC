!*==BINTK.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK BINTK
SUBROUTINE BINTK(X,Y,T,N,K,Bcoef,Q,Work)
  IMPLICIT NONE
  !*--BINTK5
  !***BEGIN PROLOGUE  BINTK
  !***PURPOSE  Compute the B-representation of a spline which interpolates
  !            given data.
  !***LIBRARY   SLATEC
  !***CATEGORY  E1A
  !***TYPE      SINGLE PRECISION (BINTK-S, DBINTK-D)
  !***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION
  !***AUTHOR  Amos, D. E., (SNLA)
  !***DESCRIPTION
  !
  !     Written by Carl de Boor and modified by D. E. Amos
  !
  !     Abstract
  !
  !         BINTK is the SPLINT routine of the reference.
  !
  !         BINTK produces the B-spline coefficients, BCOEF, of the
  !         B-spline of order K with knots T(I), I=1,...,N+K, which
  !         takes on the value Y(I) at X(I), I=1,...,N.  The spline or
  !         any of its derivatives can be evaluated by calls to BVALU.
  !         The I-th equation of the linear system A*BCOEF = B for the
  !         coefficients of the interpolant enforces interpolation at
  !         X(I)), I=1,...,N.  Hence, B(I) = Y(I), all I, and A is
  !         a band matrix with 2K-1 bands if A is invertible. The matrix
  !         A is generated row by row and stored, diagonal by diagonal,
  !         in the rows of Q, with the main diagonal going into row K.
  !         The banded system is then solved by a call to BNFAC (which
  !         constructs the triangular factorization for A and stores it
  !         again in Q), followed by a call to BNSLV (which then
  !         obtains the solution BCOEF by substitution). BNFAC does no
  !         pivoting, since the total positivity of the matrix A makes
  !         this unnecessary.  The linear system to be solved is
  !         (theoretically) invertible if and only if
  !                 T(I) .LT. X(I)) .LT. T(I+K),        all I.
  !         Equality is permitted on the left for I=1 and on the right
  !         for I=N when K knots are used at X(1) or X(N).  Otherwise,
  !         violation of this condition is certain to lead to an error.
  !
  !     Description of Arguments
  !         Input
  !           X       - vector of length N containing data point abscissa
  !                     in strictly increasing order.
  !           Y       - corresponding vector of length N containing data
  !                     point ordinates.
  !           T       - knot vector of length N+K
  !                     since T(1),..,T(K) .LE. X(1) and T(N+1),..,T(N+K)
  !                     .GE. X(N), this leaves only N-K knots (not nec-
  !                     essarily X(I)) values) interior to (X(1),X(N))
  !           N       - number of data points, N .GE. K
  !           K       - order of the spline, K .GE. 1
  !
  !         Output
  !           BCOEF   - a vector of length N containing the B-spline
  !                     coefficients
  !           Q       - a work vector of length (2*K-1)*N, containing
  !                     the triangular factorization of the coefficient
  !                     matrix of the linear system being solved.  The
  !                     coefficients for the interpolant of an
  !                     additional data set (X(I)),YY(I)), I=1,...,N
  !                     with the same abscissa can be obtained by loading
  !                     YY into BCOEF and then executing
  !                         CALL BNSLV (Q,2K-1,N,K-1,K-1,BCOEF)
  !           WORK    - work vector of length 2*K
  !
  !     Error Conditions
  !         Improper  input is a fatal error
  !         Singular system of equations is a fatal error
  !
  !***REFERENCES  D. E. Amos, Computation with splines and B-splines,
  !                 Report SAND78-1968, Sandia Laboratories, March 1979.
  !               Carl de Boor, Package for calculating with B-splines,
  !                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
  !                 pp. 441-472.
  !               Carl de Boor, A Practical Guide to Splines, Applied
  !                 Mathematics Series 27, Springer-Verlag, New York,
  !                 1978.
  !***ROUTINES CALLED  BNFAC, BNSLV, BSPVN, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   800901  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  BINTK
  !
  INTEGER iflag , iwork , K , N , i , ilp1mx , j , jj , km1 , kpkm2 , left , &
    lenq , np1
  REAL Bcoef(*) , Y(*) , Q(*) , T(*) , X(*) , xi , Work(*)
  !     DIMENSION Q(2*K-1,N), T(N+K)
  !***FIRST EXECUTABLE STATEMENT  BINTK
  IF ( K<1 ) THEN
    CALL XERMSG('SLATEC','BINTK','K DOES NOT SATISFY K.GE.1',2,1)
    RETURN
  ELSE
    IF ( N<K ) THEN
      CALL XERMSG('SLATEC','BINTK','N DOES NOT SATISFY N.GE.K',2,1)
      RETURN
    ELSE
      jj = N - 1
      IF ( jj/=0 ) THEN
        DO i = 1 , jj
          IF ( X(i)>=X(i+1) ) GOTO 50
        ENDDO
      ENDIF
      np1 = N + 1
      km1 = K - 1
      kpkm2 = 2*km1
      left = K
      !                ZERO OUT ALL ENTRIES OF Q
      lenq = N*(K+km1)
      DO i = 1 , lenq
        Q(i) = 0.0E0
      ENDDO
      !
      !  ***   LOOP OVER I TO CONSTRUCT THE  N  INTERPOLATION EQUATIONS
      DO i = 1 , N
        xi = X(i)
        ilp1mx = MIN(i+K,np1)
        !        *** FIND  LEFT  IN THE CLOSED INTERVAL (I,I+K-1) SUCH THAT
        !                T(LEFT) .LE. X(I) .LT. T(LEFT+1)
        !        MATRIX IS SINGULAR IF THIS IS NOT POSSIBLE
        left = MAX(left,i)
        IF ( xi<T(left) ) GOTO 20
        DO WHILE ( xi>=T(left+1) )
          left = left + 1
          IF ( left>=ilp1mx ) THEN
            left = left - 1
            IF ( xi<=T(left+1) ) EXIT
            GOTO 20
          ENDIF
        ENDDO
        !        *** THE I-TH EQUATION ENFORCES INTERPOLATION AT XI, HENCE
        !        A(I,J) = B(J,K,T)(XI), ALL J. ONLY THE  K  ENTRIES WITH  J =
        !        LEFT-K+1,...,LEFT ACTUALLY MIGHT BE NONZERO. THESE  K  NUMBERS
        !        ARE RETURNED, IN  BCOEF (USED FOR TEMP. STORAGE HERE), BY THE
        !        FOLLOWING
        CALL BSPVN(T,K,K,1,xi,left,Bcoef,Work,iwork)
        !        WE THEREFORE WANT  BCOEF(J) = B(LEFT-K+J)(XI) TO GO INTO
        !        A(I,LEFT-K+J), I.E., INTO  Q(I-(LEFT+J)+2*K,(LEFT+J)-K) SINCE
        !        A(I+J,J)  IS TO GO INTO  Q(I+K,J), ALL I,J,  IF WE CONSIDER  Q
        !        AS A TWO-DIM. ARRAY , WITH  2*K-1  ROWS (SEE COMMENTS IN
        !        BNFAC). IN THE PRESENT PROGRAM, WE TREAT  Q  AS AN EQUIVALENT
        !        ONE-DIMENSIONAL ARRAY (BECAUSE OF FORTRAN RESTRICTIONS ON
        !        DIMENSION STATEMENTS) . WE THEREFORE WANT  BCOEF(J) TO GO INTO
        !        ENTRY
        !            I -(LEFT+J) + 2*K + ((LEFT+J) - K-1)*(2*K-1)
        !                   =  I-LEFT+1 + (LEFT -K)*(2*K-1) + (2*K-2)*J
        !        OF  Q .
        jj = i - left + 1 + (left-K)*(K+km1)
        DO j = 1 , K
          jj = jj + kpkm2
          Q(jj) = Bcoef(j)
        ENDDO
      ENDDO
      !
      !     ***OBTAIN FACTORIZATION OF  A  , STORED AGAIN IN  Q.
      CALL BNFAC(Q,K+km1,N,km1,km1,iflag)
      IF ( iflag==2 ) THEN
        CALL XERMSG('SLATEC','BINTK',&
          'THE SYSTEM OF SOLVER DETECTS A SINGULAR SYSTEM ALTHOUGH '&
          //&
          'THE THEORETICAL CONDITIONS FOR A SOLUTION WERE SATISFIED.'&
          ,8,1)
        RETURN
      ELSE
        !     *** SOLVE  A*BCOEF = Y  BY BACKSUBSTITUTION
        DO i = 1 , N
          Bcoef(i) = Y(i)
        ENDDO
        CALL BNSLV(Q,K+km1,N,km1,km1,Bcoef)
        RETURN
      ENDIF
      !
      !
      20       CALL XERMSG('SLATEC','BINTK',&
        'SOME ABSCISSA WAS NOT IN THE SUPPORT OF THE CORRESPONDING '&
        //'BASIS FUNCTION AND THE SYSTEM IS SINGULAR.',2,1)
      RETURN
    ENDIF
    50     CALL XERMSG('SLATEC','BINTK',&
      'X(I) DOES NOT SATISFY X(I).LT.X(I+1) FOR SOME I',2,1)
  ENDIF
END SUBROUTINE BINTK
