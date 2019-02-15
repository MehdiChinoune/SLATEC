!DECK DDOGLG
SUBROUTINE DDOGLG(N,R,Lr,Diag,Qtb,Delta,X,Wa1,Wa2)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DDOGLG
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DNSQ and DNSQE
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (DOGLEG-S, DDOGLG-D)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     Given an M by N matrix A, an N by N nonsingular diagonal
  !     matrix D, an M-vector B, and a positive number DELTA, the
  !     problem is to determine the convex combination X of the
  !     Gauss-Newton and scaled gradient directions that minimizes
  !     (A*X - B) in the least squares sense, subject to the
  !     restriction that the Euclidean norm of D*X be at most DELTA.
  !
  !     This subroutine completes the solution of the problem
  !     if it is provided with the necessary information from the
  !     QR factorization of A. That is, if A = Q*R, where Q has
  !     orthogonal columns and R is an upper triangular matrix,
  !     then DDOGLG expects the full upper triangle of R and
  !     the first N components of (Q transpose)*B.
  !
  !     The subroutine statement is
  !
  !       SUBROUTINE DDOGLG(N,R,LR,DIAG,QTB,DELTA,X,WA1,WA2)
  !
  !     where
  !
  !       N is a positive integer input variable set to the order of R.
  !
  !       R is an input array of length LR which must contain the upper
  !         triangular matrix R stored by rows.
  !
  !       LR is a positive integer input variable not less than
  !         (N*(N+1))/2.
  !
  !       DIAG is an input array of length N which must contain the
  !         diagonal elements of the matrix D.
  !
  !       QTB is an input array of length N which must contain the first
  !         N elements of the vector (Q transpose)*B.
  !
  !       DELTA is a positive input variable which specifies an upper
  !         bound on the Euclidean norm of D*X.
  !
  !       X is an output array of length N which contains the desired
  !         convex combination of the Gauss-Newton direction and the
  !         scaled gradient direction.
  !
  !       WA1 and WA2 are work arrays of length N.
  !
  !***SEE ALSO  DNSQ, DNSQE
  !***ROUTINES CALLED  D1MACH, DENORM
  !***REVISION HISTORY  (YYMMDD)
  !   800301  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   900328  Added TYPE section.  (WRB)
  !***END PROLOGUE  DDOGLG
  REAL(8) :: D1MACH, DENORM
  INTEGER i, j, jj, jp1, k, l, Lr, N
  REAL(8) :: alpha, bnorm, Delta, Diag(*), epsmch, gnorm, one, &
    qnorm, Qtb(*), R(*), sgnorm, sum, temp, Wa1(*), &
    Wa2(*), X(*), zero
  SAVE one, zero
  DATA one, zero/1.0D0, 0.0D0/
  !
  !     EPSMCH IS THE MACHINE PRECISION.
  !
  !***FIRST EXECUTABLE STATEMENT  DDOGLG
  epsmch = D1MACH(4)
  !
  !     FIRST, CALCULATE THE GAUSS-NEWTON DIRECTION.
  !
  jj = (N*(N+1))/2 + 1
  DO k = 1, N
    j = N - k + 1
    jp1 = j + 1
    jj = jj - k
    l = jj + 1
    sum = zero
    IF ( N>=jp1 ) THEN
      DO i = jp1, N
        sum = sum + R(l)*X(i)
        l = l + 1
      ENDDO
    ENDIF
    temp = R(jj)
    IF ( temp==zero ) THEN
      l = j
      DO i = 1, j
        temp = MAX(temp,ABS(R(l)))
        l = l + N - i
      ENDDO
      temp = epsmch*temp
      IF ( temp==zero ) temp = epsmch
    ENDIF
    X(j) = (Qtb(j)-sum)/temp
  ENDDO
  !
  !     TEST WHETHER THE GAUSS-NEWTON DIRECTION IS ACCEPTABLE.
  !
  DO j = 1, N
    Wa1(j) = zero
    Wa2(j) = Diag(j)*X(j)
  ENDDO
  qnorm = DENORM(N,Wa2)
  IF ( qnorm>Delta ) THEN
    !
    !     THE GAUSS-NEWTON DIRECTION IS NOT ACCEPTABLE.
    !     NEXT, CALCULATE THE SCALED GRADIENT DIRECTION.
    !
    l = 1
    DO j = 1, N
      temp = Qtb(j)
      DO i = j, N
        Wa1(i) = Wa1(i) + R(l)*temp
        l = l + 1
      ENDDO
      Wa1(j) = Wa1(j)/Diag(j)
    ENDDO
    !
    !     CALCULATE THE NORM OF THE SCALED GRADIENT AND TEST FOR
    !     THE SPECIAL CASE IN WHICH THE SCALED GRADIENT IS ZERO.
    !
    gnorm = DENORM(N,Wa1)
    sgnorm = zero
    alpha = Delta/qnorm
    IF ( gnorm/=zero ) THEN
      !
      !     CALCULATE THE POINT ALONG THE SCALED GRADIENT
      !     AT WHICH THE QUADRATIC IS MINIMIZED.
      !
      DO j = 1, N
        Wa1(j) = (Wa1(j)/gnorm)/Diag(j)
      ENDDO
      l = 1
      DO j = 1, N
        sum = zero
        DO i = j, N
          sum = sum + R(l)*Wa1(i)
          l = l + 1
        ENDDO
        Wa2(j) = sum
      ENDDO
      temp = DENORM(N,Wa2)
      sgnorm = (gnorm/temp)/temp
      !
      !     TEST WHETHER THE SCALED GRADIENT DIRECTION IS ACCEPTABLE.
      !
      alpha = zero
      IF ( sgnorm<Delta ) THEN
        !
        !     THE SCALED GRADIENT DIRECTION IS NOT ACCEPTABLE.
        !     FINALLY, CALCULATE THE POINT ALONG THE DOGLEG
        !     AT WHICH THE QUADRATIC IS MINIMIZED.
        !
        bnorm = DENORM(N,Qtb)
        temp = (bnorm/gnorm)*(bnorm/qnorm)*(sgnorm/Delta)
        temp = temp - (Delta/qnorm)*(sgnorm/Delta)&
          **2 + SQRT((temp-(Delta/qnorm))**2+(one-(Delta/qnorm)**2)&
          *(one-(sgnorm/Delta)**2))
        alpha = ((Delta/qnorm)*(one-(sgnorm/Delta)**2))/temp
      ENDIF
    ENDIF
    !
    !     FORM APPROPRIATE CONVEX COMBINATION OF THE GAUSS-NEWTON
    !     DIRECTION AND THE SCALED GRADIENT DIRECTION.
    !
    temp = (one-alpha)*MIN(sgnorm,Delta)
    DO j = 1, N
      X(j) = temp*Wa1(j) + alpha*X(j)
    ENDDO
  ENDIF
  !
  !     LAST CARD OF SUBROUTINE DDOGLG.
  !
END SUBROUTINE DDOGLG
