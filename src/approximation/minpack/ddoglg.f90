!** DDOGLG
PURE SUBROUTINE DDOGLG(N,R,Lr,Diag,Qtb,Delta,X,Wa1,Wa2)
  !> Subsidiary to DNSQ and DNSQE
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (DOGLEG-S, DDOGLG-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
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
  !***
  ! **See also:**  DNSQ, DNSQE
  !***
  ! **Routines called:**  D1MACH, DENORM

  !* REVISION HISTORY  (YYMMDD)
  !   800301  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTIONsection.  (WRB)
  !   900328  Added TYPE section.  (WRB)
  USE service, ONLY : D1MACH
  !
  INTEGER, INTENT(IN) :: Lr, N
  REAL(DP), INTENT(IN) :: Delta
  REAL(DP), INTENT(IN) :: Diag(N), Qtb(N), R(Lr)
  REAL(DP), INTENT(OUT) :: Wa1(N), Wa2(N), X(N)
  !
  INTEGER :: i, j, jj, jp1, k, l
  REAL(DP) :: alpha, bnorm, epsmch, gnorm, qnorm, sgnorm, summ, temp
  !
  !     EPSMCH IS THE MACHINE PRECISION.
  !
  !* FIRST EXECUTABLE STATEMENT  DDOGLG
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
    summ = 0._DP
    IF( N>=jp1 ) THEN
      DO i = jp1, N
        summ = summ + R(l)*X(i)
        l = l + 1
      END DO
    END IF
    temp = R(jj)
    IF( temp==0._DP ) THEN
      l = j
      DO i = 1, j
        temp = MAX(temp,ABS(R(l)))
        l = l + N - i
      END DO
      temp = epsmch*temp
      IF( temp==0._DP ) temp = epsmch
    END IF
    X(j) = (Qtb(j)-summ)/temp
  END DO
  !
  !     TEST WHETHER THE GAUSS-NEWTON DIRECTION IS ACCEPTABLE.
  !
  DO j = 1, N
    Wa1(j) = 0._DP
    Wa2(j) = Diag(j)*X(j)
  END DO
  qnorm = NORM2(Wa2)
  IF( qnorm>Delta ) THEN
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
      END DO
      Wa1(j) = Wa1(j)/Diag(j)
    END DO
    !
    !     CALCULATE THE NORM OF THE SCALED GRADIENT AND TEST FOR
    !     THE SPECIAL CASE IN WHICH THE SCALED GRADIENT IS ZERO.
    !
    gnorm = NORM2(Wa1)
    sgnorm = 0._DP
    alpha = Delta/qnorm
    IF( gnorm/=0._DP ) THEN
      !
      !     CALCULATE THE POINT ALONG THE SCALED GRADIENT
      !     AT WHICH THE QUADRATIC IS MINIMIZED.
      !
      DO j = 1, N
        Wa1(j) = (Wa1(j)/gnorm)/Diag(j)
      END DO
      l = 1
      DO j = 1, N
        summ = 0._DP
        DO i = j, N
          summ = summ + R(l)*Wa1(i)
          l = l + 1
        END DO
        Wa2(j) = summ
      END DO
      temp = NORM2(Wa2)
      sgnorm = (gnorm/temp)/temp
      !
      !     TEST WHETHER THE SCALED GRADIENT DIRECTION IS ACCEPTABLE.
      !
      alpha = 0._DP
      IF( sgnorm<Delta ) THEN
        !
        !     THE SCALED GRADIENT DIRECTION IS NOT ACCEPTABLE.
        !     FINALLY, CALCULATE THE POINT ALONG THE DOGLEG
        !     AT WHICH THE QUADRATIC IS MINIMIZED.
        !
        bnorm = NORM2(Qtb)
        temp = (bnorm/gnorm)*(bnorm/qnorm)*(sgnorm/Delta)
        temp = temp - (Delta/qnorm)*(sgnorm/Delta)&
          **2 + SQRT((temp-(Delta/qnorm))**2+(1._DP-(Delta/qnorm)**2)&
          *(1._DP-(sgnorm/Delta)**2))
        alpha = ((Delta/qnorm)*(1._DP-(sgnorm/Delta)**2))/temp
      END IF
    END IF
    !
    !     FORM APPROPRIATE CONVEX COMBINATION OF THE GAUSS-NEWTON
    !     DIRECTION AND THE SCALED GRADIENT DIRECTION.
    !
    temp = (1._DP-alpha)*MIN(sgnorm,Delta)
    DO j = 1, N
      X(j) = temp*Wa1(j) + alpha*X(j)
    END DO
  END IF
  !
  !     LAST CARD OF SUBROUTINE DDOGLG.
  !
END SUBROUTINE DDOGLG