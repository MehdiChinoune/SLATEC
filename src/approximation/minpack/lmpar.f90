!** LMPAR
SUBROUTINE LMPAR(N,R,Ldr,Ipvt,Diag,Qtb,Delta,Par,X,Sigma,Wa1,Wa2)
  !>
  !  Subsidiary to SNLS1 and SNLS1E
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (LMPAR-S, DMPAR-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     Given an M by N matrix A, an N by N nonsingular DIAGONAL
  !     matrix D, an M-vector B, and a positive number DELTA,
  !     the problem is to determine a value for the parameter
  !     PAR such that if X solves the system
  !
  !           A*X = B,     SQRT(PAR)*D*X = 0 ,
  !
  !     in the least squares sense, and DXNORM is the Euclidean
  !     norm of D*X, then either PAR is zero and
  !
  !           (DXNORM-DELTA) .LE. 0.1*DELTA ,
  !
  !     or PAR is positive and
  !
  !           ABS(DXNORM-DELTA) .LE. 0.1*DELTA .
  !
  !     This subroutine completes the solution of the problem
  !     if it is provided with the necessary information from the
  !     QR factorization, with column pivoting, of A. That is, if
  !     A*P = Q*R, where P is a permutation matrix, Q has orthogonal
  !     columns, and R is an upper triangular matrix with diagonal
  !     elements of nonincreasing magnitude, then LMPAR expects
  !     the full upper triangle of R, the permutation matrix P,
  !     and the first N components of (Q TRANSPOSE)*B. On output
  !     LMPAR also provides an upper triangular matrix S such that
  !
  !            T   T                   T
  !           P *(A *A + PAR*D*D)*P = S *S .
  !
  !     S is employed within LMPAR and may be of separate interest.
  !
  !     Only a few iterations are generally needed for convergence
  !     of the algorithm. If, however, the limit of 10 iterations
  !     is reached, then the output PAR will contain the best
  !     value obtained so far.
  !
  !     The subroutine statement is
  !
  !       SUBROUTINE LMPAR(N,R,LDR,IPVT,DIAG,QTB,DELTA,PAR,X,SIGMA,
  !                        WA1,WA2)
  !
  !     where
  !
  !       N is a positive integer input variable set to the order of R.
  !
  !       R is an N by N array. On input the full upper triangle
  !         must contain the full upper triangle of the matrix R.
  !         On output the full upper triangle is unaltered, and the
  !         strict lower triangle contains the strict upper triangle
  !         (transposed) of the upper triangular matrix S.
  !
  !       LDR is a positive integer input variable not less than N
  !         which specifies the leading dimension of the array R.
  !
  !       IPVT is an integer input array of length N which defines the
  !         permutation matrix P such that A*P = Q*R. Column J of P
  !         is column IPVT(J) of the identity matrix.
  !
  !       DIAG is an input array of length N which must contain the
  !         diagonal elements of the matrix D.
  !
  !       QTB is an input array of length N which must contain the first
  !         N elements of the vector (Q TRANSPOSE)*B.
  !
  !       DELTA is a positive input variable which specifies an upper
  !         bound on the Euclidean norm of D*X.
  !
  !       PAR is a nonnegative variable. On input PAR contains an
  !         initial estimate of the Levenberg-Marquardt parameter.
  !         On output PAR contains the final estimate.
  !
  !       X is an output array of length N which contains the least
  !         squares solution of the system A*X = B, SQRT(PAR)*D*X = 0,
  !         for the output PAR.
  !
  !       SIGMA is an output array of length N which contains the
  !         diagonal elements of the upper triangular matrix S.
  !
  !       WA1 and WA2 are work arrays of length N.
  !
  !***
  ! **See also:**  SNLS1, SNLS1E
  !***
  ! **Routines called:**  ENORM, QRSOLV, R1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   800301  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   900328  Added TYPE section.  (WRB)
  USE service, ONLY : R1MACH
  INTEGER :: N, Ldr
  INTEGER :: Ipvt(N)
  REAL :: Delta, Par
  REAL :: R(Ldr,N), Diag(N), Qtb(N), X(N), Sigma(N), Wa1(N), Wa2(N)
  INTEGER :: i, iter, j, jm1, jp1, k, l, nsing
  REAL :: dxnorm, dwarf, fp, gnorm, parc, parl, paru, summ, temp
  REAL, PARAMETER :: p1 = 1.0E-1, p001 = 1.0E-3, zero = 0.0E0
  !* FIRST EXECUTABLE STATEMENT  LMPAR
  dwarf = R1MACH(1)
  !
  !     COMPUTE AND STORE IN X THE GAUSS-NEWTON DIRECTION. IF THE
  !     JACOBIAN IS RANK-DEFICIENT, OBTAIN A LEAST SQUARES SOLUTION.
  !
  nsing = N
  DO j = 1, N
    Wa1(j) = Qtb(j)
    IF ( R(j,j)==zero.AND.nsing==N ) nsing = j - 1
    IF ( nsing<N ) Wa1(j) = zero
  END DO
  IF ( nsing>=1 ) THEN
    DO k = 1, nsing
      j = nsing - k + 1
      Wa1(j) = Wa1(j)/R(j,j)
      temp = Wa1(j)
      jm1 = j - 1
      IF ( jm1>=1 ) THEN
        DO i = 1, jm1
          Wa1(i) = Wa1(i) - R(i,j)*temp
        END DO
      END IF
    END DO
  END IF
  DO j = 1, N
    l = Ipvt(j)
    X(l) = Wa1(j)
  END DO
  !
  !     INITIALIZE THE ITERATION COUNTER.
  !     EVALUATE THE FUNCTION AT THE ORIGIN, AND TEST
  !     FOR ACCEPTANCE OF THE GAUSS-NEWTON DIRECTION.
  !
  iter = 0
  DO j = 1, N
    Wa2(j) = Diag(j)*X(j)
  END DO
  dxnorm = ENORM(N,Wa2)
  fp = dxnorm - Delta
  IF ( fp<=p1*Delta ) THEN
    !
    !     TERMINATION.
    !
    IF ( iter==0 ) Par = zero
  ELSE
    !
    !     IF THE JACOBIAN IS NOT RANK DEFICIENT, THE NEWTON
    !     STEP PROVIDES A LOWER BOUND, PARL, FOR THE ZERO OF
    !     THE FUNCTION. OTHERWISE SET THIS BOUND TO ZERO.
    !
    parl = zero
    IF ( nsing>=N ) THEN
      DO j = 1, N
        l = Ipvt(j)
        Wa1(j) = Diag(l)*(Wa2(l)/dxnorm)
      END DO
      DO j = 1, N
        summ = zero
        jm1 = j - 1
        IF ( jm1>=1 ) THEN
          DO i = 1, jm1
            summ = summ + R(i,j)*Wa1(i)
          END DO
        END IF
        Wa1(j) = (Wa1(j)-summ)/R(j,j)
      END DO
      temp = ENORM(N,Wa1)
      parl = ((fp/Delta)/temp)/temp
    END IF
    !
    !     CALCULATE AN UPPER BOUND, PARU, FOR THE ZERO OF THE FUNCTION.
    !
    DO j = 1, N
      summ = zero
      DO i = 1, j
        summ = summ + R(i,j)*Qtb(i)
      END DO
      l = Ipvt(j)
      Wa1(j) = summ/Diag(l)
    END DO
    gnorm = ENORM(N,Wa1)
    paru = gnorm/Delta
    IF ( paru==zero ) paru = dwarf/MIN(Delta,p1)
    !
    !     IF THE INPUT PAR LIES OUTSIDE OF THE INTERVAL (PARL,PARU),
    !     SET PAR TO THE CLOSER ENDPOINT.
    !
    Par = MAX(Par,parl)
    Par = MIN(Par,paru)
    IF ( Par==zero ) Par = gnorm/dxnorm
    DO
      !
      !     BEGINNING OF AN ITERATION.
      !
      iter = iter + 1
      !
      !        EVALUATE THE FUNCTION AT THE CURRENT VALUE OF PAR.
      !
      IF ( Par==zero ) Par = MAX(dwarf,p001*paru)
      temp = SQRT(Par)
      DO j = 1, N
        Wa1(j) = temp*Diag(j)
      END DO
      CALL QRSOLV(N,R,Ldr,Ipvt,Wa1,Qtb,X,Sigma,Wa2)
      DO j = 1, N
        Wa2(j) = Diag(j)*X(j)
      END DO
      dxnorm = ENORM(N,Wa2)
      temp = fp
      fp = dxnorm - Delta
      !
      !        IF THE FUNCTION IS SMALL ENOUGH, ACCEPT THE CURRENT VALUE
      !        OF PAR. ALSO TEST FOR THE EXCEPTIONAL CASES WHERE PARL
      !        IS ZERO OR THE NUMBER OF ITERATIONS HAS REACHED 10.
      !
      IF ( ABS(fp)<=p1*Delta.OR.parl==zero.AND.fp<=temp.AND.temp<zero.OR.&
          iter==10 ) THEN
        IF ( iter==0 ) Par = zero
        EXIT
      ELSE
        !
        !        COMPUTE THE NEWTON CORRECTION.
        !
        DO j = 1, N
          l = Ipvt(j)
          Wa1(j) = Diag(l)*(Wa2(l)/dxnorm)
        END DO
        DO j = 1, N
          Wa1(j) = Wa1(j)/Sigma(j)
          temp = Wa1(j)
          jp1 = j + 1
          IF ( N>=jp1 ) THEN
            DO i = jp1, N
              Wa1(i) = Wa1(i) - R(i,j)*temp
            END DO
          END IF
        END DO
        temp = ENORM(N,Wa1)
        parc = ((fp/Delta)/temp)/temp
        !
        !        DEPENDING ON THE SIGN OF THE FUNCTION, UPDATE PARL OR PARU.
        !
        IF ( fp>zero ) parl = MAX(parl,Par)
        IF ( fp<zero ) paru = MIN(paru,Par)
        !
        !        COMPUTE AN IMPROVED ESTIMATE FOR PAR.
        !
        !
        !        END OF AN ITERATION.
        !
        Par = MAX(parl,Par+parc)
      END IF
    END DO
  END IF
  !
  !     LAST CARD OF SUBROUTINE LMPAR.
  !
END SUBROUTINE LMPAR
