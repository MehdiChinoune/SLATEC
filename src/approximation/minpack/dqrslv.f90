!** DQRSLV
SUBROUTINE DQRSLV(N,R,Ldr,Ipvt,Diag,Qtb,X,Sigma,Wa)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to DNLS1 and DNLS1E
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (QRSOLV-S, DQRSLV-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !  **** Double Precision version of QRSOLV ****
  !
  !     Given an M by N matrix A, an N by N diagonal matrix D,
  !     and an M-vector B, the problem is to determine an X which
  !     solves the system
  !
  !           A*X = B,     D*X = 0 ,
  !
  !     in the least squares sense.
  !
  !     This subroutine completes the solution of the problem
  !     if it is provided with the necessary information from the
  !     QR factorization, with column pivoting, of A. That is, if
  !     A*P = Q*R, where P is a permutation matrix, Q has orthogonal
  !     columns, and R is an upper triangular matrix with diagonal
  !     elements of nonincreasing magnitude, then DQRSLV expects
  !     the full upper triangle of R, the permutation matrix P,
  !     and the first N components of (Q TRANSPOSE)*B. The system
  !     A*X = B, D*X = 0, is then equivalent to
  !
  !                  T       T
  !           R*Z = Q *B,  P *D*P*Z = 0 ,
  !
  !     where X = P*Z. If this system does not have full rank,
  !     then a least squares solution is obtained. On output DQRSLV
  !     also provides an upper triangular matrix S such that
  !
  !            T   T               T
  !           P *(A *A + D*D)*P = S *S .
  !
  !     S is computed within DQRSLV and may be of separate interest.
  !
  !     The subroutine statement is
  !
  !       SUBROUTINE DQRSLV(N,R,LDR,IPVT,DIAG,QTB,X,SIGMA,WA)
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
  !       X is an output array of length N which contains the least
  !         squares solution of the system A*X = B, D*X = 0.
  !
  !       SIGMA is an output array of length N which contains the
  !         diagonal elements of the upper triangular matrix S.
  !
  !       WA is a work array of length N.
  !
  !***
  ! **See also:**  DNLS1, DNLS1E
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   800301  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   900328  Added TYPE section.  (WRB)

  INTEGER N, Ldr
  INTEGER Ipvt(*)
  REAL(8) :: R(Ldr,*), Diag(*), Qtb(*), X(*), Sigma(*), Wa(*)
  INTEGER i, j, jp1, k, kp1, l, nsing
  REAL(8) :: coss, cotan, qtbpj, sinn, summ, tann, temp
  REAL(8), PARAMETER :: p5 = 5.0D-1, p25 = 2.5D-1, zero = 0.0D0
  !* FIRST EXECUTABLE STATEMENT  DQRSLV
  DO j = 1, N
    DO i = j, N
      R(i,j) = R(j,i)
    END DO
    X(j) = R(j,j)
    Wa(j) = Qtb(j)
  END DO
  !
  !     ELIMINATE THE DIAGONAL MATRIX D USING A GIVENS ROTATION.
  !
  DO j = 1, N
    !
    !        PREPARE THE ROW OF D TO BE ELIMINATED, LOCATING THE
    !        DIAGONAL ELEMENT USING P FROM THE QR FACTORIZATION.
    !
    l = Ipvt(j)
    IF ( Diag(l)/=zero ) THEN
      DO k = j, N
        Sigma(k) = zero
      END DO
      Sigma(j) = Diag(l)
      !
      !        THE TRANSFORMATIONS TO ELIMINATE THE ROW OF D
      !        MODIFY ONLY A SINGLE ELEMENT OF (Q TRANSPOSE)*B
      !        BEYOND THE FIRST N, WHICH IS INITIALLY ZERO.
      !
      qtbpj = zero
      DO k = j, N
        !
        !           DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
        !           APPROPRIATE ELEMENT IN THE CURRENT ROW OF D.
        !
        IF ( Sigma(k)/=zero ) THEN
          IF ( ABS(R(k,k))>=ABS(Sigma(k)) ) THEN
            tann = Sigma(k)/R(k,k)
            coss = p5/SQRT(p25+p25*tann**2)
            sinn = coss*tann
          ELSE
            cotan = R(k,k)/Sigma(k)
            sinn = p5/SQRT(p25+p25*cotan**2)
            coss = sinn*cotan
          END IF
          !
          !           COMPUTE THE MODIFIED DIAGONAL ELEMENT OF R AND
          !           THE MODIFIED ELEMENT OF ((Q TRANSPOSE)*B,0).
          !
          R(k,k) = coss*R(k,k) + sinn*Sigma(k)
          temp = coss*Wa(k) + sinn*qtbpj
          qtbpj = -sinn*Wa(k) + coss*qtbpj
          Wa(k) = temp
          !
          !           ACCUMULATE THE TRANSFORMATION IN THE ROW OF S.
          !
          kp1 = k + 1
          IF ( N>=kp1 ) THEN
            DO i = kp1, N
              temp = coss*R(i,k) + sinn*Sigma(i)
              Sigma(i) = -sinn*R(i,k) + coss*Sigma(i)
              R(i,k) = temp
            END DO
          END IF
        END IF
      END DO
    END IF
    !
    !        STORE THE DIAGONAL ELEMENT OF S AND RESTORE
    !        THE CORRESPONDING DIAGONAL ELEMENT OF R.
    !
    Sigma(j) = R(j,j)
    R(j,j) = X(j)
  END DO
  !
  !     SOLVE THE TRIANGULAR SYSTEM FOR Z. IF THE SYSTEM IS
  !     SINGULAR, THEN OBTAIN A LEAST SQUARES SOLUTION.
  !
  nsing = N
  DO j = 1, N
    IF ( Sigma(j)==zero.AND.nsing==N ) nsing = j - 1
    IF ( nsing<N ) Wa(j) = zero
  END DO
  IF ( nsing>=1 ) THEN
    DO k = 1, nsing
      j = nsing - k + 1
      summ = zero
      jp1 = j + 1
      IF ( nsing>=jp1 ) THEN
        DO i = jp1, nsing
          summ = summ + R(i,j)*Wa(i)
        END DO
      END IF
      Wa(j) = (Wa(j)-summ)/Sigma(j)
    END DO
  END IF
  !
  !     PERMUTE THE COMPONENTS OF Z BACK TO COMPONENTS OF X.
  !
  DO j = 1, N
    l = Ipvt(j)
    X(l) = Wa(j)
  END DO
  !
  !     LAST CARD OF SUBROUTINE DQRSLV.
  !
END SUBROUTINE DQRSLV
