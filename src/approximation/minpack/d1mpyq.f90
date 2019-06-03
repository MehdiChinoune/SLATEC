!** D1MPYQ
SUBROUTINE D1MPYQ(M,N,A,Lda,V,W)
  !>
  !  Subsidiary to DNSQ and DNSQE
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (R1MPYQ-S, D1MPYQ-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     Given an M by N matrix A, this subroutine computes A*Q where
  !     Q is the product of 2*(N - 1) transformations
  !
  !           GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
  !
  !     and GV(I), GW(I) are Givens rotations in the (I,N) plane which
  !     eliminate elements in the I-th and N-th planes, respectively.
  !     Q itself is not given, rather the information to recover the
  !     GV, GW rotations is supplied.
  !
  !     The SUBROUTINE statement is
  !
  !       SUBROUTINE D1MPYQ(M,N,A,LDA,V,W)
  !
  !     where
  !
  !       M is a positive integer input variable set to the number
  !         of rows of A.
  !
  !       N IS a positive integer input variable set to the number
  !         of columns of A.
  !
  !       A is an M by N array. On input A must contain the matrix
  !         to be postmultiplied by the orthogonal matrix Q
  !         described above. On output A*Q has replaced A.
  !
  !       LDA is a positive integer input variable not less than M
  !         which specifies the leading dimension of the array A.
  !
  !       V is an input array of length N. V(I) must contain the
  !         information necessary to recover the Givens rotation GV(I)
  !         described above.
  !
  !       W is an input array of length N. W(I) must contain the
  !         information necessary to recover the Givens rotation GW(I)
  !         described above.
  !
  !***
  ! **See also:**  DNSQ, DNSQE
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

  INTEGER Lda, M, N
  REAL(DP) :: A(Lda,N), V(N), W(N)
  INTEGER :: i, j, nm1, nmj
  REAL(DP) :: coss, sinn, temp
  REAL(DP), PARAMETER :: one = 1.0D0
  !
  !     APPLY THE FIRST SET OF GIVENS ROTATIONS TO A.
  !
  !* FIRST EXECUTABLE STATEMENT  D1MPYQ
  nm1 = N - 1
  IF ( nm1>=1 ) THEN
    DO nmj = 1, nm1
      j = N - nmj
      IF ( ABS(V(j))>one ) coss = one/V(j)
      IF ( ABS(V(j))>one ) sinn = SQRT(one-coss**2)
      IF ( ABS(V(j))<=one ) sinn = V(j)
      IF ( ABS(V(j))<=one ) coss = SQRT(one-sinn**2)
      DO i = 1, M
        temp = coss*A(i,j) - sinn*A(i,N)
        A(i,N) = sinn*A(i,j) + coss*A(i,N)
        A(i,j) = temp
      END DO
    END DO
    !
    !     APPLY THE SECOND SET OF GIVENS ROTATIONS TO A.
    !
    DO j = 1, nm1
      IF ( ABS(W(j))>one ) coss = one/W(j)
      IF ( ABS(W(j))>one ) sinn = SQRT(one-coss**2)
      IF ( ABS(W(j))<=one ) sinn = W(j)
      IF ( ABS(W(j))<=one ) coss = SQRT(one-sinn**2)
      DO i = 1, M
        temp = coss*A(i,j) + sinn*A(i,N)
        A(i,N) = -sinn*A(i,j) + coss*A(i,N)
        A(i,j) = temp
      END DO
    END DO
  END IF
  !
  !     LAST CARD OF SUBROUTINE D1MPYQ.
  !
END SUBROUTINE D1MPYQ
