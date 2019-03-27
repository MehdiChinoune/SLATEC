!** TRIS4
SUBROUTINE TRIS4(N,A,B,C,D,U,Z)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to SEPX4
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (TRIS4-S)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     This subroutine solves for a non-zero eigenvector corresponding
  !     to the zero eigenvalue of the transpose of the rank
  !     deficient ONE matrix with subdiagonal A, diagonal B, and
  !     superdiagonal C, with A(1) in the (1,N) position, with
  !     C(N) in the (N,1) position, AND all other elements zero.
  !
  !***
  ! **See also:**  SEPX4
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  
  REAL A(*), an, B(*), bn, C(*), D(*), den, U(*), v, Z(*)
  INTEGER j, k, N, nm1, nm2
  !* FIRST EXECUTABLE STATEMENT  TRIS4
  bn = B(N)
  D(1) = A(2)/B(1)
  v = A(1)
  U(1) = C(N)/B(1)
  nm2 = N - 2
  DO j = 2, nm2
    den = B(j) - C(j-1)*D(j-1)
    D(j) = A(j+1)/den
    U(j) = -C(j-1)*U(j-1)/den
    bn = bn - v*U(j-1)
    v = -v*D(j-1)
  ENDDO
  den = B(N-1) - C(N-2)*D(N-2)
  D(N-1) = (A(N)-C(N-2)*U(N-2))/den
  an = C(N-1) - v*D(N-2)
  bn = bn - v*U(N-2)
  den = bn - an*D(N-1)
  !
  !     SET LAST COMPONENT EQUAL TO ONE
  !
  Z(N) = 1.0
  Z(N-1) = -D(N-1)
  nm1 = N - 1
  DO j = 2, nm1
    k = N - j
    Z(k) = -D(k)*Z(k+1) - U(k)*Z(N)
  ENDDO
END SUBROUTINE TRIS4
