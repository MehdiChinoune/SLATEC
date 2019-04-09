!** CMPTR3
SUBROUTINE CMPTR3(M,A,B,C,K,Y1,Y2,Y3,Tcos,D,W1,W2,W3)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to CMGNBN
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      COMPLEX (TRI3-S, CMPTR3-C)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     Subroutine to solve tridiagonal systems.
  !
  !***
  ! **See also:**  CMGNBN
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890206  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)

  INTEGER i, ip, K(4), k1, k2, k2k3k4, k3, k4, kint1, kint2, kint3, &
    l1, l2, l3, lint1, lint2, lint3, M, mm1, n
  COMPLEX A(*), B(*), C(*), Y1(*), Y2(*), Y3(*), Tcos(*), D(*), W1(*), W2(*), &
    W3(*), x, xx, z
  INTEGER k1p1, k2p1, k3p1, k4p1
  !
  !* FIRST EXECUTABLE STATEMENT  CMPTR3
  mm1 = M - 1
  k1 = K(1)
  k2 = K(2)
  k3 = K(3)
  k4 = K(4)
  k1p1 = k1 + 1
  k2p1 = k2 + 1
  k3p1 = k3 + 1
  k4p1 = k4 + 1
  k2k3k4 = k2 + k3 + k4
  IF ( k2k3k4/=0 ) THEN
    l1 = k1p1/k2p1
    l2 = k1p1/k3p1
    l3 = k1p1/k4p1
    lint1 = 1
    lint2 = 1
    lint3 = 1
    kint1 = k1
    kint2 = kint1 + k2
    kint3 = kint2 + k3
  END IF
  DO n = 1, k1
    x = Tcos(n)
    IF ( k2k3k4/=0 ) THEN
      IF ( n==l1 ) THEN
        DO i = 1, M
          W1(i) = Y1(i)
        END DO
      END IF
      IF ( n==l2 ) THEN
        DO i = 1, M
          W2(i) = Y2(i)
        END DO
      END IF
      IF ( n==l3 ) THEN
        DO i = 1, M
          W3(i) = Y3(i)
        END DO
      END IF
    END IF
    z = 1./(B(1)-x)
    D(1) = C(1)*z
    Y1(1) = Y1(1)*z
    Y2(1) = Y2(1)*z
    Y3(1) = Y3(1)*z
    DO i = 2, M
      z = 1./(B(i)-x-A(i)*D(i-1))
      D(i) = C(i)*z
      Y1(i) = (Y1(i)-A(i)*Y1(i-1))*z
      Y2(i) = (Y2(i)-A(i)*Y2(i-1))*z
      Y3(i) = (Y3(i)-A(i)*Y3(i-1))*z
    END DO
    DO ip = 1, mm1
      i = M - ip
      Y1(i) = Y1(i) - D(i)*Y1(i+1)
      Y2(i) = Y2(i) - D(i)*Y2(i+1)
      Y3(i) = Y3(i) - D(i)*Y3(i+1)
    END DO
    IF ( k2k3k4/=0 ) THEN
      IF ( n==l1 ) THEN
        i = lint1 + kint1
        xx = x - Tcos(i)
        DO i = 1, M
          Y1(i) = xx*Y1(i) + W1(i)
        END DO
        lint1 = lint1 + 1
        l1 = (lint1*k1p1)/k2p1
      END IF
      IF ( n==l2 ) THEN
        i = lint2 + kint2
        xx = x - Tcos(i)
        DO i = 1, M
          Y2(i) = xx*Y2(i) + W2(i)
        END DO
        lint2 = lint2 + 1
        l2 = (lint2*k1p1)/k3p1
      END IF
      IF ( n==l3 ) THEN
        i = lint3 + kint3
        xx = x - Tcos(i)
        DO i = 1, M
          Y3(i) = xx*Y3(i) + W3(i)
        END DO
        lint3 = lint3 + 1
        l3 = (lint3*k1p1)/k4p1
      END IF
    END IF
  END DO
END SUBROUTINE CMPTR3
