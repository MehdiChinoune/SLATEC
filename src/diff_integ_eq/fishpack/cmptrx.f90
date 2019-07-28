!** CMPTRX
PURE SUBROUTINE CMPTRX(Idegbr,Idegcr,M,A,B,C,Y,Tcos,D,W)
  !> Subsidiary to CMGNBN
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      COMPLEX (TRIX-S, CMPTRX-C)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     Subroutine to solve a system of linear equations where the
  !     coefficient matrix is a rational function in the matrix given by
  !     tridiagonal  ( . . ., A(I), B(I), C(I), . . . ).
  !
  !***
  ! **See also:**  CMGNBN
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)

  INTEGER, INTENT(IN) :: Idegbr, Idegcr, M
  COMPLEX(SP), INTENT(IN) :: A(M), B(M), C(M), Tcos(:)
  COMPLEX(SP), INTENT(INOUT) :: Y(:)
  COMPLEX(SP), INTENT(OUT) :: D(M), W(M)
  !
  INTEGER :: i, ip, k, l, lint, mm1, kb, kc
  COMPLEX(SP) :: x, xx, z
  !* FIRST EXECUTABLE STATEMENT  CMPTRX
  mm1 = M - 1
  kb = Idegbr + 1
  kc = Idegcr + 1
  l = kb/kc
  lint = 1
  DO k = 1, Idegbr
    x = Tcos(k)
    IF( k==l ) THEN
      i = Idegbr + lint
      xx = x - Tcos(i)
      DO i = 1, M
        W(i) = Y(i)
        Y(i) = xx*Y(i)
      END DO
    END IF
    z = 1._SP/(B(1)-x)
    D(1) = C(1)*z
    Y(1) = Y(1)*z
    DO i = 2, mm1
      z = 1._SP/(B(i)-x-A(i)*D(i-1))
      D(i) = C(i)*z
      Y(i) = (Y(i)-A(i)*Y(i-1))*z
    END DO
    z = B(M) - x - A(M)*D(mm1)
    IF( ABS(z)/=0. ) THEN
      Y(M) = (Y(M)-A(M)*Y(mm1))/z
    ELSE
      Y(M) = (0._SP,0._SP)
    END IF
    DO ip = 1, mm1
      i = M - ip
      Y(i) = Y(i) - D(i)*Y(i+1)
    END DO
    IF( k==l ) THEN
      DO i = 1, M
        Y(i) = Y(i) + W(i)
      END DO
      lint = lint + 1
      l = (lint*kb)/kc
    END IF
  END DO
  !
END SUBROUTINE CMPTRX