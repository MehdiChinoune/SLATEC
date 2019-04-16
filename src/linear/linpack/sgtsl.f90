!** SGTSL
SUBROUTINE SGTSL(N,C,D,E,B,Info)
  !>
  !***
  !  Solve a tridiagonal linear system.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D2A2A
  !***
  ! **Type:**      SINGLE PRECISION (SGTSL-S, DGTSL-D, CGTSL-C)
  !***
  ! **Keywords:**  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE, TRIDIAGONAL
  !***
  ! **Author:**  Dongarra, J., (ANL)
  !***
  ! **Description:**
  !
  !     SGTSL given a general tridiagonal matrix and a right hand
  !     side will find the solution.
  !
  !     On Entry
  !
  !        N       INTEGER
  !                is the order of the tridiagonal matrix.
  !
  !        C       REAL(N)
  !                is the subdiagonal of the tridiagonal matrix.
  !                C(2) through C(N) should contain the subdiagonal.
  !                On output, C is destroyed.
  !
  !        D       REAL(N)
  !                is the diagonal of the tridiagonal matrix.
  !                On output, D is destroyed.
  !
  !        E       REAL(N)
  !                is the superdiagonal of the tridiagonal matrix.
  !                E(1) through E(N-1) should contain the superdiagonal.
  !                On output, E is destroyed.
  !
  !        B       REAL(N)
  !                is the right hand side vector.
  !
  !     On Return
  !
  !        B       is the solution vector.
  !
  !        INFO    INTEGER
  !                = 0 normal value.
  !                = K if the K-th element of the diagonal becomes
  !                    exactly zero.  The subroutine returns when
  !                    this is detected.
  !
  !***
  ! **References:**  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER N, Info
  REAL C(*), D(*), E(*), B(*)
  !
  INTEGER k, kb, kp1, nm1, nm2
  REAL t
  !* FIRST EXECUTABLE STATEMENT  SGTSL
  Info = 0
  C(1) = D(1)
  nm1 = N - 1
  IF ( nm1>=1 ) THEN
    D(1) = E(1)
    E(1) = 0.0E0
    E(N) = 0.0E0
    !
    DO k = 1, nm1
      kp1 = k + 1
      !
      !              FIND THE LARGEST OF THE TWO ROWS
      !
      IF ( ABS(C(kp1))>=ABS(C(k)) ) THEN
        !
        !                 INTERCHANGE ROW
        !
        t = C(kp1)
        C(kp1) = C(k)
        C(k) = t
        t = D(kp1)
        D(kp1) = D(k)
        D(k) = t
        t = E(kp1)
        E(kp1) = E(k)
        E(k) = t
        t = B(kp1)
        B(kp1) = B(k)
        B(k) = t
      END IF
      !
      !              ZERO ELEMENTS
      !
      IF ( C(k)/=0.0E0 ) THEN
        t = -C(kp1)/C(k)
        C(kp1) = D(kp1) + t*D(k)
        D(kp1) = E(kp1) + t*E(k)
        E(kp1) = 0.0E0
        B(kp1) = B(kp1) + t*B(k)
      ELSE
        Info = k
        RETURN
      END IF
    END DO
  END IF
  IF ( C(N)/=0.0E0 ) THEN
    !
    !           BACK SOLVE
    !
    nm2 = N - 2
    B(N) = B(N)/C(N)
    IF ( N/=1 ) THEN
      B(nm1) = (B(nm1)-D(nm1)*B(N))/C(nm1)
      IF ( nm2>=1 ) THEN
        DO kb = 1, nm2
          k = nm2 - kb + 1
          B(k) = (B(k)-D(k)*B(k+1)-E(k)*B(k+2))/C(k)
        END DO
      END IF
    END IF
  ELSE
    Info = N
  END IF
  !
  RETURN
END SUBROUTINE SGTSL
