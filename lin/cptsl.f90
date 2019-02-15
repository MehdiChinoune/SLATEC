!DECK CPTSL
SUBROUTINE CPTSL(N,D,E,B)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CPTSL
  !***PURPOSE  Solve a positive definite tridiagonal linear system.
  !***LIBRARY   SLATEC (LINPACK)
  !***CATEGORY  D2D2A
  !***TYPE      COMPLEX (SPTSL-S, DPTSL-D, CPTSL-C)
  !***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, POSITIVE DEFINITE, SOLVE,
  !             TRIDIAGONAL
  !***AUTHOR  Dongarra, J., (ANL)
  !***DESCRIPTION
  !
  !     CPTSL given a positive definite tridiagonal matrix and a right
  !     hand side will find the solution.
  !
  !     On Entry
  !
  !        N        INTEGER
  !                 is the order of the tridiagonal matrix.
  !
  !        D        COMPLEX(N)
  !                 is the diagonal of the tridiagonal matrix.
  !                 On output D is destroyed.
  !
  !        E        COMPLEX(N)
  !                 is the offdiagonal of the tridiagonal matrix.
  !                 E(1) through E(N-1) should contain the
  !                 offdiagonal.
  !
  !        B        COMPLEX(N)
  !                 is the right hand side vector.
  !
  !     On Return
  !
  !        B        contains the solution.
  !
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890505  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  CPTSL
  INTEGER N
  COMPLEX D(*), E(*), B(*)
  !
  INTEGER k, kbm1, ke, kf, kp1, nm1, nm1d2
  COMPLEX t1, t2
  !
  !     CHECK FOR 1 X 1 CASE
  !
  !***FIRST EXECUTABLE STATEMENT  CPTSL
  IF ( N/=1 ) THEN
    nm1 = N - 1
    nm1d2 = nm1/2
    IF ( N/=2 ) THEN
      kbm1 = N - 1
      !
      !           ZERO TOP HALF OF SUBDIAGONAL AND BOTTOM HALF OF
      !           SUPERDIAGONAL
      !
      DO k = 1, nm1d2
        t1 = CONJG(E(k))/D(k)
        D(k+1) = D(k+1) - t1*E(k)
        B(k+1) = B(k+1) - t1*B(k)
        t2 = E(kbm1)/D(kbm1+1)
        D(kbm1) = D(kbm1) - t2*CONJG(E(kbm1))
        B(kbm1) = B(kbm1) - t2*B(kbm1+1)
        kbm1 = kbm1 - 1
      ENDDO
    ENDIF
    kp1 = nm1d2 + 1
    !
    !        CLEAN UP FOR POSSIBLE 2 X 2 BLOCK AT CENTER
    !
    IF ( MOD(N,2)==0 ) THEN
      t1 = CONJG(E(kp1))/D(kp1)
      D(kp1+1) = D(kp1+1) - t1*E(kp1)
      B(kp1+1) = B(kp1+1) - t1*B(kp1)
      kp1 = kp1 + 1
    ENDIF
    !
    !        BACK SOLVE STARTING AT THE CENTER, GOING TOWARDS THE TOP
    !        AND BOTTOM
    !
    B(kp1) = B(kp1)/D(kp1)
    IF ( N/=2 ) THEN
      k = kp1 - 1
      ke = kp1 + nm1d2 - 1
      DO kf = kp1, ke
        B(k) = (B(k)-E(k)*B(k+1))/D(k)
        B(kf+1) = (B(kf+1)-CONJG(E(kf))*B(kf))/D(kf+1)
        k = k - 1
      ENDDO
    ENDIF
    IF ( MOD(N,2)==0 ) B(1) = (B(1)-E(1)*B(2))/D(1)
  ELSE
    B(1) = B(1)/D(1)
  ENDIF
END SUBROUTINE CPTSL
