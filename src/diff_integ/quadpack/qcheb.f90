!** QCHEB
SUBROUTINE QCHEB(X,Fval,Cheb12,Cheb24)
  !>
  !  This routine computes the CHEBYSHEV series expansion
  !            of degrees 12 and 24 of a function using A
  !            FAST FOURIER TRANSFORM METHOD
  !            F(X) = SUM(K=1,..,13) (CHEB12(K)*T(K-1,X)),
  !            F(X) = SUM(K=1,..,25) (CHEB24(K)*T(K-1,X)),
  !            Where T(K,X) is the CHEBYSHEV POLYNOMIAL OF DEGREE K.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (QCHEB-S, DQCHEB-D)
  !***
  ! **Keywords:**  CHEBYSHEV SERIES EXPANSION, FAST FOURIER TRANSFORM
  !***
  ! **Author:**  Piessens, Robert
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !           de Doncker, Elise
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !***
  ! **Description:**
  !
  !        Chebyshev Series Expansion
  !        Standard Fortran Subroutine
  !        Real version
  !
  !        PARAMETERS
  !          ON ENTRY
  !           X      - Real
  !                    Vector of dimension 11 containing the
  !                    Values COS(K*PI/24), K = 1, ..., 11
  !
  !           FVAL   - Real
  !                    Vector of dimension 25 containing the
  !                    function values at the points
  !                    (B+A+(B-A)*COS(K*PI/24))/2, K = 0, ...,24,
  !                    where (A,B) is the approximation interval.
  !                    FVAL(1) and FVAL(25) are divided by two
  !                    (these values are destroyed at output).
  !
  !          ON RETURN
  !           CHEB12 - Real
  !                    Vector of dimension 13 containing the
  !                    CHEBYSHEV coefficients for degree 12
  !
  !           CHEB24 - Real
  !                    Vector of dimension 25 containing the
  !                    CHEBYSHEV Coefficients for degree 24
  !
  !***
  ! **See also:**  QC25C, QC25F, QC25S
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   810101  DATE WRITTEN
  !   830518  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)

  !
  REAL(SP) alam, alam1, alam2, Cheb12(13), Cheb24(25), Fval(25), part1, part2, &
    part3, v(12), X(11)
  INTEGER i, j
  !
  !* FIRST EXECUTABLE STATEMENT  QCHEB
  DO i = 1, 12
    j = 26 - i
    v(i) = Fval(i) - Fval(j)
    Fval(i) = Fval(i) + Fval(j)
  END DO
  alam1 = v(1) - v(9)
  alam2 = X(6)*(v(3)-v(7)-v(11))
  Cheb12(4) = alam1 + alam2
  Cheb12(10) = alam1 - alam2
  alam1 = v(2) - v(8) - v(10)
  alam2 = v(4) - v(6) - v(12)
  alam = X(3)*alam1 + X(9)*alam2
  Cheb24(4) = Cheb12(4) + alam
  Cheb24(22) = Cheb12(4) - alam
  alam = X(9)*alam1 - X(3)*alam2
  Cheb24(10) = Cheb12(10) + alam
  Cheb24(16) = Cheb12(10) - alam
  part1 = X(4)*v(5)
  part2 = X(8)*v(9)
  part3 = X(6)*v(7)
  alam1 = v(1) + part1 + part2
  alam2 = X(2)*v(3) + part3 + X(10)*v(11)
  Cheb12(2) = alam1 + alam2
  Cheb12(12) = alam1 - alam2
  alam = X(1)*v(2) + X(3)*v(4) + X(5)*v(6) + X(7)*v(8) + X(9)*v(10) + X(11)*v(12)
  Cheb24(2) = Cheb12(2) + alam
  Cheb24(24) = Cheb12(2) - alam
  alam = X(11)*v(2) - X(9)*v(4) + X(7)*v(6) - X(5)*v(8) + X(3)*v(10) - X(1)*v(12)
  Cheb24(12) = Cheb12(12) + alam
  Cheb24(14) = Cheb12(12) - alam
  alam1 = v(1) - part1 + part2
  alam2 = X(10)*v(3) - part3 + X(2)*v(11)
  Cheb12(6) = alam1 + alam2
  Cheb12(8) = alam1 - alam2
  alam = X(5)*v(2) - X(9)*v(4) - X(1)*v(6) - X(11)*v(8) + X(3)*v(10) + X(7)*v(12)
  Cheb24(6) = Cheb12(6) + alam
  Cheb24(20) = Cheb12(6) - alam
  alam = X(7)*v(2) - X(3)*v(4) - X(11)*v(6) + X(1)*v(8) - X(9)*v(10) - X(5)*v(12)
  Cheb24(8) = Cheb12(8) + alam
  Cheb24(18) = Cheb12(8) - alam
  DO i = 1, 6
    j = 14 - i
    v(i) = Fval(i) - Fval(j)
    Fval(i) = Fval(i) + Fval(j)
  END DO
  alam1 = v(1) + X(8)*v(5)
  alam2 = X(4)*v(3)
  Cheb12(3) = alam1 + alam2
  Cheb12(11) = alam1 - alam2
  Cheb12(7) = v(1) - v(5)
  alam = X(2)*v(2) + X(6)*v(4) + X(10)*v(6)
  Cheb24(3) = Cheb12(3) + alam
  Cheb24(23) = Cheb12(3) - alam
  alam = X(6)*(v(2)-v(4)-v(6))
  Cheb24(7) = Cheb12(7) + alam
  Cheb24(19) = Cheb12(7) - alam
  alam = X(10)*v(2) - X(6)*v(4) + X(2)*v(6)
  Cheb24(11) = Cheb12(11) + alam
  Cheb24(15) = Cheb12(11) - alam
  DO i = 1, 3
    j = 8 - i
    v(i) = Fval(i) - Fval(j)
    Fval(i) = Fval(i) + Fval(j)
  END DO
  Cheb12(5) = v(1) + X(8)*v(3)
  Cheb12(9) = Fval(1) - X(8)*Fval(3)
  alam = X(4)*v(2)
  Cheb24(5) = Cheb12(5) + alam
  Cheb24(21) = Cheb12(5) - alam
  alam = X(8)*Fval(2) - Fval(4)
  Cheb24(9) = Cheb12(9) + alam
  Cheb24(17) = Cheb12(9) - alam
  Cheb12(1) = Fval(1) + Fval(3)
  alam = Fval(2) + Fval(4)
  Cheb24(1) = Cheb12(1) + alam
  Cheb24(25) = Cheb12(1) - alam
  Cheb12(13) = v(1) - v(3)
  Cheb24(13) = Cheb12(13)
  alam = 0.1E+01/0.6E+01
  DO i = 2, 12
    Cheb12(i) = Cheb12(i)*alam
  END DO
  alam = 0.5E+00*alam
  Cheb12(1) = Cheb12(1)*alam
  Cheb12(13) = Cheb12(13)*alam
  DO i = 2, 24
    Cheb24(i) = Cheb24(i)*alam
  END DO
  Cheb24(1) = 0.5E+00*alam*Cheb24(1)
  Cheb24(25) = 0.5E+00*alam*Cheb24(25)
END SUBROUTINE QCHEB
